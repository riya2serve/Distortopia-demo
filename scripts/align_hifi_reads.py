# scripts/align_hifi_reads_cross_app.py
import shutil
import subprocess
from pathlib import Path
from typing import Tuple, Optional
import pandas as pd
import streamlit as st

# ------------------ helpers ------------------
def check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def run_pipe_to_file(p1, p2, outfile):
    """Run p1 | p2 -o outfile  (minimap2 -> samtools sort -> BAM)."""
    st.code("$ " + " ".join(p1) + " | " + " ".join(p2 + ["-o", outfile]))
    p_minimap = subprocess.Popen(p1, stdout=subprocess.PIPE)
    p_sort = subprocess.Popen(p2 + ["-o", outfile], stdin=p_minimap.stdout)
    p_minimap.stdout.close()
    rc = p_sort.wait()
    if rc != 0:
        raise RuntimeError("samtools sort failed")

def run_cmd_text(cmd) -> str:
    """Run a command and return stdout; on failure, surface stderr (incl. pipes)."""
    st.code("$ " + " ".join(cmd))
    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return out.stdout
    except subprocess.CalledProcessError as e:
        msg = []
        if e.stdout:
            msg.append("STDOUT:\n" + e.stdout.strip())
        if e.stderr:
            msg.append("STDERR:\n" + e.stderr.strip())
        st.error("Command failed:\n" + ("\n\n".join(msg) or "<no output>"))
        raise

def parse_flagstat(txt: str) -> dict:
    d = {"total": None, "mapped": None, "primary_mapped": None}
    for line in txt.splitlines():
        ls = line.strip()
        if " in total" in ls and d["total"] is None:
            d["total"] = int(ls.split()[0])
        elif " mapped (" in ls and "primary" not in ls and d["mapped"] is None:
            d["mapped"] = int(ls.split()[0])
        elif "primary mapped" in ls and d["primary_mapped"] is None:
            d["primary_mapped"] = int(ls.split()[0])
    return d

def auto_out_bam(reads_path: Path, ref_path: Path, outdir: Path) -> Path:
    """Generate BAM output name for a reads–reference pair."""
    rname = reads_path.name
    if rname.endswith(".fq.gz"):
        r = rname[:-6]
    elif rname.endswith(".fastq.gz"):
        r = rname[:-9]
    else:
        r = reads_path.stem
    ref = ref_path.stem
    return outdir / f"{r}__to__{ref}.bam"

def auto_out_vcfs(bam_path: Path) -> Tuple[Path, Path, Path]:
    """VCF paths go to repo's data/ folder; ensure it exists."""
    data_dir = Path("data")
    data_dir.mkdir(parents=True, exist_ok=True)
    base = bam_path.stem
    vcf_gz      = data_dir / f"{base}.vcf.gz"
    norm_vcf_gz = data_dir / f"{base}.norm.vcf.gz"
    filt_vcf_gz = data_dir / f"{base}.filtered.vcf.gz"
    return vcf_gz, norm_vcf_gz, filt_vcf_gz

def ensure_fasta_index(ref_fa: Path):
    """Ensure reference FASTA has .fai index."""
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if not fai.exists():
        run_cmd_text(["samtools", "faidx", str(ref_fa)])

def ensure_bam_index(bam_p: Path):
    """Ensure BAM is indexed; quickcheck the BAM."""
    bai = Path(str(bam_p) + ".bai")
    if not bai.exists():
        run_cmd_text(["samtools", "index", str(bam_p)])
    # prints nothing if OK; raises if not
    run_cmd_text(["bash", "-lc", f"set -e; samtools quickcheck -v {bam_p} >/dev/null || exit 1"])

def vcf_count_records(vcf_path: Path) -> int:
    """Count variant records quickly."""
    cmd = ["bash", "-lc", f"bcftools view -H {vcf_path} | wc -l"]
    out = run_cmd_text(cmd)
    try:
        return int(out.strip())
    except Exception:
        return 0

def bcftools_call_and_filter(
    ref_fa: Path,
    bam_p: Path,
    mapq_min: int,
    baseq_min: int,
    normalize: bool,
    ploidy_mode: str,   # "diploid" or "haploid"
    qual_min: int,
    require_pass: bool,
) -> Tuple[Path, Optional[Path], Path, int, int]:
    """
    Call variants and produce three outputs into data/:
      - *.vcf.gz (all calls)
      - *.norm.vcf.gz (optional)
      - *.filtered.vcf.gz (SNPs only, biallelic, QUAL>qual_min, het if diploid, optional PASS)
    Returns: vcf_gz, norm_vcf_gz_or_None, filt_vcf_gz, n_all, n_filtered
    """
    if not check_tool("bcftools"):
        st.error("`bcftools` not found. Try: `conda install -c bioconda bcftools`")
        st.stop()

    ensure_fasta_index(ref_fa)
    ensure_bam_index(bam_p)
    vcf_gz, norm_vcf_gz, filt_vcf_gz = auto_out_vcfs(bam_p)

    # 1) mpileup -> call
    with st.spinner(f"bcftools call → {vcf_gz.name}"):
        ploidy_flag = "" if ploidy_mode == "diploid" else "--ploidy 1"
        cmd = [
            "bash", "-lc",
            (
                "set -o pipefail; "
                f"bcftools mpileup -f {ref_fa} -q {mapq_min} -Q {baseq_min} "
                f"-Ou -a AD,DP,SP,MQ {bam_p} "
                f"| bcftools call {ploidy_flag} -mv -Oz -o {vcf_gz} "
                f"&& bcftools index -t {vcf_gz}"
            )
        ]
        run_cmd_text(cmd)

    n_all = vcf_count_records(vcf_gz)

    # 2) normalize (optional)
    src_vcf = vcf_gz
    if normalize:
        with st.spinner(f"bcftools norm → {norm_vcf_gz.name}"):
            run_cmd_text(["bcftools", "norm", "-f", str(ref_fa), "-Oz", "-o", str(norm_vcf_gz), str(vcf_gz)])
            run_cmd_text(["bcftools", "index", "-t", str(norm_vcf_gz)])
        src_vcf = norm_vcf_gz

    # 3) strict filtered VCF
    with st.spinner(f"Filter HQ biallelic het SNPs → {filt_vcf_gz.name}"):
        clauses = [f"QUAL>{qual_min}"]
        if require_pass:
            clauses.append('FILTER="PASS"')
        if ploidy_mode == "diploid":
            clauses.append('GT="het"')  # haploid cannot be het
        expr = " && ".join(clauses)

        run_cmd_text([
            "bcftools", "view",
            "-v", "snps",       # exclude INDELs & MNPs
            "-m2", "-M2",       # biallelic only
            "-i", expr,         # QUAL / PASS / GT
            "-Oz", "-o", str(filt_vcf_gz), str(src_vcf)
        ])
        run_cmd_text(["bcftools", "index", "-t", str(filt_vcf_gz)])

    n_filtered = vcf_count_records(filt_vcf_gz)
    return vcf_gz, (norm_vcf_gz if normalize else None), filt_vcf_gz, n_all, n_filtered

# ------------------ UI ------------------
st.title("Distortopia: Cross-map HiFi reads + bcftools calling → data/")

st.markdown("""
Align simulated HiFi reads to two references (2×2 matrix), then call variants and
write VCFs into the repo’s **data/** folder. Filtered VCFs keep **biallelic SNPs only** with **QUAL>30**
and **heterozygous GT** (in diploid mode).
""")

col1, col2 = st.columns(2)
with col1:
    ref1 = st.text_input("Reference FASTA #1", value="raw_data/A_thaliana.fna")
    reads1 = st.text_input("Reads FASTQ #1", value="long_reads/Athaliana_gametes_hifi.fq.gz")
with col2:
    ref2 = st.text_input("Reference FASTA #2", value="raw_data/A_lyrata.fna")
    reads2 = st.text_input("Reads FASTQ #2", value="long_reads/Alyrata_gametes_hifi.fq.gz")

threads = st.number_input("Threads (minimap2)", min_value=1, max_value=64, value=8)
outdir = Path(st.text_input("BAM output directory", value="alignments"))
index_bam = st.checkbox("Index BAMs (samtools index)", value=True)
show_flagstat = st.checkbox("Show per-pair flagstat", value=False)

st.markdown("---")
st.subheader("bcftools settings")
call_variants = st.checkbox("Call variants for each BAM", value=True)

c1, c2, c3 = st.columns(3)
with c1:
    ploidy_mode = st.selectbox("Genotype ploidy", ["diploid", "haploid"], index=0)
with c2:
    normalize_vcf = st.checkbox("Normalize VCF (bcftools norm)", value=True)
with c3:
    require_pass = st.checkbox('Require FILTER="PASS"', value=True)

c4, c5, c6 = st.columns(3)
with c4:
    mapq_min = st.number_input("Read MAPQ ≥", 0, 60, 30)
with c5:
    baseq_min = st.number_input("Base QUAL ≥", 0, 60, 30)
with c6:
    qual_min = st.number_input("VCF QUAL >", 0, 500, 30)

# ------------------ main ------------------
if st.button("Run cross-alignments", type="primary"):
    # tool checks
    if not check_tool("minimap2") or not check_tool("samtools"):
        st.error("Please install minimap2 and samtools (`conda install -c bioconda minimap2 samtools`)")
        st.stop()
    if call_variants and not check_tool("bcftools"):
        st.error("Please install bcftools (`conda install -c bioconda bcftools`)")
        st.stop()

    # paths
    ref_paths = [Path(ref1), Path(ref2)]
    read_paths = [Path(reads1), Path(reads2)]
    missing = [str(p) for p in (ref_paths + read_paths) if not p.exists()]
    if missing:
        st.error("Missing file(s):\n" + "\n".join(missing))
        st.stop()

    outdir.mkdir(parents=True, exist_ok=True)
    jobs = [(r, ref, auto_out_bam(r, ref, outdir)) for r in read_paths for ref in ref_paths]
    results = []

    for reads_p, ref_p, bam_p in jobs:
        st.subheader(f"Map {reads_p.name} → {ref_p.name}")
        with st.spinner(f"minimap2 map-hifi → {bam_p.name}"):
            p1 = ["minimap2", "-t", str(int(threads)), "-ax", "map-hifi", str(ref_p), str(reads_p)]
            p2 = ["samtools", "sort"]
            run_pipe_to_file(p1, p2, str(bam_p))

        if index_bam:
            with st.spinner(f"Index {bam_p.name}"):
                run_cmd_text(["samtools", "index", str(bam_p)])

        fs_text = run_cmd_text(["samtools", "flagstat", str(bam_p)])
        fs = parse_flagstat(fs_text)
        total = fs.get("total") or 0
        primary = fs.get("primary_mapped") or 0
        pct_primary = (primary / total * 100) if total > 0 else 0.0

        if show_flagstat:
            st.caption("`samtools flagstat` (raw)")
            st.code(fs_text)

        row = {
            "reads": reads_p.name,
            "ref": ref_p.name,
            "bam": str(bam_p),
            "pct_mapped_primary": round(pct_primary, 2),
        }

        # ---- bcftools pipeline ----
        if call_variants:
            vcf_gz, norm_vcf_gz, filt_vcf_gz, n_all, n_filt = bcftools_call_and_filter(
                ref_fa=ref_p,
                bam_p=bam_p,
                mapq_min=int(mapq_min),
                baseq_min=int(baseq_min),
                normalize=bool(normalize_vcf),
                ploidy_mode=ploidy_mode,
                qual_min=int(qual_min),
                require_pass=bool(require_pass),
            )

            st.caption("VCF outputs written to data/")
            st.code(f"{vcf_gz.name} (all) → {n_all} sites\n"
                    f"{(Path(norm_vcf_gz).name if norm_vcf_gz else '')}\n"
                    f"{filt_vcf_gz.name} (filtered) → {n_filt} HQ SNPs")

            row.update({
                "vcf_all": str(vcf_gz),
                "vcf_norm": str(norm_vcf_gz) if norm_vcf_gz else "",
                "vcf_filtered": str(filt_vcf_gz),
                "variants_all": n_all,
                "variants_filtered": n_filt,
            })

        results.append(row)

    # tables
    df = pd.DataFrame(results)
    st.subheader("Per-pair summary")
    st.dataframe(df, use_container_width=True)

    st.subheader("Primary-mapping % matrix (reads × ref)")
    mat = df.pivot(index="reads", columns="ref", values="pct_mapped_primary").fillna(0)
    st.dataframe(mat, use_container_width=True)

    st.success("Cross-alignments + bcftools variant calling complete.")
