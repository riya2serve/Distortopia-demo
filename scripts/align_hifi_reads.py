# scripts/align_hifi_reads.py
import shutil, subprocess
from pathlib import Path
from typing import Tuple, Optional, List, Dict
import pandas as pd
import streamlit as st

def check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def run_pipe_to_file(p1, p2, outfile):
    st.code("$ " + " ".join(p1) + " | " + " ".join(p2 + ["-o", outfile]))
    p_minimap = subprocess.Popen(p1, stdout=subprocess.PIPE)
    p_sort = subprocess.Popen(p2 + ["-o", outfile], stdin=p_minimap.stdout)
    p_minimap.stdout.close()
    rc = p_sort.wait()
    if rc != 0:
        raise RuntimeError("samtools sort failed")

def run_cmd_text(cmd) -> str:
    st.code("$ " + " ".join(cmd))
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out.stdout

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
    data_dir = Path("data")
    data_dir.mkdir(parents=True, exist_ok=True)
    base = bam_path.stem
    return (
        data_dir / f"{base}.vcf",
        data_dir / f"{base}.norm.vcf",
        data_dir / f"{base}.filtered.vcf",
    )

def ensure_fasta_index(ref_fa: Path):
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if not fai.exists():
        run_cmd_text(["samtools", "faidx", str(ref_fa)])

def ensure_bam_index(bam_p: Path):
    bai = Path(str(bam_p) + ".bai")
    if not bai.exists():
        run_cmd_text(["samtools", "index", str(bam_p)])
    run_cmd_text(["bash", "-lc", f"set -e; samtools quickcheck -v {bam_p} >/dev/null || exit 1"])

def vcf_count_records(vcf_path: Path) -> int:
    out = run_cmd_text(["bash", "-lc", f"grep -vc '^#' {vcf_path} || true"])
    try:
        return int(out.strip())
    except:
        return 0

def bcftools_call_and_filter(
    ref_fa: Path,
    bam_p: Path,
    mapq_min: int,
    baseq_min: int,
    normalize: bool,
    ploidy_n: int,            # numeric ploidy
    qual_min: int,
    require_pass: bool,
    keep_indels: bool = False # include INDELs in filtering step
) -> Tuple[Path, Optional[Path], Path, int, int]:
    if not check_tool("bcftools"):
        st.error("`bcftools` not found"); st.stop()
    ensure_fasta_index(ref_fa); ensure_bam_index(bam_p)
    vcf, norm_vcf, filt_vcf = auto_out_vcfs(bam_p)

    # mpileup → call → uncompressed VCF (-Ov), explicit ploidy and -m model
    with st.spinner(f"bcftools call → {vcf.name}"):
        anno = "FORMAT/AD,FORMAT/DP,FORMAT/SP"
        run_cmd_text([
            "bash","-lc",
            "set -o pipefail; "
            f"bcftools mpileup -f {ref_fa} -q {mapq_min} -Q {baseq_min} -Ou -a {anno} {bam_p} "
            f"| bcftools call -m --ploidy {ploidy_n} -v -Ov -o {vcf}"
        ])

    n_all = vcf_count_records(vcf)
    src_vcf = vcf

    if normalize:
        with st.spinner(f"bcftools norm → {norm_vcf.name}"):
            run_cmd_text(["bcftools", "norm", "-f", str(ref_fa), "-Ov", "-o", str(norm_vcf), str(vcf)])
        src_vcf = norm_vcf

    # filter
    with st.spinner(f"Filter high-quality variants → {filt_vcf.name}"):
        clauses = [f"QUAL>{qual_min}"]
        if require_pass:
            clauses.append('FILTER="PASS"')
        # keep diploid hets only (comment out if you don’t want this)
        if ploidy_n == 2:
            clauses.append('GT="het"')
        expr = " && ".join(clauses)

        variant_types = ["snps", "indels"] if keep_indels else ["snps"]

        run_cmd_text([
            "bcftools","view",
            "-v", ",".join(variant_types),
            "-m2","-M2",
            "-i", expr,
            "-Ov","-o", str(filt_vcf), str(src_vcf)
        ])

    n_filtered = vcf_count_records(filt_vcf)
    return vcf, (norm_vcf if normalize else None), filt_vcf, n_all, n_filtered

def run(
    ref1: Path, ref2: Path, reads1: Path, reads2: Path, *,
    threads: int, outdir: Path, index_bam: bool, call_variants: bool,
    ploidy_n: int, normalize_vcf: bool, require_pass: bool,
    mapq_min: int, baseq_min: int, qual_min: int, keep_indels: bool,
    show_flagstat: bool=False
) -> pd.DataFrame:
    for tool in ["minimap2","samtools"]:
        if not check_tool(tool):
            raise RuntimeError(f"`{tool}` not found")
    if call_variants and not check_tool("bcftools"):
        raise RuntimeError("`bcftools` not found")
    for p in [ref1,ref2,reads1,reads2]:
        if not Path(p).exists():
            raise FileNotFoundError(p)

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    ref_paths  = [Path(ref1), Path(ref2)]
    read_paths = [Path(reads1), Path(reads2)]
    jobs = [(r, ref, auto_out_bam(r, ref, outdir)) for r in read_paths for ref in ref_paths]
    rows: List[Dict] = []
    for reads_p, ref_p, bam_p in jobs:
        run_pipe_to_file(
            ["minimap2","-t",str(int(threads)),"-ax","map-hifi", str(ref_p), str(reads_p)],
            ["samtools","sort"], str(bam_p)
        )
        if index_bam:
            run_cmd_text(["samtools","index",str(bam_p)])
        fs_text = run_cmd_text(["samtools","flagstat",str(bam_p)])
        fs = parse_flagstat(fs_text)
        total = fs.get("total") or 0
        primary = fs.get("primary_mapped") or 0
        row = {
            "reads": reads_p.name,
            "ref": ref_p.name,
            "bam": str(bam_p),
            "pct_mapped_primary": round((primary/total*100) if total>0 else 0.0, 2)
        }
        if show_flagstat:
            row["flagstat"] = fs_text

        if call_variants:
            vcf, norm_vcf, filt_vcf, n_all, n_filt = bcftools_call_and_filter(
                ref_fa=ref_p,
                bam_p=bam_p,
                mapq_min=int(mapq_min),
                baseq_min=int(baseq_min),
                normalize=bool(normalize_vcf),
                ploidy_n=int(ploidy_n),
                qual_min=int(qual_min),
                require_pass=bool(require_pass),
                keep_indels=bool(keep_indels),
            )
            row.update({
                "vcf_all": str(vcf),
                "vcf_norm": str(norm_vcf) if norm_vcf else "",
                "vcf_filtered": str(filt_vcf),
                "variants_all": n_all,
                "variants_filtered": n_filt,
            })
        rows.append(row)
    return pd.DataFrame(rows)

def streamlit_panel(state):
    st.header("3) Cross-map long reads & call variants → data/ (uncompressed VCF)")
    col1,col2 = st.columns(2)
    with col1:
        ref1 = st.text_input("Reference FASTA #1", value=str(state.raw_data / "A_thaliana.fna"))
        reads1 = st.text_input("Reads FASTQ #1", value=str(state.long_reads / "thaliana_gametes_hifi.fq"))
    with col2:
        ref2 = st.text_input("Reference FASTA #2", value=str(state.raw_data / "A_lyrata.fna"))
        reads2 = st.text_input("Reads FASTQ #2", value=str(state.long_reads / "lyrata_gametes_hifi.fq"))

    threads = st.number_input("Threads", 1, 64, int(state.threads))
    outdir = Path(st.text_input("BAM output dir", value=str(state.alignments)))
    index_bam = st.checkbox("Index BAMs", True)
    show_flagstat = st.checkbox("Include flagstat text", False)

    st.subheader("bcftools settings")
    call_variants = st.checkbox("Call variants", True)
    c1,c2,c3 = st.columns(3)
    with c1:
        ploidy_n = st.selectbox("Ploidy", [1, 2], index=1)  # 2 = diploid
    with c2:
        normalize_vcf = st.checkbox("Normalize (bcftools norm)", True)
    with c3:
        require_pass = st.checkbox('Require FILTER="PASS"', True)

    c4,c5,c6 = st.columns(3)
    with c4:
        mapq_min = st.number_input("Read MAPQ ≥", 0, 60, 30)
    with c5:
        baseq_min = st.number_input("Base QUAL ≥", 0, 60, 30)
    with c6:
        qual_min = st.number_input("VCF QUAL >", 0, 500, 30)

    keep_indels = st.checkbox("Keep INDELs (not only SNPs)", value=False)

    if st.button("Run align+call", type="primary"):
        try:
            df = run(
                Path(ref1), Path(ref2), Path(reads1), Path(reads2),
                threads=int(threads), outdir=outdir, index_bam=bool(index_bam),
                call_variants=bool(call_variants), ploidy_n=int(ploidy_n),
                normalize_vcf=bool(normalize_vcf), require_pass=bool(require_pass),
                mapq_min=int(mapq_min), baseq_min=int(baseq_min), qual_min=int(qual_min),
                keep_indels=bool(keep_indels),
                show_flagstat=bool(show_flagstat),
            )
            st.subheader("Per-pair summary")
            st.dataframe(df, use_container_width=True)
            st.subheader("Primary-mapping % matrix")
            mat = df.pivot(index="reads", columns="ref", values="pct_mapped_primary").fillna(0)
            st.dataframe(mat, use_container_width=True)
            st.success("Done.")
        except Exception as e:
            st.error(str(e))
