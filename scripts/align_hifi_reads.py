# scripts/align_hifi_reads_cross_app.py
import shutil
import subprocess
from pathlib import Path
import pandas as pd
import streamlit as st

# ------------------ helpers ------------------
def check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def run_pipe_to_file(p1, p2, outfile):
    """minimap2 SAM -> samtools sort -> BAM"""
    st.code("$ " + " ".join(p1) + " | " + " ".join(p2 + ["-o", outfile]))
    p_minimap = subprocess.Popen(p1, stdout=subprocess.PIPE)
    p_sort = subprocess.Popen(p2 + ["-o", outfile], stdin=p_minimap.stdout)
    p_minimap.stdout.close()
    rc = p_sort.wait()
    if rc != 0:
        raise RuntimeError("samtools sort failed")

def run_cmd_text(cmd):
    st.code("$ " + " ".join(cmd))
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out.stdout

def parse_flagstat(txt: str) -> dict:
    """
    Robust-ish parser across samtools versions.
    We want:
      - total reads
      - mapped (includes secondary/supplementary)
      - primary mapped (unique primary alignments)
    """
    d = {"total": None, "mapped": None, "primary_mapped": None}
    for line in txt.splitlines():
        ls = line.strip()
        if " in total" in ls and d["total"] is None:
            # e.g. "1016 + 0 in total ..."
            d["total"] = int(ls.split()[0])
        elif " mapped (" in ls and "primary" not in ls and d["mapped"] is None:
            # e.g. "1016 + 0 mapped (100.00% : N/A)"
            d["mapped"] = int(ls.split()[0])
        elif "primary mapped" in ls and d["primary_mapped"] is None:
            # e.g. "1000 + 0 primary mapped (98.43% : N/A)"
            d["primary_mapped"] = int(ls.split()[0])
    return d

def auto_out_bam(reads_path: Path, ref_path: Path, outdir: Path) -> Path:
    r = reads_path.stem  # .fq or .fq.gz → "stem" will give "fq" if .gz; handle nicer:
    if reads_path.name.endswith(".fq.gz") or reads_path.name.endswith(".fastq.gz"):
        r = reads_path.name.rsplit(".fq.gz", 1)[0].rsplit(".fastq.gz", 1)[0]
    elif reads_path.suffix in (".fq", ".fastq"):
        r = reads_path.stem
    ref = ref_path.stem
    return outdir / f"{r}__to__{ref}.bam"

# ------------------ UI ------------------
st.title("Distortopia: Cross-map HiFi reads to references (2×2 matrix)")

with st.expander("Requirements", expanded=False):
    st.markdown("""
- **minimap2** (preset `map-hifi`)
- **samtools** (for `sort`, `index`, `flagstat`)
- Inputs may be plain `.fq/.fastq` or gzipped `.fq.gz/.fastq.gz`.
""")

# Inputs
col1, col2 = st.columns(2)
with col1:
    ref1 = st.text_input("Reference FASTA #1", value="raw_data/A_thaliana.fna")
    reads1 = st.text_input("Reads FASTQ #1", value="reads/Athal_gametes_hifi.fq")
with col2:
    ref2 = st.text_input("Reference FASTA #2", value="raw_data/A_lyrata.fna")
    reads2 = st.text_input("Reads FASTQ #2", value="reads/Aly_gametes_hifi.fq")

threads = st.number_input("Threads", min_value=1, max_value=64, value=8)
outdir  = Path(st.text_input("Output directory for BAMs", value="alignments"))
index_bam = st.checkbox("Index BAMs (samtools index)", value=True)
show_flagstat = st.checkbox("Show per-pair flagstat", value=False)

if st.button("Run cross-alignments", type="primary"):
    # Tool checks
    if not check_tool("minimap2"):
        st.error("`minimap2` not found on PATH. Try: `conda install -c bioconda minimap2`")
        st.stop()
    if not check_tool("samtools"):
        st.error("`samtools` not found on PATH. Try: `conda install -c bioconda samtools`")
        st.stop()

    # Resolve paths and existence
    ref_paths  = [Path(ref1), Path(ref2)]
    read_paths = [Path(reads1), Path(reads2)]
    for p in ref_paths + read_paths:
        if not p.exists():
            st.error(f"Missing file: {p}")
            st.stop()

    outdir.mkdir(parents=True, exist_ok=True)

    # Build the 2×2 jobs: each reads vs each ref
    jobs = []
    for r in read_paths:
        for ref in ref_paths:
            bam = auto_out_bam(r, ref, outdir)
            jobs.append((r, ref, bam))

    results = []  # one row per pair

    for reads_p, ref_p, bam_p in jobs:
        st.subheader(f"Map {reads_p.name} → {ref_p.name}")
        with st.spinner(f"minimap2 map-hifi → {bam_p.name}"):
            p1 = ["minimap2", "-t", str(int(threads)), "-ax", "map-hifi", str(ref_p), str(reads_p)]
            p2 = ["samtools", "sort"]
            run_pipe_to_file(p1, p2, str(bam_p))

        if index_bam:
            with st.spinner(f"Index {bam_p.name}"):
                run_cmd_text(["samtools", "index", str(bam_p)])

        # flagstat
        fs = parse_flagstat(run_cmd_text(["samtools", "flagstat", str(bam_p) ]))
        total = fs.get("total") or 0
        mapped = fs.get("mapped") or 0
        primary = fs.get("primary_mapped") or 0
        pct_primary = (primary / total * 100) if total > 0 else 0.0

        if show_flagstat:
            st.caption("`samtools flagstat` (raw)")
            st.code(run_cmd_text(["samtools", "flagstat", str(bam_p)]))

        results.append({
            "reads": reads_p.name,
            "ref": ref_p.name,
            "bam": str(bam_p),
            "total_reads": total,
            "mapped_reads": mapped,
            "mapped_primary": primary,
            "pct_mapped_primary": round(pct_primary, 2),
        })

    # Detailed per-pair table
    st.subheader("Per-pair summary")
    df = pd.DataFrame(results)
    st.dataframe(df, use_container_width=True)

    # 2×2 matrix (rows: reads; cols: ref; values: % primary mapped)
    st.subheader("Primary-mapping % matrix (reads × reference)")
    mat = df.pivot(index="reads", columns="ref", values="pct_mapped_primary").fillna(0)
    st.dataframe(mat, use_container_width=True)

    # Download buttons
    st.subheader("Download BAMs")
    for row in results:
        bam = Path(row["bam"])
        if bam.exists():
            with open(bam, "rb") as fh:
                st.download_button(f"Download {bam.name}", fh, file_name=bam.name)
        bai = bam.with_suffix(bam.suffix + ".bai")
        if bai.exists():
            with open(bai, "rb") as fh:
                st.download_button(f"Download {bai.name}", fh, file_name=bai.name)

    st.success("Cross-alignments complete.")

