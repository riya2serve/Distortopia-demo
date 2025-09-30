# scripts/align_gamete_reads_app.py
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import streamlit as st

# ------------------ helpers ------------------
def check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def run_and_capture(cmd):
    """Run a command and capture stdout as text (raises on non-zero)."""
    st.code("$ " + " ".join(cmd))
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out.stdout

def run_pipe_to_file(p1, p2, outfile):
    """Pipe minimap2 SAM to 'samtools sort -o outfile'."""
    st.code("$ " + " ".join(p1) + " | " + " ".join(p2 + ["-o", outfile]))
    p_minimap = subprocess.Popen(p1, stdout=subprocess.PIPE)
    p_sort = subprocess.Popen(p2 + ["-o", outfile], stdin=p_minimap.stdout)
    p_minimap.stdout.close()
    rc = p_sort.wait()
    if rc != 0:
        raise RuntimeError("samtools sort failed")

def flagstat(bam_path: str) -> dict:
    txt = run_and_capture(["samtools", "flagstat", bam_path])
    # Return whole text plus a few extracted numbers
    lines = txt.strip().splitlines()
    out = {"flagstat_raw": txt}
    try:
        # first line example: "10000 + 0 in total (QC-passed reads + QC-failed reads)"
        parts0 = lines[0].split()
        out["total_reads"] = int(parts0[0])
        # second line example: "9999 + 0 mapped (99.99% : N/A)"
        out["mapped_reads"] = int(lines[4].split()[0]) if len(lines) > 4 else None
    except Exception:
        pass
    return out

# ------------------ UI ------------------
st.title("Distortopia: Map Gamete HiFi Reads to Reference")

with st.expander("Requirements", expanded=False):
    st.markdown("- **minimap2** (use preset `map-hifi`)\n- **samtools**\n- Reads are FASTQ (plain or .gz)")

col1, col2 = st.columns(2)
with col1:
    ref1 = st.text_input("Reference FASTA #1", value="raw_data/A_thaliana.fna")
    reads1 = st.text_input("Reads FASTQ #1", value="reads/Athal_gametes_hifi.fq")
    out1 = st.text_input("Output BAM #1", value="alignments/Athal_gametes.bam")
with col2:
    ref2 = st.text_input("Reference FASTA #2 (optional)", value="raw_data/A_lyrata.fna")
    reads2 = st.text_input("Reads FASTQ #2 (optional)", value="reads/Aly_gametes_hifi.fq")
    out2 = st.text_input("Output BAM #2 (optional)", value="alignments/Aly_gametes.bam")

threads = st.number_input("Threads", min_value=1, max_value=64, value=8)
make_index = st.checkbox("Index BAM (samtools index)", value=True)
show_flagstat = st.checkbox("Show samtools flagstat summary", value=True)

if st.button("Run alignments", type="primary"):
    # tool checks
    if not check_tool("minimap2"):
        st.error("`minimap2` not found on PATH. Try: `conda install -c bioconda minimap2`")
        st.stop()
    if not check_tool("samtools"):
        st.error("`samtools` not found on PATH. Try: `conda install -c bioconda samtools`")
        st.stop()

    jobs = []
    if ref1 and reads1 and out1:
        jobs.append((ref1, reads1, out1))
    if ref2 and reads2 and out2:
        jobs.append((ref2, reads2, out2))

    if not jobs:
        st.warning("No jobs configured.")
        st.stop()

    # ensure output dirs
    for _, _, bam_out in jobs:
        Path(bam_out).parent.mkdir(parents=True, exist_ok=True)

    summaries = []
    for ref, reads, bam_out in jobs:
        ref_p, reads_p, bam_p = Path(ref), Path(reads), Path(bam_out)
        if not ref_p.exists():
            st.error(f"Missing reference: {ref}")
            st.stop()
        if not reads_p.exists():
            st.error(f"Missing reads: {reads}")
            st.stop()

        name = bam_p.name
        st.subheader(f"Mapping: {reads_p.name} → {ref_p.name}")
        with st.spinner(f"Running minimap2 (map-hifi) → {name}"):
            # minimap2 to samtools sort
            p1 = ["minimap2", "-t", str(int(threads)), "-ax", "map-hifi", str(ref_p), str(reads_p)]
            p2 = ["samtools", "sort"]
            run_pipe_to_file(p1, p2, str(bam_p))

        if make_index:
            with st.spinner(f"Indexing {name}"):
                st.code("$ samtools index " + str(bam_p))
                subprocess.run(["samtools", "index", str(bam_p)], check=True)

        summary = {"bam": str(bam_p), "ref": ref_p.name, "reads": reads_p.name}
        if show_flagstat:
            with st.spinner(f"samtools flagstat {name}"):
                fs = flagstat(str(bam_p))
            summary.update({
                "total_reads": fs.get("total_reads"),
                "mapped_reads": fs.get("mapped_reads"),
            })
            st.text_area(f"flagstat: {name}", fs["flagstat_raw"], height=180)
        summaries.append(summary)

    # show small table and download buttons
    df = pd.DataFrame(summaries)
    st.subheader("Summary")
    st.dataframe(df, use_container_width=True)

    for s in summaries:
        with open(s["bam"], "rb") as fh:
            st.download_button(f"Download {Path(s['bam']).name}", fh, file_name=Path(s["bam"]).name)
        bai = s["bam"] + ".bai"
        if Path(bai).exists():
            with open(bai, "rb") as fh:
                st.download_button(f"Download {Path(bai).name}", fh, file_name=Path(bai).name)

    st.success("All alignments completed.")

