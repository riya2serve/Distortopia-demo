import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import streamlit as st

# ---------- Config ----------
DEFAULT_INDIR = "alt_refs"
DEFAULT_OUTDIR = "alignments"

EXPECTED = {
    "ath_ref": "A_thaliana.fna",
    "ath_h1":  "A_thaliana_hap1.fa",
    "ath_h2":  "A_thaliana_hap2.fa",
    "aly_ref": "A_lyrata.fna",
    "aly_h1":  "A_lyrata_hap1.fa",
    "aly_h2":  "A_lyrata_hap2.fa",
}

# ---------- Helpers ----------
def check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def run(cmd, out_path=None):
    st.write("`$ " + " ".join(cmd) + "`")
    if out_path is None:
        subprocess.run(cmd, check=True)
    else:
        with open(out_path, "w") as fout:
            subprocess.run(cmd, check=True, stdout=fout)

def paf_summary(paf_path: str) -> dict:
    """
    PAF format mandatory cols:
      col10 = number of residue matches
      col11 = alignment block length
    We report totals and mean identity.
    """
    matches = 0
    alen = 0
    lines = 0
    with open(paf_path) as f:
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            try:
                m = int(parts[9])
                a = int(parts[10])
            except ValueError:
                continue
            matches += m
            alen += a
            lines += 1
    ident = (matches / alen) * 100 if alen > 0 else 0.0
    return {"records": lines, "aligned_bp": alen, "matches": matches, "identity_%": ident}

# ---------- UI ----------
st.title("Distortopia: Run Minimap2 Alignments")

st.subheader("Inputs")
indir = st.text_input("Directory with FASTAs", value=DEFAULT_INDIR)
outdir = st.text_input("Output directory", value=DEFAULT_OUTDIR)
threads = st.number_input("Threads", min_value=1, max_value=64, value=6)
preset = st.selectbox("Minimap2 preset", options=["asm5", "asm10", "asm20"], index=0)

st.subheader("File names (override if needed)")
colL, colR = st.columns(2)
with colL:
    ath_ref = st.text_input("A. thaliana reference", value=EXPECTED["ath_ref"])
    ath_h1  = st.text_input("A. thaliana hap1",     value=EXPECTED["ath_h1"])
    ath_h2  = st.text_input("A. thaliana hap2",     value=EXPECTED["ath_h2"])
with colR:
    aly_ref = st.text_input("A. lyrata reference",  value=EXPECTED["aly_ref"])
    aly_h1  = st.text_input("A. lyrata hap1",       value=EXPECTED["aly_h1"])
    aly_h2  = st.text_input("A. lyrata hap2",       value=EXPECTED["aly_h2"])

st.subheader("Which alignments?")
do_ref_h1 = st.checkbox("Reference vs hap1 (both species)", value=True)
do_ref_h2 = st.checkbox("Reference vs hap2 (both species)", value=True)
do_h2_h1  = st.checkbox("hap2 vs hap1 (both species)", value=True)
do_cross  = st.checkbox("Cross-species: A_thaliana_ref vs A_lyrata_ref (PAF)", value=False)

make_bam = st.checkbox("Also produce sorted+indexed BAM (requires samtools)", value=False)

if st.button("Run alignments", type="primary"):
    # Check tools
    if not check_tool("minimap2"):
        st.error("minimap2 not found on PATH. Install: `conda install -c bioconda minimap2`")
        st.stop()
    if make_bam and not check_tool("samtools"):
        st.error("samtools not found on PATH. Install: `conda install -c bioconda samtools`")
        st.stop()

    # Resolve paths
    indir_p = Path(indir)
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    files = {
        "ath_ref": str(indir_p / ath_ref),
        "ath_h1":  str(indir_p / ath_h1),
        "ath_h2":  str(indir_p / ath_h2),
        "aly_ref": str(indir_p / aly_ref),
        "aly_h1":  str(indir_p / aly_h1),
        "aly_h2":  str(indir_p / aly_h2),
    }
    missing = [k for k,v in files.items() if not Path(v).exists()]
    if missing:
        st.error(f"Missing files: {', '.join(missing)}")
        st.stop()

    # Build job list
    paf_jobs = []
    bam_jobs = []

    if do_h2_h1:
        paf_jobs += [
            (files["ath_h2"], files["ath_h1"], outdir_p / "ath_hap2_x_hap1.paf"),
            (files["aly_h2"], files["aly_h1"], outdir_p / "aly_hap2_x_hap1.paf"),
        ]
    if do_ref_h1:
        paf_jobs += [
            (files["ath_ref"], files["ath_h1"], outdir_p / "ath_ref_x_hap1.paf"),
            (files["aly_ref"], files["aly_h1"], outdir_p / "aly_ref_x_hap1.paf"),
        ]
        if make_bam:
            bam_jobs += [
                (files["ath_ref"], files["ath_h1"], outdir_p / "ath_ref_x_hap1.bam"),
                (files["aly_ref"], files["aly_h1"], outdir_p / "aly_ref_x_hap1.bam"),
            ]
    if do_ref_h2:
        paf_jobs += [
            (files["ath_ref"], files["ath_h2"], outdir_p / "ath_ref_x_hap2.paf"),
            (files["aly_ref"], files["aly_h2"], outdir_p / "aly_ref_x_hap2.paf"),
        ]
        if make_bam:
            bam_jobs += [
                (files["ath_ref"], files["ath_h2"], outdir_p / "ath_ref_x_hap2.bam"),
                (files["aly_ref"], files["aly_h2"], outdir_p / "aly_ref_x_hap2.bam"),
            ]

    # Run PAFs
    st.subheader("Running minimap2 (PAF)")
    results = []
    for ref, qry, paf_out in paf_jobs:
        st.write(f"**PAF:** `{Path(paf_out).name}`")
        with st.spinner(f"Aligning {Path(qry).name} → {Path(ref).name}"):
            run(["minimap2", "-x", preset, "-t", str(int(threads)), ref, qry], out_path=str(paf_out))
        summary = paf_summary(str(paf_out))
        summary.update({"paf": str(paf_out), "ref": Path(ref).name, "qry": Path(qry).name})
        results.append(summary)

    # Optional cross-species
    if do_cross:
        st.write("**PAF (cross-species):** `ath_ref_x_aly_ref.paf`")
        cross_preset = preset if preset != "asm5" else "asm10"
        with st.spinner("Aligning A_thaliana_ref ↔ A_lyrata_ref (PAF)"):
            outp = outdir_p / "ath_ref_x_aly_ref.paf"
            run(["minimap2", "-x", cross_preset, "-t", str(int(threads)), files["ath_ref"], files["aly_ref"]],
                out_path=str(outp))
        summary = paf_summary(str(outp))
        summary.update({"paf": str(outp), "ref": Path(files['ath_ref']).name, "qry": Path(files['aly_ref']).name})
        results.append(summary)

    # Run BAMs
    if make_bam and bam_jobs:
        st.subheader("Running minimap2+samtools (BAM)")
        for ref, qry, bam_out in bam_jobs:
            name = Path(bam_out).name
            with st.spinner(f"Aligning {Path(qry).name} → {Path(ref).name} -> {name}"):
                p1 = subprocess.Popen(["minimap2", "-ax", preset, "-t", str(int(threads)), ref, qry], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["samtools", "sort", "-o", str(bam_out)], stdin=p1.stdout)
                p1.stdout.close()
                rc = p2.wait()
                if rc != 0:
                    st.error(f"samtools sort failed for {name}")
                    continue
                run(["samtools", "index", str(bam_out)])

    # Show results table + downloads
    st.subheader("PAF summary")
    if results:
        df = pd.DataFrame(results)[["ref", "qry", "records", "aligned_bp", "matches", "identity_%", "paf"]]
        df["identity_%"] = df["identity_%"].map(lambda x: round(x, 4))
        st.dataframe(df, use_container_width=True)
        for _, row in df.iterrows():
            with open(row["paf"], "rb") as f:
                st.download_button(label=f"Download {Path(row['paf']).name}", data=f,
                                   file_name=Path(row["paf"]).name)

    st.success(f"Done. Outputs in: {outdir_p.resolve()}")

