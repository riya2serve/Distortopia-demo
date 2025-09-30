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

# --- NEW: parse PAF optional tags like dv:f:0.001, tp:A:P, cm:i:..., etc.
def parse_tags(fields):
    tags = {}
    for f in fields:
        parts = f.split(":", 2)
        if len(parts) < 3:
            continue
        tag, typ, val = parts
        if typ == "i":
            try: val = int(val)
            except: pass
        elif typ == "f":
            try: val = float(val)
            except: pass
        # A/Z/H stay as strings
        tags[tag] = val
    return tags

# --- UPDATED: prefer dv (divergence) when available for identity
def paf_summary(paf_path: str) -> dict:
    """
    PAF format mandatory cols:
      col10 = number of residue matches
      col11 = alignment block length
    We report totals and mean identity (prefer dv:f: tag when available).
    """
    matches = 0
    alen = 0
    lines = 0
    pid_sum = 0.0
    pid_n   = 0

    with open(paf_path) as f:
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            try:
                m = int(parts[9])     # nmatch
                a = int(parts[10])    # aln_len
            except ValueError:
                continue

            matches += m
            alen    += a
            lines   += 1

            # optional tags for identity
            tags = parse_tags(parts[12:])
            dv = tags.get("dv", None)
            if isinstance(dv, float):
                pid_sum += max(0.0, min(1.0, 1.0 - dv))
                pid_n   += 1

    # mean identity: prefer dv-derived; else nmatch/aln_len
    if pid_n > 0:
        ident = (pid_sum / pid_n) * 100.0
    else:
        ident = (matches / alen) * 100.0 if alen > 0 else 0.0

    return {"records": lines, "aligned_bp": alen, "matches": matches, "identity_%": ident}

# --- NEW: write a labeled TSV for every PAF row (human-friendly)
def paf_to_tsv(paf_path: str, out_tsv: str):
    headers = [
        "qname","qlen","qstart","qend","strand",
        "tname","tlen","tstart","tend",
        "nmatch","aln_len","mapq",
        "tp","cm","s1","s2","dv","rl",
        "pct_identity_est"
    ]
    n = 0
    with open(paf_path) as fin, open(out_tsv, "w") as fout:
        fout.write("\t".join(headers) + "\n")
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue

            qname  = f[0]
            qlen   = int(f[1])
            qstart = int(f[2])
            qend   = int(f[3])
            strand = f[4]
            tname  = f[5]
            tlen   = int(f[6])
            tstart = int(f[7])
            tend   = int(f[8])
            nmatch = int(f[9])
            alen   = int(f[10])
            mapq   = int(f[11])

            tags = parse_tags(f[12:])
            tp = tags.get("tp","")
            cm = tags.get("cm","")
            s1 = tags.get("s1","")
            s2 = tags.get("s2","")
            dv = tags.get("dv","")
            rl = tags.get("rl","")

            if isinstance(dv, float):
                pid = max(0.0, min(100.0, (1.0 - dv) * 100.0))
            else:
                pid = (nmatch / alen * 100.0) if alen > 0 else ""

            row = [
                qname, qlen, qstart, qend, strand,
                tname, tlen, tstart, tend,
                nmatch, alen, mapq,
                tp, cm, s1, s2, dv, rl,
                f"{pid:.4f}" if pid != "" else ""
            ]
            fout.write("\t".join(map(str,row)) + "\n")
            n += 1
    return n

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

# --- NEW toggles for tables
emit_row_tsv   = st.checkbox("Write row-level TSV for each PAF (…paf.tsv)", value=True)
emit_rollup_tsv= st.checkbox("Also write a roll-up TSV (paf_summary.tsv)", value=True)

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
    rollup_path = outdir_p / "paf_summary.tsv"
    if emit_rollup_tsv and rollup_path.exists():
        # start a fresh roll-up
        rollup_path.unlink()

    for ref, qry, paf_out in paf_jobs:
        st.write(f"**PAF:** `{Path(paf_out).name}`")
        with st.spinner(f"Aligning {Path(qry).name} → {Path(ref).name}"):
            run(["minimap2", "-x", preset, "-t", str(int(threads)), ref, qry], out_path=str(paf_out))

        # row-level TSV
        if emit_row_tsv:
            tsv_path = str(paf_out) + ".tsv"
            nrows = paf_to_tsv(str(paf_out), tsv_path)
            st.caption(f"→ wrote {nrows} rows to `{Path(tsv_path).name}`")

        # per-file summary (in-memory)
        summary = paf_summary(str(paf_out))
        summary.update({"paf": str(paf_out), "ref": Path(ref).name, "qry": Path(qry).name})
        results.append(summary)

        # append to roll-up TSV on disk
        if emit_rollup_tsv:
            hdr = "paf\tref\tqry\trecords\taligned_bp\tmatches\tidentity_%\n"
            line = "\t".join([
                str(paf_out), Path(ref).name, Path(qry).name,
                str(summary["records"]), str(summary["aligned_bp"]),
                str(summary["matches"]), f"{summary['identity_%']:.6f}"
            ]) + "\n"
            write_header = not rollup_path.exists()
            with open(rollup_path, "a") as fo:
                if write_header:
                    fo.write(hdr)
                fo.write(line)

    # Optional cross-species
    if do_cross:
        st.write("**PAF (cross-species):** `ath_ref_x_aly_ref.paf`")
        cross_preset = preset if preset != "asm5" else "asm10"
        with st.spinner("Aligning A_thaliana_ref ↔ A_lyrata_ref (PAF)"):
            outp = outdir_p / "ath_ref_x_aly_ref.paf"
            run(["minimap2", "-x", cross_preset, "-t", str(int(threads)), files["ath_ref"], files["aly_ref"]],
                out_path=str(outp))

        if emit_row_tsv:
            tsv_path = str(outp) + ".tsv"
            nrows = paf_to_tsv(str(outp), tsv_path)
            st.caption(f"→ wrote {nrows} rows to `{Path(tsv_path).name}`")

        summary = paf_summary(str(outp))
        summary.update({"paf": str(outp), "ref": Path(files['ath_ref']).name, "qry": Path(files['aly_ref']).name})
        results.append(summary)

        if emit_rollup_tsv:
            hdr = "paf\tref\tqry\trecords\taligned_bp\tmatches\tidentity_%\n"
            line = "\t".join([
                str(outp), Path(files['ath_ref']).name, Path(files['aly_ref']).name,
                str(summary["records"]), str(summary["aligned_bp"]),
                str(summary["matches"]), f"{summary['identity_%']:.6f}"
            ]) + "\n"
            write_header = not rollup_path.exists()
            with open(rollup_path, "a") as fo:
                if write_header:
                    fo.write(hdr)
                fo.write(line)

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

        # download buttons for PAFs and any TSVs we wrote
        for _, row in df.iterrows():
            paf = row["paf"]
            with open(paf, "rb") as f:
                st.download_button(label=f"Download {Path(paf).name}", data=f, file_name=Path(paf).name)
            tsv = paf + ".tsv"
            if emit_row_tsv and Path(tsv).exists():
                with open(tsv, "rb") as tf:
                    st.download_button(label=f"Download {Path(tsv).name}", data=tf, file_name=Path(tsv).name)

        # roll-up download
        if emit_rollup_tsv and rollup_path.exists():
            with open(rollup_path, "rb") as rf:
                st.download_button(label=f"Download {rollup_path.name}", data=rf, file_name=rollup_path.name)

    st.success(f"Done. Outputs in: {outdir_p.resolve()}")

