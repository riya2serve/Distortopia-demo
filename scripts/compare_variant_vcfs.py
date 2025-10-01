import shutil
import subprocess
from pathlib import Path
import tempfile

import pandas as pd
import streamlit as st

# ---------- helpers ----------
def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run(cmd, check=True):
    st.code("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=check)

def ensure_bgzip_tbi(vcf_path: str) -> str:
    """
    Ensure VCF is bgzipped and tabix-indexed. Return .vcf.gz path.
    Accepts .vcf/.vcf.gz/.bcf; if BCF, leaves as-is but ensures .csi exists.
    """
    p = Path(vcf_path)
    if not p.exists():
        raise FileNotFoundError(vcf_path)

    if p.suffix == ".bcf":
        # Ensure BCF index
        csi = p.with_suffix(p.suffix + ".csi")
        if not csi.exists():
            run(["bcftools", "index", "-f", str(p)])
        return str(p)

    if p.suffixes[-2:] == [".vcf", ".gz"]:
        gz = p
    elif p.suffix == ".vcf":
        gz = p.with_suffix(".vcf.gz")
        run(["bgzip", "-c", str(p)], check=True).stdout
        # Above returns None; do with shell redir alternative:
        with open(gz, "wb") as out:
            subprocess.run(["bgzip", "-c", str(p)], check=True, stdout=out)
    else:
        # something else; try to pass through
        gz = p

    # tabix if missing
    tbi = Path(str(gz) + ".tbi")
    if not tbi.exists():
        run(["tabix", "-p", "vcf", str(gz)])
    return str(gz)

def count_records(vcf_gz: str) -> int:
    p = subprocess.run(["bcftools", "view", "-H", vcf_gz],
                       check=True, capture_output=True, text=True)
    return 0 if not p.stdout else len(p.stdout.splitlines())

# ---------- UI ----------
st.title("Distortopia: Compare two VCFs (shared vs unique)")

with st.expander("What this does", expanded=False):
    st.markdown(
        "- Ensures both VCFs are **bgzipped + indexed**\n"
        "- Runs **`bcftools isec`** to create:\n"
        "  - unique to File1\n"
        "  - unique to File2\n"
        "  - shared (present in both)\n"
        "- Shows counts and lets you **download the result VCFs**"
    )

col1, col2 = st.columns(2)
with col1:
    vcf1 = st.text_input("VCF #1 (Athal)", value="variants/Athal_bcftools.vcf")
    label1 = st.text_input("Label for VCF #1", value="A_thaliana")
with col2:
    vcf2 = st.text_input("VCF #2 (Aly)", value="variants/Aly_bcftools.vcf")
    label2 = st.text_input("Label for VCF #2", value="A_lyrata")

out_dir = st.text_input("Output directory", value="compare_vcfs")
go = st.button("Run comparison", type="primary")

if go:
    # tool checks
    for tool in ["bcftools", "bgzip", "tabix"]:
        if not have(tool):
            st.error(f"`{tool}` not found on PATH. Install via conda (bioconda).")
            st.stop()

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Prepare inputs
    try:
        v1 = ensure_bgzip_tbi(vcf1)
        v2 = ensure_bgzip_tbi(vcf2)
    except Exception as e:
        st.error(str(e))
        st.stop()

    # bcftools isec
    work = Path(out_dir) / f"{Path(v1).stem}__vs__{Path(v2).stem}"
    work.mkdir(parents=True, exist_ok=True)
    st.subheader("Running bcftools isec")
    run(["bcftools", "isec", "-p", str(work), "-Oz", v1, v2])

    # Map masks → labels
    mask_to_label = {
        "0010.vcf.gz": f"unique_to_{label1}",
        "0001.vcf.gz": f"unique_to_{label2}",
        "0011.vcf.gz": "shared",
    }

    rows = []
    downloads = []
    for fn, nice in mask_to_label.items():
        path = work / fn
        if path.exists():
            # rename a copy with nicer name
            out_vcf = work / f"{nice}.vcf.gz"
            shutil.copyfile(path, out_vcf)
            # index to be nice
            run(["tabix", "-p", "vcf", str(out_vcf)])
            n = count_records(str(out_vcf))
            rows.append({"subset": nice, "records": n, "vcf": str(out_vcf)})
            downloads.append(out_vcf)

    if not rows:
        st.warning("No records produced. Are the two VCFs in the same coordinate/contig space?")
        st.stop()

    st.subheader("Summary")
    df = pd.DataFrame(rows)[["subset", "records", "vcf"]]
    st.dataframe(df, use_container_width=True)

    # Downloads
    st.subheader("Download VCFs")
    for v in downloads:
        with open(v, "rb") as fh:
            st.download_button(f"Download {v.name}", fh, file_name=v.name)
    st.success(f"Done. Results in: {work.resolve()}")

