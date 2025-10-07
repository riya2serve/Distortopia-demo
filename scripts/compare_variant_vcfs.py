# scripts/compare_variant_vcfs.py
from __future__ import annotations
import shutil
import subprocess
from pathlib import Path
from typing import Tuple, List

import pandas as pd
import streamlit as st

# ---------- helpers ----------
def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    st.code("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=check, capture_output=True, text=True)

def ensure_bgzip_tbi(vcf_path: str) -> str:
    """
    Accept .vcf/.vcf.gz/.bcf and return a bgzipped (and indexed) path usable by bcftools/tabix.
    If input is BCF, ensure .csi exists and return the .bcf path.
    """
    p = Path(vcf_path)
    if not p.exists():
        raise FileNotFoundError(vcf_path)

    # If BCF, ensure CSI and return
    if p.suffix == ".bcf":
        csi = p.with_suffix(p.suffix + ".csi")
        if not csi.exists():
            run(["bcftools", "index", "-f", str(p)])
        return str(p)

    # If .vcf.gz already, just ensure .tbi
    if p.suffixes[-2:] == [".vcf", ".gz"]:
        gz = p
    elif p.suffix == ".vcf":
        # bgzip the plain VCF
        gz = p.with_suffix(".vcf.gz")
        with open(gz, "wb") as out:
            subprocess.run(["bgzip", "-c", str(p)], check=True, stdout=out)
    else:
        # Unexpected suffix; try to pass through
        gz = p

    # tabix index if missing
    tbi = Path(str(gz) + ".tbi")
    if not tbi.exists():
        run(["tabix", "-p", "vcf", str(gz)])
    return str(gz)

def count_records(vcf_like: str) -> int:
    # Works for .vcf.gz and .bcf
    cp = subprocess.run(["bcftools", "view", "-H", vcf_like],
                        check=True, capture_output=True, text=True)
    return 0 if not cp.stdout else len(cp.stdout.splitlines())

# ---------- unified-app panel ----------
def streamlit_panel(state) -> None:
    """
    Render the 'Compare VCFs' tab inside the unified Streamlit app.
    `state` is the Namespace-like object you pass around with common paths.
    """
    st.header("5) Compare two VCFs (shared vs unique)")

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
        vcf1 = st.text_input(
            "VCF #1",
            value=str(state.data / "Athaliana_gametes_hifi__to__A_thaliana.filtered.vcf"),
            key="isec_vcf1",
        )
        label1 = st.text_input("Label for VCF #1", value="A_thaliana", key="isec_label1")
    with col2:
        vcf2 = st.text_input(
            "VCF #2",
            value=str(state.data / "Alyrata_gametes_hifi__to__A_lyrata.filtered.vcf"),
            key="isec_vcf2",
        )
        label2 = st.text_input("Label for VCF #2", value="A_lyrata", key="isec_label2")

    out_dir = Path(st.text_input("Output directory", value="compare_vcfs", key="isec_outdir"))
    go = st.button("Run comparison", type="primary", key="isec_go")

    if not go:
        return

    # tool checks
    for tool in ["bcftools", "bgzip", "tabix"]:
        if not have(tool):
            st.error(f"`{tool}` not found on PATH. Install via conda (bioconda).")
            return

    out_dir.mkdir(parents=True, exist_ok=True)

    # Prepare inputs
    try:
        v1 = ensure_bgzip_tbi(vcf1)
        v2 = ensure_bgzip_tbi(vcf2)
    except Exception as e:
        st.error(str(e))
        return

    # bcftools isec
    work = out_dir / f"{Path(v1).stem}__vs__{Path(v2).stem}"
    work.mkdir(parents=True, exist_ok=True)
    st.subheader("Running bcftools isec")
    rc = run(["bcftools", "isec", "-p", str(work), "-Oz", v1, v2])

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
            out_vcf = work / f"{nice}.vcf.gz"
            shutil.copyfile(path, out_vcf)
            run(["tabix", "-p", "vcf", str(out_vcf)])
            n = count_records(str(out_vcf))
            rows.append({"subset": nice, "records": n, "vcf": str(out_vcf)})
            downloads.append(out_vcf)

    if not rows:
        st.warning("No records produced. Are the two VCFs in the same contig/coordinate space?")
        return

    st.subheader("Summary")
    df = pd.DataFrame(rows)[["subset", "records", "vcf"]]
    st.dataframe(df, use_container_width=True)

    st.subheader("Download VCFs")
    for v in downloads:
        with open(v, "rb") as fh:
            st.download_button(f"Download {v.name}", fh, file_name=v.name, key=f"dl_{v.name}")

    st.success(f"Done. Results in: {work.resolve()}")

# ---------- keep standalone mode working ----------
def _standalone():
    class _State:
        # sensible defaults if run directly
        data = Path("data")
    streamlit_panel(_State())

if __name__ == "__main__":
    _standalone()

