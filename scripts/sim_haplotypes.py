import streamlit as st
import subprocess
import os
import pathlib

def run_haptools(input_fasta, output_prefix, snp_rate=0.001, indel_rate=0.0, seeds=(42, 84)):
    """Generate two haploid genomes from a reference FASTA using haptools mutate."""
    out1 = f"{output_prefix}_hap1.fa"
    out2 = f"{output_prefix}_hap2.fa"

    # Hap1
    subprocess.run([
        "haptools", "mutate",
        "--input", input_fasta,
        "--output", out1,
        "--snp-rate", str(snp_rate),
        "--indel-rate", str(indel_rate),
        "--seed", str(seeds[0])
    ], check=True)

    # Hap2
    subprocess.run([
        "haptools", "mutate",
        "--input", input_fasta,
        "--output", out2,
        "--snp-rate", str(snp_rate),
        "--indel-rate", str(indel_rate),
        "--seed", str(seeds[1])
    ], check=True)

    return out1, out2


# ---------------------- Streamlit UI ---------------------- #
st.title("Distortopia: Simulate Haplotypes with Haptools")

# Inputs
ref1 = st.file_uploader("Upload Reference Genome 1 (FASTA)", type=["fa", "fasta", "fna"])
ref2 = st.file_uploader("Upload Reference Genome 2 (FASTA)", type=["fa", "fasta", "fna"])

snp_rate = st.number_input("SNP mutation rate", min_value=0.0, max_value=0.01,
                           value=0.001, step=0.0001)
indel_rate = st.number_input("Indel mutation rate", min_value=0.0, max_value=0.01,
                             value=0.0, step=0.0001)
outdir = st.text_input("Output directory", value="alt_refs")

if st.button("Generate Haplotypes"):
    os.makedirs(outdir, exist_ok=True)

    for ref in [ref1, ref2]:
        if ref is not None:
            # Base output name from uploaded file (strip extension)
            prefix = pathlib.Path(ref.name).stem
            ref_path = os.path.join(outdir, ref.name)
            with open(ref_path, "wb") as f:
                f.write(ref.getbuffer())

            h1, h2 = run_haptools(ref_path, os.path.join(outdir, prefix),
                                  snp_rate, indel_rate)

            st.success(f"Generated haplotypes for {ref.name}: {h1}, {h2}")

            # Add download buttons
            with open(h1, "rb") as f1, open(h2, "rb") as f2:
                st.download_button("Download hap1", f1, file_name=os.path.basename(h1))
                st.download_button("Download hap2", f2, file_name=os.path.basename(h2))

