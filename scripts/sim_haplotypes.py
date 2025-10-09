# scripts/sim_haplotypes.py
from pathlib import Path
import streamlit as st
from scripts.mutate import mutate_fasta

def run(input_fasta: str, output_prefix: str, snp_rate: float = 0.001,
        indel_rate: float = 0.0, seeds=(42, 84)):
    out1 = f"{output_prefix}_hap1.fa"
    out2 = f"{output_prefix}_hap2.fa"
    mutate_fasta(input_fasta, out1, snp_rate=float(snp_rate), indel_rate=float(indel_rate), seed=int(seeds[0]))
    mutate_fasta(input_fasta, out2, snp_rate=float(snp_rate), indel_rate=float(indel_rate), seed=int(seeds[1]))
    return out1, out2

def streamlit_panel(state):
    st.header("1) Simulate haplotypes with Haptools")
    outdir = Path(st.text_input("Output directory", value=str(state.alt_refs)))
    snp_rate = st.number_input("SNP mutation rate", 0.0, 0.01, 0.001, 0.0001, format="%.4f")
    indel_rate = st.number_input("Indel mutation rate", 0.0, 0.01, 0.0, 0.0001, format="%.4f")
    seeds_str = st.text_input("Seeds (hap1,hap2)", value="42,84")
    seeds = tuple(int(x.strip()) for x in seeds_str.split(","))

    colA, colB = st.columns(2)
    with colA:
        ref1_path = st.text_input("Reference 1 FASTA", value=str(state.raw_data / "A_thaliana.fna"))
    with colB:
        ref2_path = st.text_input("Reference 2 FASTA", value=str(state.raw_data / "A_lyrata.fna"))

    if st.button("Generate haplotypes"):
        outdir.mkdir(parents=True, exist_ok=True)
        for label, ref in [("Ref1", ref1_path), ("Ref2", ref2_path)]:
            prefix = outdir / Path(ref).stem
            try:
                h1, h2 = run(ref, str(prefix), snp_rate, indel_rate, seeds)
                st.success(f"{label}: {Path(h1).name}, {Path(h2).name}")
            except Exception as e:
                st.error(f"{label} failed — {e}")

