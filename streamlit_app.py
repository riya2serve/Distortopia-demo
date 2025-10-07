# streamlit_app.py
import streamlit as st
from pathlib import Path

import scripts.sim_haplotypes as sim_haplotypes
import scripts.sim_gamete_long_reads as sim_reads
import scripts.align_hifi_reads as align_hifi_reads
import scripts.haplotag_helper as haplotag_helper
import scripts.find_crossovers as find_crossovers
import scripts.visualize_crossovers as visualize_crossovers
try:
    import scripts.compare_variant_vcfs as compare_variant_vcfs
except Exception:
    compare_variant_vcfs = None

st.set_page_config(page_title="Distortopia pipeline", layout="wide")
st.title("Distortopia: one-window pipeline")

state = st.session_state
if "root" not in state:
    state.root       = Path(".").resolve()
    state.raw_data   = state.root / "raw_data"
    state.alt_refs   = state.root / "alt_refs"
    state.long_reads = state.root / "long_reads"
    state.alignments = state.root / "alignments"
    state.data       = state.root / "data"
    state.plots      = state.root / "plots"
    for p in (state.raw_data, state.alt_refs, state.long_reads, state.alignments, state.data, state.plots):
        p.mkdir(parents=True, exist_ok=True)

with st.sidebar:
    st.header("Global settings")
    state.threads     = st.number_input("Threads", 1, 64, value=8)
    state.read_len    = st.number_input("Read length (bp)", 1_000, 2_000_000, value=100_000, step=1_000)
    state.target_cov  = st.number_input("Target coverage (×)", 1.0, 200.0, value=50.0, step=1.0)
    state.mean_xovers = st.number_input("Mean crossovers per read", 0.0, 5.0, value=0.1, step=0.1)
    state.min_snps    = st.number_input("Min phased SNPs/read", 1, 100, value=5, step=1)

tabs = st.tabs([
    "1) Haplotypes",
    "2) Simulate reads",
    "3) Align & call",
    "4) Haplotag",
    "5) Crossovers",
    "6) Visualize",
    "7) (optional) Compare VCFs"
])

with tabs[0]: sim_haplotypes.streamlit_panel(state)
with tabs[1]: sim_reads.streamlit_panel(state)
with tabs[2]: align_hifi_reads.streamlit_panel(state)
with tabs[3]: haplotag_helper.streamlit_panel(state)
with tabs[4]: find_crossovers.streamlit_panel(state)
with tabs[5]: visualize_crossovers.streamlit_panel(state)
with tabs[6]:
    if compare_variant_vcfs:
        compare_variant_vcfs.streamlit_panel(state)
    else:
        st.info("compare_variant_vcfs.py not available.")

