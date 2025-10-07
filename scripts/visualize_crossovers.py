#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
sns.set_context("talk")

def plot_crossovers(tsv_path, species_name, outdir="plots"):
    df = pd.read_csv(tsv_path, sep="\t")
    outdir = Path(outdir); outdir.mkdir(exist_ok=True)
    df["crossover_left"]  = pd.to_numeric(df["crossover_left"], errors="coerce")
    df["crossover_right"] = pd.to_numeric(df["crossover_right"], errors="coerce")
    counts = (df.groupby("scaff")["crossover_left"].apply(lambda x: x.notna().sum()).reset_index())
    fig1 = plt.figure(figsize=(8,5)); sns.barplot(x="scaff", y="crossover_left", data=counts)
    plt.title(f"{species_name}: Crossovers per scaffold"); plt.xticks(rotation=45); plt.tight_layout()
    out1 = outdir / f"{species_name}_crossovers_per_scaffold.png"; fig1.savefig(out1, dpi=300); plt.close(fig1)
    fig2 = plt.figure(figsize=(10,5)); sns.histplot(df["crossover_left"].dropna(), bins=50)
    plt.title(f"{species_name}: Breakpoint distribution"); plt.tight_layout()
    out2 = outdir / f"{species_name}_crossover_density.png"; fig2.savefig(out2, dpi=300); plt.close(fig2)
    df["has_CO"] = df["crossover_left"].notna(); n_co = int(df["has_CO"].sum())
    return {"summary": f"{species_name}: {n_co}/{len(df)} reads ({n_co/len(df):.2%}) with CO",
            "per_scaffold_png": str(out1), "density_png": str(out2), "counts_df": counts}

def run(tsv_path: Path, species_name: str, outdir: Path) -> dict:
    return plot_crossovers(str(tsv_path), species_name, str(outdir))

def streamlit_panel(state):
    st.header("6) Visualize crossover TSVs")
    col1,col2 = st.columns(2)
    with col1:
        th_tsv = Path(st.text_input("Thaliana TSV", value=str(state.data / "A_thaliana_gamete_crossovers.tsv")))
        th_label = st.text_input("Thaliana label", value="A. thaliana")
    with col2:
        ly_tsv = Path(st.text_input("Lyrata TSV", value=str(state.data / "A_lyrata_gamete_crossovers.tsv")))
        ly_label = st.text_input("Lyrata label", value="A. lyrata")
    outdir = Path(st.text_input("Output plots folder", value=str(state.plots)))
    if st.button("Generate plots"):
        try:
            out_th = run(th_tsv, th_label, outdir); out_ly = run(ly_tsv, ly_label, outdir)
            st.success("Saved plots to ./plots/")
            for out in (out_th, out_ly):
                st.caption(out["summary"]); st.image(out["per_scaffold_png"]); st.image(out["density_png"])
        except Exception as e:
            st.error(str(e))

