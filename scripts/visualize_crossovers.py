#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st

sns.set_context("talk")


# ---------- core plotting ----------
def plot_crossovers(tsv_path: str, species_name: str, outdir: str) -> dict:
    df = pd.read_csv(tsv_path, sep="\t")
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    # Clean numeric columns (turn "NA" into NaN)
    df["crossover_left"] = pd.to_numeric(df["crossover_left"], errors="coerce")
    df["crossover_right"] = pd.to_numeric(df["crossover_right"], errors="coerce")

    # Helper name for filenames
    species_safe = species_name.replace(" ", "_").replace(".", "")

    # 1) Crossovers per scaffold
    # Count a crossover if either left or right is present on the read
    has_any = df["crossover_left"].notna() | df["crossover_right"].notna()
    counts = (
        df.assign(has_any=has_any)
          .groupby("scaff", as_index=False)["has_any"]
          .sum(numeric_only=True)
          .rename(columns={"has_any": "crossovers"})
    )

    fig1 = plt.figure(figsize=(9, 5))
    sns.barplot(x="scaff", y="crossovers", data=counts, color="steelblue")
    plt.title(f"{species_name}: Crossovers per scaffold")
    plt.xlabel("Scaffold")
    plt.ylabel("Crossover count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    out1 = outdir_p / f"{species_safe}_crossovers_per_scaffold.png"
    fig1.savefig(out1, dpi=300)
    plt.close(fig1)

    # 2) Breakpoint density — overlay left vs right (Option B)
    left_vals = df["crossover_left"].dropna()
    right_vals = df["crossover_right"].dropna()

    fig2 = plt.figure(figsize=(10, 5))
    sns.histplot(left_vals, bins=60, alpha=0.6, label="Left breakpoints", color="steelblue")
    sns.histplot(right_vals, bins=60, alpha=0.6, label="Right breakpoints", color="tomato")
    plt.title(f"{species_name}: Breakpoint distribution (left vs right)")
    plt.xlabel("Genomic position (bp)")
    plt.ylabel("Count")
    plt.legend(title="Breakpoint side")
    plt.tight_layout()
    out2 = outdir_p / f"{species_safe}_breakpoint_density_overlay.png"
    fig2.savefig(out2, dpi=300)
    plt.close(fig2)

    # 3) Small summary
    n_reads = len(df)
    n_with_co = int(has_any.sum())
    summary = f"{species_name}: {n_with_co}/{n_reads} reads ({n_with_co/n_reads:.2%}) show ≥1 crossover"

    return {
        "summary": summary,
        "per_scaffold_png": str(out1),
        "density_png": str(out2),
        "counts_df": counts,
    }


# ---------- programmatic entry ----------
def run(tsv_path: Path, species_name: str, outdir: Path) -> dict:
    return plot_crossovers(str(tsv_path), species_name, str(outdir))


# ---------- Streamlit panel ----------
def streamlit_panel(state):
    st.header("6) Visualize crossover TSVs")

    col1, col2 = st.columns(2)
    with col1:
        th_tsv = Path(
            st.text_input(
                "Thaliana TSV",
                value=str(state.data / "A_thaliana_gamete_crossovers.tsv"),
                key="viz_th_tsv",
            )
        )
        th_label = st.text_input("Thaliana label", value="A. thaliana", key="viz_th_label")

    with col2:
        ly_tsv = Path(
            st.text_input(
                "Lyrata TSV",
                value=str(state.data / "A_lyrata_gamete_crossovers.tsv"),
                key="viz_ly_tsv",
            )
        )
        ly_label = st.text_input("Lyrata label", value="A. lyrata", key="viz_ly_label")

    outdir = Path(
        st.text_input("Output plots folder", value=str(state.plots), key="viz_outdir")
    )

    if st.button("Generate plots", key="viz_go"):
        try:
            out_th = run(th_tsv, th_label, outdir)
            out_ly = run(ly_tsv, ly_label, outdir)

            st.success(f"Saved plots to: {outdir.resolve()}")
            for out in (out_th, out_ly):
                st.caption(out["summary"])
                st.image(out["per_scaffold_png"], caption="Crossovers per scaffold")
                st.image(out["density_png"], caption="Breakpoint density (left vs right)")
        except Exception as e:
            st.error(str(e))

