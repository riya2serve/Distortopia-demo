#!/usr/bin/env python3
# scripts/visualize_crossovers.py
from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Optional smoothing (fallback if SciPy not present)
try:
    from scipy.ndimage import gaussian_filter1d  # type: ignore
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

sns.set_context("talk")
sns.set_style("whitegrid")


# -------------------- helpers --------------------
def _species_to_safe(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


def _smooth(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    if y.size == 0:
        return y
    if _HAVE_SCIPY and sigma_bins > 0:
        return gaussian_filter1d(y.astype(float), sigma=sigma_bins)
    # Fallback: 5-bin moving average
    w = 5
    if y.size < w:
        return y.astype(float)
    k = np.ones(w, dtype=float) / w
    return np.convolve(y.astype(float), k, mode="same")


def _one_breakpoint_positions(df: pd.DataFrame) -> np.ndarray:
    """
    Return one genomic position per read:
      - If both left & right present -> midpoint
      - Else whichever side is present
      - Else NaN (dropped)
    """
    l = pd.to_numeric(df["crossover_left"], errors="coerce") if "crossover_left" in df else pd.Series(dtype=float)
    r = pd.to_numeric(df["crossover_right"], errors="coerce") if "crossover_right" in df else pd.Series(dtype=float)

    # Align indices if needed
    if l.shape[0] != df.shape[0]: l = l.reindex(df.index)
    if r.shape[0] != df.shape[0]: r = r.reindex(df.index)

    pos = l.copy()
    both = l.notna() & r.notna()
    pos[both] = (l[both] + r[both]) / 2.0
    only_r = l.isna() & r.notna()
    pos[only_r] = r[only_r]

    return pos.dropna().to_numpy()


# -------------------- main plotting --------------------
def recombination_map(
    tsv_path: str | Path,
    species_name: str,
    outdir: str | Path = "plots",
    *,
    genome_bin_size: int = 1_000_000,  # 1 Mb
    smooth_sigma_bins: float = 2.0,
) -> Dict[str, str]:
    """
    Build a single smooth recombination map (crossovers per Mb) from a TSV.

    Saves one PNG: <species>_recombination_map.png
    Returns dict with summary and path.
    """
    tsv_path = str(tsv_path)
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    species_safe = _species_to_safe(species_name)

    df = pd.read_csv(tsv_path, sep="\t")
    if not (("crossover_left" in df.columns) or ("crossover_right" in df.columns)):
        raise ValueError("TSV must include 'crossover_left' and/or 'crossover_right' columns.")

    # One position per read (midpoint if both sides available)
    pos = _one_breakpoint_positions(df)
    if pos.size == 0:
        # Create a tiny blank plot to avoid crashing
        fig = plt.figure(figsize=(10, 4))
        plt.title(f"{species_name}: Recombination map (no crossovers)")
        plt.xlabel("Genomic position (Mb)")
        plt.ylabel("Crossovers per Mb")
        plt.tight_layout()
        out_png = outdir_p / f"{species_safe}_recombination_map.png"
        fig.savefig(out_png, dpi=300)
        plt.close(fig)
        return {"summary": f"{species_name}: 0 crossovers found", "recomb_map_png": str(out_png)}

    # Axis range from observed positions
    gmax = max(int(np.nanmax(pos)), genome_bin_size)
    bins = np.arange(0, gmax + genome_bin_size, genome_bin_size)
    counts, _ = np.histogram(pos, bins=bins)

    # Convert to rate per Mb
    denom_mb = genome_bin_size / 1e6
    rate = counts.astype(float) / denom_mb

    # Smooth in bin space
    rate_smooth = _smooth(rate, sigma_bins=smooth_sigma_bins)

    # Bin centers as Mb
    x = (bins[:-1] + genome_bin_size / 2.0) / 1e6

    # Plot single smooth map
    fig = plt.figure(figsize=(11, 4.8))
    plt.plot(x, rate_smooth, lw=2.5, color="slateblue")
    plt.fill_between(x, rate_smooth, color="slateblue", alpha=0.2)
    plt.title(f"{species_name}: Recombination map (smoothed crossovers per Mb)")
    plt.xlabel("Genomic position (Mb)")
    plt.ylabel("Crossovers per Mb")
    plt.tight_layout()

    out_png = outdir_p / f"{species_safe}_recombination_map.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    summary = f"{species_name}: n={pos.size} crossover positions → map written"
    return {"summary": summary, "recomb_map_png": str(out_png)}


# -------------------- streamlit panel --------------------
def streamlit_panel(state) -> None:
    import streamlit as st

    st.header("6) Recombination map (single smooth track)")

    c1, c2 = st.columns(2)
    with c1:
        th_tsv = Path(st.text_input(
            "Thaliana TSV", value=str(state.data / "A_thaliana_gamete_crossovers.tsv"), key="rec_th_tsv"))
        th_label = st.text_input("Thaliana label", value="A. thaliana", key="rec_th_label")
    with c2:
        ly_tsv = Path(st.text_input(
            "Lyrata TSV", value=str(state.data / "A_lyrata_gamete_crossovers.tsv"), key="rec_ly_tsv"))
        ly_label = st.text_input("Lyrata label", value="A. lyrata", key="rec_ly_label")

    outdir = Path(st.text_input("Output plots folder", value=str(state.plots), key="rec_outdir"))

    c3, c4 = st.columns(2)
    with c3:
        bin_size = int(st.number_input("Genome bin size (bp)", 50_000, 5_000_000, 1_000_000, 50_000, key="rec_bins"))
    with c4:
        sigma_bins = float(st.number_input("Smoothing (sigma in bins)", 0.0, 10.0, 2.0, 0.5, key="rec_sigma"))

    if st.button("Make recombination maps", type="primary", key="rec_go"):
        try:
            out_th = recombination_map(th_tsv, th_label, outdir, genome_bin_size=bin_size, smooth_sigma_bins=sigma_bins)
            out_ly = recombination_map(ly_tsv, ly_label, outdir, genome_bin_size=bin_size, smooth_sigma_bins=sigma_bins)

            st.success("Saved recombination maps to ./plots/")
            for tag, out in (("Thaliana", out_th), ("Lyrata", out_ly)):
                st.subheader(tag)
                st.caption(out["summary"])
                st.image(out["recomb_map_png"], use_container_width=True)
        except Exception as e:
            st.error(f"Visualization failed: {e}")


# -------------------- CLI --------------------
if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Create a single smooth recombination map from a crossover TSV.")
    ap.add_argument("--tsv", required=True, help="Input TSV from find_crossovers.py")
    ap.add_argument("--species", required=True, help="Species label")
    ap.add_argument("--outdir", default="plots", help="Output folder for PNG")
    ap.add_argument("--bin-size", type=int, default=1_000_000, help="Genome bin size (bp)")
    ap.add_argument("--sigma", type=float, default=2.0, help="Smoothing sigma in bins")
    args = ap.parse_args()

    outs = recombination_map(args.tsv, args.species, args.outdir,
                             genome_bin_size=args.bin_size, smooth_sigma_bins=args.sigma)
    print("\n".join(f"{k}: {v}" for k, v in outs.items()))

