#!/usr/bin/env python3
# scripts/visualize_crossovers.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Optional smoothing (fallback if SciPy not present)
try:
    from scipy.ndimage import gaussian_filter1d  # type: ignore
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

PURPLE = "#7B66B4"  # poster purple


# -------------------- helpers --------------------
def _species_to_safe(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


def _smooth(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    if y.size == 0:
        return y
    if _HAVE_SCIPY and sigma_bins > 0:
        return gaussian_filter1d(y.astype(float), sigma=sigma_bins)
    # Fallback: simple 5-bin moving average
    w = 5
    if y.size < w:
        return y.astype(float)
    k = np.ones(w, dtype=float) / w
    return np.convolve(y.astype(float), k, mode="same")


def _chrom_sort_key(ch: str) -> tuple:
    import re
    s = str(ch).replace("chr", "").replace("Chr", "")
    m = re.search(r"(\d+)$", s)
    return (0, int(m.group(1))) if m else (1, str(ch))


def _xo_midpoints(df: pd.DataFrame) -> np.ndarray:
    """Return one XO position per recombinant read (midpoint if both sides)."""
    l = pd.to_numeric(df["crossover_left"], errors="coerce")
    r = pd.to_numeric(df["crossover_right"], errors="coerce")
    pos = l.copy()
    both = l.notna() & r.notna()
    pos[both] = (l[both] + r[both]) / 2.0
    only_r = l.isna() & r.notna()
    pos[only_r] = r[only_r]
    return pos.dropna().to_numpy()


def _bin_edges(gmax: int, bin_size: int) -> np.ndarray:
    gmax = max(int(gmax), bin_size)
    return np.arange(0, gmax + bin_size, bin_size)


def _coverage_per_bin_from_intervals(starts: np.ndarray, ends: np.ndarray, bins: np.ndarray) -> np.ndarray:
    """
    Coverage per bin using a difference-array sweep:
      - starts, ends: 1-based inclusive from TSV
      - bins: 0-based edges [0, B, 2B, ...]
    """
    # Convert to 0-based half-open intervals
    s0 = np.maximum(starts.astype(np.int64) - 1, 0)
    e0 = ends.astype(np.int64)

    # Map to bin indices
    li = np.clip(np.searchsorted(bins, s0, side="right") - 1, 0, len(bins) - 2)
    ri = np.clip(np.searchsorted(bins, e0, side="left") - 1, 0, len(bins) - 2)

    cov = np.zeros(len(bins) - 1, dtype=np.int64)
    np.add.at(cov, li, 1)
    # decrement after rightmost bin touched
    dec_idx = np.minimum(ri + 1, len(cov) - 1)
    np.add.at(cov, dec_idx, -1)
    cov = np.cumsum(cov)
    cov[cov < 0] = 0
    return cov


# -------------------- per-chromosome RATE maps --------------------
def recombination_maps_by_chrom_rate(
    tsv_path: str | Path,
    species_name: str,
    outdir: str | Path = "plots",
    *,
    genome_bin_size: int = 1_000_000,    # 1 Mb
    smooth_sigma_bins: float = 2.0,
    chrom_lengths: Optional[Dict[str, int]] = None,  # optional real lengths for nice axes
) -> Dict[str, str]:
    """
    Create **per-chromosome recombination rate** maps (one PNG per chromosome).
    Rate(bin) = XO_count_in_bin / informative_reads_overlapping_bin  (per-Mb bins)
    Returns a dict of {chrom: path_to_png}.
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    species_safe = _species_to_safe(species_name)

    df = pd.read_csv(tsv_path, sep="\t")

    # Mask for recombinant reads (numerator)
    xo_mask = (df["crossover_left"].astype(str) != "NA") | (df["crossover_right"].astype(str) != "NA")

    chroms = sorted(df["scaff"].astype(str).unique().tolist(), key=_chrom_sort_key)
    pngs: Dict[str, str] = {}

    for chrom in chroms:
        sub = df[df["scaff"].astype(str) == chrom]
        if sub.empty:
            continue

        # Bin edges
        if chrom_lengths and chrom in chrom_lengths:
            gmax = int(chrom_lengths[chrom])
        else:
            gmax = int(np.nanmax(sub["end"]))
        bins = _bin_edges(gmax, genome_bin_size)

        # Denominator per bin: informative read coverage (XO + non-XO)
        coverage = _coverage_per_bin_from_intervals(
            sub["start"].to_numpy(), sub["end"].to_numpy(), bins
        )

        # Numerator per bin: XO placements (midpoints)
        xo_sub = sub[xo_mask & (df["scaff"].astype(str) == chrom)]
        xo_pos = _xo_midpoints(xo_sub) if not xo_sub.empty else np.array([], dtype=float)

        xo_counts = np.zeros(len(bins) - 1, dtype=np.int64)
        if xo_pos.size:
            idx = np.searchsorted(bins, xo_pos, side="right") - 1
            idx = idx[(idx >= 0) & (idx < len(xo_counts))]
            if idx.size:
                xo_counts = np.bincount(idx, minlength=len(xo_counts))

        # Rate per bin (fraction); avoid divide-by-zero
        denom = coverage.astype(float)
        denom[denom == 0.0] = np.nan
        rate = xo_counts.astype(float) / denom
        rate = np.nan_to_num(rate, nan=0.0)

        # Smooth and plot
        rate_smooth = _smooth(rate, sigma_bins=smooth_sigma_bins)
        x_mb = (bins[:-1] + genome_bin_size / 2.0) / 1e6

        fig = plt.figure(figsize=(11, 3.5))
        plt.plot(x_mb, rate_smooth, lw=2.2, color=PURPLE)
        plt.fill_between(x_mb, rate_smooth, color=PURPLE, alpha=0.25)
        plt.title(f"{species_name} — {chrom}: recombination rate (XO / informative reads per Mb)")
        plt.xlabel("Position (Mb)")
        plt.ylabel("Rate/Mb")
        plt.tight_layout()

        out_png = outdir / f"{species_safe}_{chrom}_recomb_rate_perMb.png"
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

        pngs[chrom] = str(out_png)

    return pngs


# -------------------- Streamlit panel (optional) --------------------
def streamlit_panel(state) -> None:
    import streamlit as st
    st.header("6) Recombination rate maps (per chromosome)")

    tsv = Path(st.text_input("Crossover TSV", value=str(state.data / "Athal_gamete_CO.tsv")))
    label = st.text_input("Species label", value="A. thaliana")
    outdir = Path(st.text_input("Output plots folder", value=str(state.plots)))

    c1, c2 = st.columns(2)
    with c1:
        bin_size = int(st.number_input("Bin size (bp)", 50_000, 5_000_000, 1_000_000, 50_000))
    with c2:
        sigma_bins = float(st.number_input("Smoothing sigma (bins)", 0.0, 10.0, 2.0, 0.5))

    with st.expander("Chromosome lengths (optional)"):
        chr_lengths = {}
        for ch in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
            val = st.text_input(f"{ch} length (bp)", value="", key=f"len_{ch}")
            if val.strip():
                try:
                    chr_lengths[ch] = int(val.replace(",", ""))
                except:
                    st.warning(f"Ignoring invalid length for {ch}")

    if st.button("Make per-chromosome rate maps", type="primary"):
        try:
            outs = recombination_maps_by_chrom_rate(
                tsv, label, outdir,
                genome_bin_size=bin_size,
                smooth_sigma_bins=sigma_bins,
                chrom_lengths=chr_lengths if chr_lengths else None,
            )
            st.success(f"Saved {len(outs)} PNGs to {outdir}")
            # Optionally show first few thumbnails
            for ch in sorted(outs, key=_chrom_sort_key)[:3]:
                st.image(outs[ch], caption=ch, use_container_width=True)
        except Exception as e:
            st.error(f"Failed: {e}")


# -------------------- CLI --------------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Make per-chromosome recombination RATE maps (XO / informative reads per Mb) from a TSV."
    )
    ap.add_argument("--tsv", required=True, help="Input TSV from find_crossovers.py")
    ap.add_argument("--species", required=True, help="Species label")
    ap.add_argument("--outdir", default="plots", help="Output folder for PNGs")
    ap.add_argument("--bin-size", type=int, default=1_000_000, help="Genome bin size (bp)")
    ap.add_argument("--sigma", type=float, default=2.0, help="Smoothing sigma (in bins)")
    args = ap.parse_args()

    outs = recombination_maps_by_chrom_rate(
        args.tsv, args.species, args.outdir,
        genome_bin_size=args.bin_size, smooth_sigma_bins=args.sigma
    )
    print("\n".join(f"{k}: {v}" for k, v in outs.items()))
