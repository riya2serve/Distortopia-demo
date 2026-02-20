#!/usr/bin/env python
"""
Plot crossover (recombination) frequency along chromosomes from disto TSVs.

Supports TWO input types:

1) Read-level TSV (from infer.py) columns include:
    scaff, start, end, crossover_left, crossover_right, read
   - Uses crossover midpoint = (crossover_left + crossover_right)/2
   - End-masks using read alignment span columns (start/end) with read_length bp
   - Normalizes by unique reads per chromosome (or nreads_total if provided)
   - Y-axis: crossover frequency (per read per bin)

2) Locus-level TSV (from disto collapse) columns include:
    scaff, locus_left, locus_right
   Optional: cx_mid, n_reads
   - Uses cx_mid if present else midpoint(locus_left/locus_right)
   - End-masks using cx_mid relative to inferred chrom end (max locus_right) with read_length bp
   - Histograms are WEIGHTED by n_reads if present
   - Normalizes by total support (sum n_reads) if present, else by number of loci
   - Y-axis: crossover frequency (weighted support per bin)

Writes ONE PDF PER CHROMOSOME:
    <prefix>.<chrom>.pdf
If you pass chroms=[...], it will only write those.

Optional extras:
- errorbars=True + bin_bp=<int>: add SD error bars computed across reads per bin (read-mode only),
  using the same logic as plot_errorbars.py (fixed bp binning, per-read per-bin counts).
- rolling_stats=True: overlay rolling mean and rolling median curves (centered window in bins).

Color scheme (more distinguishable):
- Bars (hist): light gray
- Raw line (same as bars): medium gray
- SD error bars: dark gray
- Rolling mean: blue
- Rolling median: orange/red
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Literal, Tuple

import numpy as np
import pandas as pd
import toyplot
import toyplot.pdf
from loguru import logger


READ_REQUIRED = {"scaff", "start", "end", "crossover_left", "crossover_right", "read"}
LOCUS_REQUIRED = {"scaff", "locus_left", "locus_right"}


def _to_numeric(series: pd.Series) -> pd.Series:
    """Convert series that may contain 'NA' strings to numeric; invalid -> NaN."""
    return pd.to_numeric(series.replace("NA", np.nan), errors="coerce")


def _detect_mode(cols: set[str]) -> Literal["read", "locus"]:
    """Detect whether TSV is read-level (infer) or locus-level (collapse)."""
    if READ_REQUIRED.issubset(cols):
        return "read"
    if LOCUS_REQUIRED.issubset(cols):
        return "locus"
    raise ValueError(
        "TSV is neither read-level (infer) nor locus-level (collapse).\n"
        f"Read-level requires: {sorted(READ_REQUIRED)}\n"
        f"Locus-level requires: {sorted(LOCUS_REQUIRED)}\n"
        f"Found columns: {sorted(cols)}"
    )


def _hist_mean_sd_by_bin_readmode(
    sub: pd.DataFrame,
    *,
    bin_bp: int,
    n_reads: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    EXACT SD logic from plot_errorbars.py (read-mode only, fixed bp binning).

    For each bin:
      - For each read, k = # crossover calls that read contributes to that bin
      - mean = E[k] across reads
      - sd   = sqrt(E[k^2] - (E[k])^2) across reads
    Returns: (mean, sd, bin_edges)
    """
    cxs = sub["cx_mid"].to_numpy()

    left = int((cxs.min() // bin_bp) * bin_bp)
    right = int(((cxs.max() // bin_bp) + 1) * bin_bp)
    bin_edges = np.arange(left, right + bin_bp, bin_bp, dtype=int)

    nbins = len(bin_edges) - 1
    if nbins <= 0:
        raise ValueError("No bins created — check bin_bp and data range.")

    # assign each crossover event to a bin index
    bidx = np.digitize(cxs, bin_edges, right=False) - 1
    bidx = np.clip(bidx, 0, nbins - 1)

    tmp = pd.DataFrame({"bin": bidx, "read": sub["read"].astype(str).to_numpy()})

    # per-read per-bin counts (k)
    per = tmp.groupby(["bin", "read"]).size().rename("k").reset_index()

    # sum(k) and sum(k^2) per bin (missing bins => 0)
    sum_k = (
        per.groupby("bin")["k"]
        .sum()
        .reindex(range(nbins), fill_value=0)
        .to_numpy(float)
    )
    sum_k2 = (
        per.groupby("bin")["k"]
        .apply(lambda s: (s * s).sum())
        .reindex(range(nbins), fill_value=0)
        .to_numpy(float)
    )

    mean = sum_k / float(n_reads)
    ex2 = sum_k2 / float(n_reads)
    var = ex2 - mean * mean
    var[var < 0] = 0.0  # numerical safety
    sd = np.sqrt(var)

    return mean, sd, bin_edges


def _rolling_stats(y: np.ndarray, window_bins: int) -> Tuple[np.ndarray, np.ndarray]:
    """Centered rolling mean and median with min_periods=1."""
    if window_bins <= 1:
        return y.copy(), y.copy()

    s = pd.Series(y)
    r = s.rolling(window=window_bins, center=True, min_periods=1)
    return r.mean().to_numpy(), r.median().to_numpy()


def run_plot(
    tsv: Path,
    outdir: Path,
    prefix: str,
    chroms: Optional[List[str]] = None,
    bins: int = 50,
    read_length: int = 100_000,
    nreads_total: Optional[int] = None,
    *,
    # NEW
    errorbars: bool = False,
    bin_bp: Optional[int] = None,
    rolling_stats: bool = False,
    rolling_window_bins: int = 7,
) -> None:
    """
    Parameters
    ----------
    errorbars:
        If True, add SD error bars (read-mode only). Requires bin_bp.
    bin_bp:
        Fixed bin size in bp for errorbar mode. If errorbars=True and bin_bp is None -> error.
    rolling_stats:
        If True, overlay rolling mean + rolling median lines.
    rolling_window_bins:
        Window size (in number of bins) for rolling mean/median (centered).
    """
    tsv = Path(tsv).expanduser().resolve()
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(tsv, sep="\t")
    if data.empty:
        raise ValueError("Input TSV is empty. Nothing to plot.")

    mode = _detect_mode(set(data.columns))
    logger.info(f"Detected plot mode: {mode}")

    # Normalize scaff column to string
    data["scaff"] = data["scaff"].astype(str)

    # Optional filter by chromosome(s)
    if chroms:
        chrom_set = set(map(str, chroms))
        data = data[data["scaff"].isin(chrom_set)].copy()
        if data.empty:
            raise ValueError(f"No rows found for chrom(s): {chroms}. Check contig names in TSV.")

    # ---- Mode-specific normalization / midpoint / weights ----
    weights_col: Optional[str] = None
    ylab: str
    title_suffix: str

    if mode == "read":
        # Convert numeric-like columns
        data["start"] = pd.to_numeric(data["start"], errors="coerce")
        data["end"] = pd.to_numeric(data["end"], errors="coerce")
        data["crossover_left"] = _to_numeric(data["crossover_left"])
        data["crossover_right"] = _to_numeric(data["crossover_right"])

        # Keep only rows with inferred crossovers
        data = data.dropna(subset=["start", "end", "crossover_left", "crossover_right"])
        if data.empty:
            raise ValueError("No crossover rows found (all crossover_left/right are NA). Nothing to plot.")

        data["cx_mid"] = ((data["crossover_left"] + data["crossover_right"]) / 2.0).astype(int)

        # No weights in read-level mode (each read contributes 1)
        weights_col = None
        ylab = "Crossover frequency (per read per bin)"
        title_suffix = "reads"

    else:
        # locus mode
        if errorbars:
            raise ValueError("errorbars=True is only supported for read-level (infer) TSVs.")

        data["locus_left"] = pd.to_numeric(data["locus_left"], errors="coerce")
        data["locus_right"] = pd.to_numeric(data["locus_right"], errors="coerce")
        data = data.dropna(subset=["locus_left", "locus_right"]).copy()
        if data.empty:
            raise ValueError("No locus rows found (locus_left/right missing). Nothing to plot.")

        data["locus_left"] = data["locus_left"].astype(int)
        data["locus_right"] = data["locus_right"].astype(int)

        if "cx_mid" in data.columns:
            data["cx_mid"] = pd.to_numeric(data["cx_mid"], errors="coerce")
            data = data.dropna(subset=["cx_mid"]).copy()
            data["cx_mid"] = data["cx_mid"].astype(int)
        else:
            data["cx_mid"] = ((data["locus_left"] + data["locus_right"]) / 2.0).astype(int)

        # Optional weighting by read support
        if "n_reads" in data.columns:
            data["n_reads"] = pd.to_numeric(data["n_reads"], errors="coerce").fillna(1).astype(int)
            data.loc[data["n_reads"] < 1, "n_reads"] = 1
            weights_col = "n_reads"
            ylab = "Crossover frequency (weighted support per bin)"
            title_suffix = "weighted support"
        else:
            weights_col = None
            ylab = "Crossover frequency (per locus per bin)"
            title_suffix = "loci"

    chrom_list = list(data["scaff"].unique())
    wrote_any = False

    for chrom in chrom_list:
        sub = data.loc[data["scaff"] == chrom].copy()
        if sub.empty:
            continue

        # ---- End-masking ----
        if mode == "read":
            sub = sub.loc[sub["start"] > read_length]
            chrom_end = sub["end"].max()
            sub = sub.loc[sub["end"] < (chrom_end - read_length)]
        else:
            chrom_end = sub["locus_right"].max()
            sub = sub.loc[sub["cx_mid"] > read_length]
            sub = sub.loc[sub["cx_mid"] < (chrom_end - read_length)]

        if sub.empty:
            logger.warning(f"{chrom}: no rows remain after end-masking; skipping.")
            continue

        # ---- Normalization ----
        if mode == "read":
            if nreads_total is None:
                n_norm = int(sub["read"].nunique())
                if n_norm <= 0:
                    logger.warning(f"{chrom}: zero unique reads after filtering; skipping.")
                    continue
            else:
                n_norm = max(1, int(nreads_total))
        else:
            if weights_col == "n_reads":
                n_norm = int(sub["n_reads"].sum())
            else:
                n_norm = int(len(sub))
            n_norm = max(1, n_norm)

        # ---- Histogram (+ optional SD for read-mode) ----
        sd = None

        if mode == "read" and errorbars:
            if bin_bp is None:
                raise ValueError("errorbars=True requires bin_bp (fixed bp bins).")

            # SD is "across reads", so use actual unique reads in this chrom after masking
            n_reads = int(sub["read"].nunique())
            if n_reads <= 0:
                logger.warning(f"{chrom}: zero unique reads after filtering; skipping.")
                continue

            mags, sd, bin_edges = _hist_mean_sd_by_bin_readmode(
                sub,
                bin_bp=int(bin_bp),
                n_reads=n_reads,
            )
            n_norm = n_reads  # title shows actual read count used for SD

        else:
            cxs = sub["cx_mid"].to_numpy()
            w = sub[weights_col].to_numpy() if weights_col else None
            mags, bin_edges = np.histogram(cxs, bins=bins, weights=w)
            mags = mags / float(n_norm)

        # Convert x-axis to Mb
        bin_edges_mb = bin_edges / 1e6
        mids_mb = (bin_edges_mb[:-1] + bin_edges_mb[1:]) / 2.0

        # Rolling stats (mean/median) over bins
        roll_mean = roll_med = None
        if rolling_stats:
            roll_mean, roll_med = _rolling_stats(mags, window_bins=int(rolling_window_bins))

        # Canvas + axes
        canvas = toyplot.Canvas(width=600, height=350)
        axes = canvas.cartesian(
            xlabel="Chromosome position (Mb)",
            ylabel=ylab,
            margin=60,
        )

        # Title
        title_extra = ""
        if mode == "read" and errorbars:
            title_extra = f", bin={int(bin_bp)//1000}kb"
        elif rolling_stats:
            title_extra = f", roll={int(rolling_window_bins)} bins"

        canvas.text(
            300,
            20,
            f"{prefix} — {chrom} (N={n_norm} {title_suffix}){title_extra}",
            style={"font-size": "16px", "font-weight": "bold", "text-anchor": "middle"},
        )

        # Styling
        axes.x.ticks.show = True
        axes.x.ticks.far = 0
        axes.x.ticks.near = 5
        axes.x.ticks.labels.offset = 10
        axes.x.label.offset = 28
        axes.x.label.style["font-size"] = 14

        axes.y.ticks.show = True
        axes.y.ticks.far = 0
        axes.y.ticks.near = 5
        axes.y.ticks.labels.offset = 10
        axes.y.label.offset = 28
        axes.y.label.style["font-size"] = 14

        # -----------------------------
        # PLOTTING (cleaner colors)
        # -----------------------------

        # Bars: light gray
        m0 = axes.bars((mags, bin_edges_mb), opacity=0.85)
        m0._fill.color = "#d9d9d9"

        # Raw line: medium gray (optional but helps guide the eye)
        m_raw = axes.plot(mids_mb, mags, opacity=0.9)
        m_raw._stroke.color = "#7a7a7a"
        m_raw._stroke.width = 1.5

        # Error bars (SD): dark gray vertical segments at each bin midpoint
        if sd is not None:
            for x, y, e in zip(mids_mb, mags, sd):
                y0 = max(0.0, y - e)
                y1 = y + e
                seg = axes.plot([x, x], [y0, y1], opacity=0.85)
                seg._stroke.color = "#333333"
                seg._stroke.width = 1.0

        # Rolling mean / median overlays
        if roll_mean is not None and roll_med is not None:
            rm = axes.plot(mids_mb, roll_mean, opacity=1.0)
            rd = axes.plot(mids_mb, roll_med, opacity=1.0)

            # Rolling mean: blue
            rm._stroke.color = "#1f77b4"
            rm._stroke.width = 2.5

            # Rolling median: orange/red
            rd._stroke.color = "#d62728"
            rd._stroke.width = 2.5

        out = outdir / f"{prefix}.{chrom}.pdf"
        toyplot.pdf.render(canvas, str(out))
        logger.info(f"wrote {out}")
        wrote_any = True

    if not wrote_any:
        raise ValueError("No plots generated (all chromosomes empty after filtering/masking).")


if __name__ == "__main__":
    raise SystemExit("Run via the CLI, e.g. `disto plot --tsv ... --out ... --prefix ...`")
