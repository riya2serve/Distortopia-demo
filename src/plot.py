#!/usr/bin/env python
"""
Plot crossover (recombination) frequency along chromosomes from disto infer TSV.

Expected TSV columns (from infer.py):
    scaff, start, end, nsnps, phased_snps, crossover_left, crossover_right, read

This plot:
- Uses crossover midpoint = (crossover_left + crossover_right)/2 for rows with a crossover.
- Masks ends based on read alignment start/end columns (read_length bp from each end).
- Writes ONE PDF PER CHROMOSOME:
    <prefix>.<chrom>.pdf
If you pass chroms=[...], it will only write those.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import toyplot
import toyplot.pdf
from loguru import logger


REQUIRED_COLS = {"scaff", "start", "end", "crossover_left", "crossover_right", "read"}


def _to_numeric(series: pd.Series) -> pd.Series:
    """Convert series that may contain 'NA' strings to numeric; invalid -> NaN."""
    return pd.to_numeric(series.replace("NA", np.nan), errors="coerce")


def run_plot(
    tsv: Path,
    outdir: Path,
    prefix: str,
    chroms: Optional[List[str]] = None,
    bins: int = 50,
    read_length: int = 100_000,
    nreads_total: Optional[int] = None,
) -> None:
    tsv = Path(tsv).expanduser().resolve()
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # infer.py writes TAB-separated TSV
    data = pd.read_csv(tsv, sep="\t")

    missing = REQUIRED_COLS - set(data.columns)
    if missing:
        raise ValueError(f"TSV missing required columns: {sorted(missing)}")

    # Normalize scaff column to string
    data["scaff"] = data["scaff"].astype(str)

    # Optional filter by chromosome(s)
    if chroms:
        chrom_set = set(map(str, chroms))
        data = data[data["scaff"].isin(chrom_set)].copy()
        if data.empty:
            raise ValueError(f"No rows found for chrom(s): {chroms}. Check contig names in TSV.")

    # Convert numeric-like columns
    data["start"] = pd.to_numeric(data["start"], errors="coerce")
    data["end"] = pd.to_numeric(data["end"], errors="coerce")
    data["crossover_left"] = _to_numeric(data["crossover_left"])
    data["crossover_right"] = _to_numeric(data["crossover_right"])

    # Keep only rows with inferred crossovers
    data = data.dropna(subset=["start", "end", "crossover_left", "crossover_right"])
    if data.empty:
        raise ValueError("No crossover rows found (all crossover_left/right are NA). Nothing to plot.")

    # Midpoint position for each crossover (bp)
    data["cx_mid"] = ((data["crossover_left"] + data["crossover_right"]) / 2.0).astype(int)

    chrom_list = list(data["scaff"].unique())

    wrote_any = False

    for chrom in chrom_list:
        sub = data.loc[data["scaff"] == chrom].copy()
        if sub.empty:
            continue

        # End-masking based on read span columns
        sub = sub.loc[sub["start"] > read_length]
        chrom_end = sub["end"].max()
        sub = sub.loc[sub["end"] < (chrom_end - read_length)]

        if sub.empty:
            logger.warning(f"{chrom}: no crossover rows remain after end-masking; skipping.")
            continue

        # Normalization (default: unique reads for this chrom)
        if nreads_total is None:
            n_norm = int(sub["read"].nunique())
            if n_norm <= 0:
                logger.warning(f"{chrom}: zero unique reads after filtering; skipping.")
                continue
        else:
            n_norm = max(1, int(nreads_total))

        cxs = sub["cx_mid"].to_numpy()
        mags, bin_edges = np.histogram(cxs, bins=bins)

        # Normalize: frequency per read per bin
        mags = mags / float(n_norm)

        # Convert x-axis to Mb
        bin_edges_mb = bin_edges / 1e6
        mids_mb = (bin_edges_mb[:-1] + bin_edges_mb[1:]) / 2.0

        # Canvas + axes
        canvas = toyplot.Canvas(width=600, height=350)
        axes = canvas.cartesian(
            xlabel="Chromosome position (Mb)",
            ylabel="Crossover frequency (per read per bin)",
            margin=60,
        )

        # Title (Toyplot: use canvas.text)
        canvas.text(
            300,
            20,
            f"{prefix} — {chrom} (N={n_norm} reads)",
            style={"font-size": "16px", "font-weight": "bold", "text-anchor": "middle"},
        )

        # Styling (kept close to your original)
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

        # Bars + line
        m0 = axes.bars((mags, bin_edges_mb), opacity=0.6)
        m1 = axes.plot(mids_mb, mags, opacity=0.8)

        # Recolor (your palette)
        m0._fill.color = "#d1cef6"
        m1._stroke.color = "#7b6de2"

        out = outdir / f"{prefix}.{chrom}.pdf"
        toyplot.pdf.render(canvas, str(out))
        logger.info(f"wrote {out}")
        wrote_any = True

    if not wrote_any:
        raise ValueError("No plots generated (all chromosomes empty after filtering/masking).")


if __name__ == "__main__":
    raise SystemExit("Run via the CLI, e.g. `disto plot --tsv ... --out ... --prefix ...`")

