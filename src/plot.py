#!/usr/bin/env python
"""
Plot recombination (crossover) frequency along chromosomes from disto infer TSV.

Input TSV columns expected (from infer.py):
    scaff, start, end, nsnps, phased_snps, crossover_left, crossover_right, read

This plot uses crossover midpoint = (crossover_left + crossover_right)/2 for rows
where crossover_left/right are defined (not NA).

Outputs:
- By default: a multi-page PDF with one page per chromosome: <prefix>.pdf
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import toyplot
import toyplot.pdf
from loguru import logger


def _to_numeric_series(s: pd.Series) -> pd.Series:
    """Convert a series that may contain 'NA' strings to numeric (float) with NaN for invalid."""
    return pd.to_numeric(s.replace("NA", np.nan), errors="coerce")


def run_plot(
    tsv: Path,
    outdir: Path,
    prefix: str,
    chroms: Optional[List[str]] = None,
    bins: int = 100,
    read_length: int = 100_000,
    nreads_total: Optional[int] = None,
    one_file_per_chrom: bool = False,
):
    """
    Plot crossover frequency along chromosomes.

    Args:
        tsv: Path to infer output TSV.
        outdir: Output directory.
        prefix: Output prefix.
        chroms: Optional list of scaffolds to plot (if None, plot all found).
        bins: Histogram bins.
        read_length: End-masking in bp (exclude reads/windows within this distance of ends).
        nreads_total: Optional normalization constant. If None, uses number of unique reads per chrom.
        one_file_per_chrom: If True, write <prefix>.<chrom>.pdf per chrom; else one multi-page <prefix>.pdf.
    """
    tsv = Path(tsv).expanduser().resolve()
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Read TSV from infer.py (tab-separated)
    data = pd.read_csv(tsv, sep="\t")

    # Prelim alidation
    required = {"scaff", "start", "end", "crossover_left", "crossover_right", "read"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"TSV missing required columns: {sorted(missing)}")

    # Optional chromosome filter
    if chroms:
        chrom_set = set(chroms)
        data = data[data["scaff"].isin(chrom_set)].copy()
        if data.empty:
            raise ValueError(f"No rows found for chrom(s): {chroms}. Check contig names in TSV.")

    # Convert numeric columns
    data["start"] = pd.to_numeric(data["start"], errors="coerce")
    data["end"] = pd.to_numeric(data["end"], errors="coerce")
    data["crossover_left"] = _to_numeric_series(data["crossover_left"])
    data["crossover_right"] = _to_numeric_series(data["crossover_right"])

    # Keep only rows with inferred crossover
    data = data.dropna(subset=["crossover_left", "crossover_right", "start", "end"])
    if data.empty:
        raise ValueError("No crossover rows found in TSV (all crossover_left/right are NA). Nothing to plot.")

    # Compute crossover midpoint (bp)
    data["cx_mid"] = ((data["crossover_left"] + data["crossover_right"]) / 2.0).astype(int)

    # Sort chromosomes for stable output
    chrom_list = list(data["scaff"].unique())

    # Output handling
    if one_file_per_chrom:
        pdf_out = None
    else:
        pdf_out = outdir / f"{prefix}.pdf"
        # We'll accumulate pages and write once at end
        pages = []

    for chrom in chrom_list:
        sub = data.loc[data["scaff"] == chrom].copy()
        if sub.empty:
            continue

        # Mask sites in the first and last read_length regions using read alignment span
        sub = sub.loc[sub["start"] > read_length]
        chrom_end = sub["end"].max()
        sub = sub.loc[sub["end"] < (chrom_end - read_length)]

        if sub.empty:
            logger.warning(f"{chrom}: no crossovers remain after end-masking; skipping.")
            continue

        # Determine normalization
        if nreads_total is None:
            # Default than 1e7: normalize by unique reads that contributed on this chrom
            n_norm = sub["read"].nunique()
            if n_norm == 0:
                logger.warning(f"{chrom}: zero unique reads; skipping.")
                continue
        else:
            n_norm = int(nreads_total)

        # Histogram of crossover midpoints
        cxs = sub["cx_mid"].to_numpy()
        mags, bin_edges = np.histogram(cxs, bins=bins)

        # Normalize per read (frequency per bin per read)
        mags = mags / float(n_norm)

        # Convert x-axis to Mb
        bin_edges_mb = bin_edges / 1e6
        mids_mb = (bin_edges_mb[:-1] + bin_edges_mb[1:]) / 2.0

        # Build toyplot canvas
        canvas = toyplot.Canvas(width=600, height=350)
        axes = canvas.cartesian(
            xlabel="Chromosome position (Mb)",
            ylabel="Crossover frequency (per read per bin)",
            margin=60,
        )

        # Styling (kept from your script)
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

        axes.title.text = f"{prefix} — {chrom} (N={n_norm} reads)"

        # Plot bars + line
        bars = axes.bars((mags, bin_edges_mb), opacity=0.6)
        line = axes.plot(mids_mb, mags, opacity=0.8)

        # Colors (your custom palette)
        bars._fill.color = "#d1cef6"
        line._stroke.color = "#7b6de2"

        if one_file_per_chrom:
            out = outdir / f"{prefix}.{chrom}.pdf"
            toyplot.pdf.render(canvas, str(out))
            logger.info(f"wrote {out}")
        else:
            pages.append(canvas)

    # Write multi-page PDF if requested
    if not one_file_per_chrom:
        if not pages:
            raise ValueError("No plots were generated (all chromosomes empty after filtering/masking).")
        toyplot.pdf.render(pages, str(pdf_out))
        logger.info(f"wrote {pdf_out}")


if __name__ == "__main__":
    # This file is meant to be called through the disto CLI.

