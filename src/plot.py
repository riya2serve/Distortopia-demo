#!/usr/bin/env python

"""
"""

from pathlib import Path
import numpy as np
import pandas as pd
import toyplot
import toyplot.pdf
from loguru import logger


def run_plot(
    tsv: Path,
    outdir: Path,
    prefix: str,
):
    out = outdir / f"{prefix}.pdf"
    # mean_rate = 1e-8
    bins = 100
    nreads_total = 1e7
    nreads_per_bin = nreads_total / bins
    read_length = 100_000
    data = pd.read_csv(tsv, sep=r"\s+")

    chroms = data.scaff.unique()
    for chrom in chroms:
        # select only reads on this chromosome
        sub = data.loc[data.scaff == chrom]

        # mask sites in the first and last read_length regions
        sub = sub.loc[sub.start > read_length]
        sub = sub.loc[sub.end < sub.end.max() - read_length]

        # select the true recomb positions
        cxs = sub.read.apply(lambda x: int(x.split("cx=")[1].split(":")[0]))
        mags, bins = np.histogram(cxs, bins=bins)

        # divide by total number of reads
        mags = mags / (nreads_total)
        # convert to Mb units
        bins = bins / 1e6

        canvas = toyplot.Canvas(width=600, height=350)
        axes = canvas.cartesian(
            xlabel="Chromosome position (Mb)", ylabel="Crossover frequency (per site per read)",
            margin=60,  
        )
        axes.x.ticks.show = True
        axes.x.ticks.far = 0
        axes.x.ticks.near = 5
        axes.x.ticks.labels.offset = 10
        axes.x.label.offset = 28
        axes.x.label.style['font-size'] = 14

        axes.y.ticks.show = True
        axes.y.ticks.far = 0
        axes.y.ticks.near = 5
        axes.y.ticks.labels.offset = 10
        axes.y.label.offset = 28
        axes.y.label.style['font-size'] = 14

        # 
        m0 = axes.bars((mags, bins), opacity=0.6)
        mids = (bins[:-1] + (bins[1] - bins[0]) / 2.)
        m1 = axes.plot(mids, mags, opacity=0.8)
        
        # --- added lines to recolor plots ---
        m0._fill.color = "#d1cef6"     # bars (fill color)
        m1.stroke.color = "#7b6de2"   # line (stroke color)
        # ---------------

        # axes.hlines([mean_rate])
        toyplot.pdf.render(canvas, str(out))
        logger.info(f"wrote {out}")


if __name__ == "__main__":
    main()
