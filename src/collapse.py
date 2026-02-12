#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path
from typing import Optional, List, Tuple, Dict

import numpy as np
import pandas as pd
from loguru import logger


REQ = {"scaff", "crossover_left", "crossover_right", "read"}


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("NA", np.nan), errors="coerce")


def _merge_intervals(intervals: List[Tuple[int, int]], merge_gap: int = 0) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + merge_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(a, b) for a, b in merged]


def run_collapse(
    tsv: Path,
    out_tsv: Path,
    merge_gap: int = 0,
    keep_na: bool = False,
) -> None:
    tsv = Path(tsv).expanduser().resolve()
    out_tsv = Path(out_tsv).expanduser().resolve()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"reading {tsv}")
    df = pd.read_csv(tsv, sep="\t")

    missing = REQ - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing required columns: {sorted(missing)}")

    df["scaff"] = df["scaff"].astype(str)
    df["crossover_left"] = _to_numeric(df["crossover_left"])
    df["crossover_right"] = _to_numeric(df["crossover_right"])

    if not keep_na:
        df = df.dropna(subset=["crossover_left", "crossover_right"]).copy()

    if df.empty:
        logger.warning("no crossover rows to collapse; writing empty output")
        pd.DataFrame(columns=["scaff","locus_left","locus_right","cx_mid","n_reads","window_size"]).to_csv(
            out_tsv, sep="\t", index=False
        )
        return

    df["crossover_left"] = df["crossover_left"].astype(int)
    df["crossover_right"] = df["crossover_right"].astype(int)

    out_rows: List[Dict[str, object]] = []

    for scaff, sub in df.groupby("scaff", sort=False):
        intervals = list(zip(sub["crossover_left"].tolist(), sub["crossover_right"].tolist()))
        merged = _merge_intervals(intervals, merge_gap=merge_gap)

        # count read support per merged locus by overlap
        for L, R in merged:
            mask = (sub["crossover_left"] <= R) & (sub["crossover_right"] >= L)
            n_reads = int(sub.loc[mask, "read"].nunique())

            out_rows.append({
                "scaff": scaff,
                "locus_left": int(L),
                "locus_right": int(R),
                "cx_mid": int((L + R) // 2),
                "n_reads": n_reads,
                "window_size": int(R - L + 1),
            })

    out = pd.DataFrame(out_rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"wrote {out_tsv} ({len(out)} loci)")


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(prog="disto collapse")
    p.add_argument("--tsv", required=True, type=Path, help="infer output TSV")
    p.add_argument("--out", required=True, type=Path, help="locus TSV")
    p.add_argument("--merge-gap", type=int, default=0)
    p.add_argument("--keep-na", action="store_true")
    args = p.parse_args(argv)

    run_collapse(args.tsv, args.out, merge_gap=args.merge_gap, keep_na=args.keep_na)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

