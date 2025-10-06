#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import pysam
from cyvcf2 import VCF
import pandas as pd


def load_phased_biallelic_snps(vcf_path: Path):
    """
    Return dict: chrom -> [(pos, hap0_base, hap1_base)], sorted by pos.
    Robust to cyvcf2 versions (uses len(rec.ALT) instead of num_alt/num_alts).
    """
    phased = {}
    vcf = VCF(str(vcf_path))
    for rec in vcf:
        # keep only biallelic SNPs
        if (not rec.is_snp) or (len(rec.ALT) != 1):
            continue

        gts = rec.genotypes
        if not gts:
            continue
        gt = gts[0]  # [allele0, allele1, phased_flag] for first/only sample

        # some cyvcf2 builds can return shorter lists; guard carefully
        if len(gt) < 3:
            continue
        phased_flag = bool(gt[2])
        if not phased_flag:
            continue

        a0, a1 = gt[0], gt[1]
        # heterozygous 0|1 or 1|0 only
        if a0 == a1 or (a0 not in (0, 1)) or (a1 not in (0, 1)):
            continue

        alleles = [rec.REF] + rec.ALT  # e.g., ['A', 'G']
        hap0, hap1 = alleles[a0], alleles[a1]
        phased.setdefault(rec.CHROM, []).append((rec.POS, hap0, hap1))

    # sort by position
    for chrom in list(phased):
        phased[chrom].sort(key=lambda x: x[0])
    return phased


def read_crossovers(bam_path: Path, phased_dict, min_snps: int) -> pd.DataFrame:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    rows = []

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        chrom = bam.get_reference_name(read.reference_id)
        if chrom not in phased_dict:
            continue

        # map reference positions -> read bases
        ref2base = {}
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if qpos is None or rpos is None:
                continue
            ref2base[rpos + 1] = read.query_sequence[qpos]

        rs, re = read.reference_start + 1, read.reference_end
        snp_pos, bits = [], []

        for pos, a0, a1 in phased_dict[chrom]:
            if pos < rs or pos > re:
                continue
            b = ref2base.get(pos)
            if b == a0:
                snp_pos.append(pos)
                bits.append(0)
            elif b == a1:
                snp_pos.append(pos)
                bits.append(1)

        if len(bits) < min_snps:
            continue

        # find hap switches
        changes = [i for i in range(1, len(bits)) if bits[i] != bits[i - 1]]
        phased_str = "".join(map(str, bits))

        if len(changes) == 0:
            xl, xr = "NA", "NA"
        elif len(changes) == 1:
            i = changes[0]
            xl, xr = snp_pos[i - 1], snp_pos[i]
        else:
            # >1 switch: discard per your spec
            continue

        rows.append({
            "scaff": chrom,
            "start": rs,
            "end": re,
            "nsnps": len(bits),
            "phased_snps": phased_str,
            "crossover_left": xl,
            "crossover_right": xr,
            "read": read.query_name,
        })

    return pd.DataFrame(
        rows,
        columns=["scaff", "start", "end", "nsnps", "phased_snps",
                 "crossover_left", "crossover_right", "read"],
    )


def main():
    ap = argparse.ArgumentParser(
        description="Detect single crossovers per read from haplotagged BAM + phased VCF"
    )
    ap.add_argument("--bam", required=True, help="haplotagged BAM (or plain BAM)")
    ap.add_argument("--vcf", required=True, help="phased VCF.gz (0|1 sites)")
    ap.add_argument("-o", "--out", required=True, help="output TSV path")
    ap.add_argument("--min-snps", type=int, default=5,
                    help="min phased SNPs on a read (default 5)")
    args = ap.parse_args()

    phased = load_phased_biallelic_snps(Path(args.vcf))
    df = read_crossovers(Path(args.bam), phased, args.min_snps)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}  (n={len(df)})")


if __name__ == "__main__":
    main()

