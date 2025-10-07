#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import pandas as pd
import pysam
from cyvcf2 import VCF
import streamlit as st

def load_phased_biallelic_snps(vcf_path: Path):
    phased = {}
    vcf = VCF(str(vcf_path))
    for rec in vcf:
        if (not rec.is_snp) or (len(rec.ALT) != 1): continue
        gts = rec.genotypes
        if not gts: continue
        gt = gts[0]
        if len(gt) < 3 or not bool(gt[2]): continue
        a0, a1 = gt[0], gt[1]
        if a0 == a1 or a0 not in (0,1) or a1 not in (0,1): continue
        alleles = [rec.REF] + rec.ALT
        phased.setdefault(rec.CHROM, []).append((rec.POS, alleles[a0], alleles[a1]))
    for c in list(phased): phased[c].sort(key=lambda x:x[0])
    return phased

def read_crossovers(bam_path: Path, phased_dict, min_snps: int) -> pd.DataFrame:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    rows = []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
        chrom = bam.get_reference_name(read.reference_id)
        if chrom not in phased_dict: continue
        ref2base = {}
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if qpos is None or rpos is None: continue
            ref2base[rpos+1] = read.query_sequence[qpos]
        rs, re = read.reference_start+1, read.reference_end
        snp_pos, bits = [], []
        for pos, a0, a1 in phased_dict[chrom]:
            if pos < rs or pos > re: continue
            b = ref2base.get(pos)
            if b == a0: snp_pos.append(pos); bits.append(0)
            elif b == a1: snp_pos.append(pos); bits.append(1)
        if len(bits) < min_snps: continue
        changes = [i for i in range(1,len(bits)) if bits[i]!=bits[i-1]]
        phased_str = "".join(map(str,bits))
        if len(changes)==0: xl,xr="NA","NA"
        elif len(changes)==1: i=changes[0]; xl,xr=snp_pos[i-1],snp_pos[i]
        else: continue
        rows.append({"scaff":chrom,"start":rs,"end":re,"nsnps":len(bits),
                     "phased_snps":phased_str,"crossover_left":xl,"crossover_right":xr,"read":read.query_name})
    return pd.DataFrame(rows, columns=["scaff","start","end","nsnps","phased_snps","crossover_left","crossover_right","read"])

def run(bam: Path, vcf_path: Path, out_tsv: Path, min_snps: int = 5) -> Path:
    phased = load_phased_biallelic_snps(Path(vcf_path))
    df = read_crossovers(Path(bam), phased, int(min_snps))
    out_tsv = Path(out_tsv); out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    return out_tsv

def streamlit_panel(state):
    st.header("5) Find crossovers per read")
    col1,col2 = st.columns(2)
    with col1:
        bam = Path(st.text_input("Haplotagged BAM", value=str(state.alignments / "thaliana_gametes_hifi__to__A_thaliana.haplotagged.bam")))
        vcf = Path(st.text_input("Phased VCF",      value=str(state.data       / "A_thaliana_hap1_vs_hap2.phased.vcf")))
    with col2:
        out_tsv = Path(st.text_input("Output TSV",  value=str(state.data / "A_thaliana_gamete_crossovers.tsv")))
        min_snps = st.number_input("Min phased SNPs/read", 1, 100, int(getattr(state,"min_snps",5)))
    if st.button("Run XO detection"):
        try:
            out = run(bam, vcf, out_tsv, int(min_snps))
            st.success(f"Wrote {out}")
            st.dataframe(pd.read_csv(out, sep="\t").head(), use_container_width=True)
        except Exception as e:
            st.error(str(e))

def _cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("-o","--out", required=True)
    ap.add_argument("--min-snps", type=int, default=5)
    args = ap.parse_args()
    print(run(Path(args.bam), Path(args.vcf), Path(args.out), int(args.min_snps)))

if __name__ == "__main__":
    _cli()

