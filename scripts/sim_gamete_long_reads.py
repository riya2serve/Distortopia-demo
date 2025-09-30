"""
Simulate PacBio HiFi long reads from a pool of gametes made by recombining two haplotypes.

CLI example
-----------
python scripts/sim_gamete_pool_long_reads.py \
  --hap1 alt_refs/A_thaliana_hap1.fa \
  --hap2 alt_refs/A_thaliana_hap2.fa \
  --read-len 100000 \
  --n-reads 1000 \
  --mean-xovers 1.0 \
  --seed 42 \
  --out reads/athal_gametes_hifi.fq.gz
"""

import argparse
import gzip
import os
import random
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import streamlit as st
from Bio import SeqIO   # biopython

# -------------------- FASTA helpers --------------------

def read_fasta(path: str) -> Dict[str, str]:
    seqs = {}
    with open(path, "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            seqs[rec.id] = str(rec.seq).upper()
    if not seqs:
        raise RuntimeError(f"No sequences read from {path}")
    return seqs

def contig_index(h1: Dict[str, str], h2: Dict[str, str]) -> List[Tuple[str, int]]:
    shared = sorted(set(h1.keys()) & set(h2.keys()))
    out = []
    for c in shared:
        L = min(len(h1[c]), len(h2[c]))
        if L > 0:
            out.append((c, L))
    if not out:
        raise RuntimeError("No shared contigs between haplotypes.")
    return out

def weighted_choice(contigs: List[Tuple[str, int]], rng: random.Random) -> Tuple[str, int]:
    total = sum(L for _, L in contigs)
    r = rng.randrange(total)
    acc = 0
    for name, L in contigs:
        acc += L
        if r < acc:
            return name, L
    return contigs[-1]

# -------------------- recombination model --------------------

def _poisson_knuth(mean_xovers: float, rng: random.Random) -> int:
    if mean_xovers <= 0:
        return 0
    L = pow(2.718281828459045, -mean_xovers)
    k, p = 0, 1.0
    while True:
        p *= rng.random()
        if p <= L:
            return k
        k += 1

def breakpoints_for_read(read_len: int, mean_xovers: float, rng: random.Random) -> List[int]:
    k = _poisson_knuth(mean_xovers, rng)
    if k == 0:
        return []
    bps = sorted(set(rng.randrange(1, read_len) for _ in range(k)))
    return [b for b in bps if 0 < b < read_len]

def mosaic_read(h1_seq: str, h2_seq: str, start: int, read_len: int,
                mean_xovers: float, rng: random.Random) -> str:
    bps = breakpoints_for_read(read_len, mean_xovers, rng)
    cuts = [0] + bps + [read_len]
    use_h1_first = rng.choice([True, False])
    out_chunks = []
    for i in range(len(cuts) - 1):
        seg_len = cuts[i+1] - cuts[i]
        s = start + cuts[i]
        if use_h1_first ^ (i % 2 == 1):
            out_chunks.append(h1_seq[s:s+seg_len])
        else:
            out_chunks.append(h2_seq[s:s+seg_len])
    return "".join(out_chunks)

# -------------------- core simulator --------------------

def simulate_reads(h1_fa: str, h2_fa: str, read_len: int, n_reads: int,
                   mean_xovers: float, rng: random.Random, out_path: str):
    h1 = read_fasta(h1_fa)
    h2 = read_fasta(h2_fa)
    contigs = [(c, L) for c, L in contig_index(h1, h2) if L >= read_len]
    if not contigs:
        raise RuntimeError(f"No contigs >= read_len={read_len} in both haplotypes.")

    os.makedirs(Path(out_path).parent, exist_ok=True)
    out_fh = gzip.open(out_path, "wt") if out_path.endswith(".gz") else open(out_path, "w")
    qline = "I" * read_len  # Q40 (HiFi-like)

    try:
        prog = st.progress(0, text="Writing reads...")
        for i in range(1, n_reads + 1):
            cname, L = weighted_choice(contigs, rng)
            start = rng.randrange(0, L - read_len + 1)
            seq = mosaic_read(h1[cname], h2[cname], start, read_len, mean_xovers, rng)
            out_fh.write(f"@simRead_{i}_{cname}:{start}-{start+read_len}\n{seq}\n+\n{qline}\n")
            if i % max(1, n_reads // 100) == 0:
                prog.progress(min(i / n_reads, 1.0))
        prog.progress(1.0)
    finally:
        out_fh.close()

# -------------------- Streamlit UI --------------------

st.title("Distortopia: Simulate HiFi long reads from recombining haplotypes")

st.markdown(
    "Upload **two haplotype FASTA files** (or point to local paths), "
    "choose a read length and number of reads, and generate a FASTQ.gz of simulated HiFi reads."
)

colA, colB = st.columns(2)
with colA:
    up1 = st.file_uploader("Haplotype 1 FASTA (.fa/.fasta/.fna)", type=["fa", "fasta", "fna"])
    hap1_path = st.text_input("OR path to haplotype 1 FASTA", value="alt_refs/A_thaliana_hap1.fa")
with colB:
    up2 = st.file_uploader("Haplotype 2 FASTA (.fa/.fasta/.fna)", type=["fa", "fasta", "fna"])
    hap2_path = st.text_input("OR path to haplotype 2 FASTA", value="alt_refs/A_thaliana_hap2.fa")

read_len    = st.number_input("Read length (bp)", min_value=1_000, max_value=2_000_000, value=100_000, step=1_000)
n_reads     = st.number_input("Number of reads", min_value=1, max_value=5_000_000, value=10_000, step=1_000)
mean_xovers = st.number_input("Mean crossovers per read", min_value=0.0, max_value=10.0, value=1.0, step=0.1)

repro = st.checkbox("Set seed for reproducibility", value=False)
seed  = st.number_input("Seed (integer)", min_value=0, max_value=2**31-1, value=42) if repro else None

out_name = st.text_input("Output filename", value="gametes_hifi.fq.gz")

if st.button("Simulate", type="primary"):
    with tempfile.TemporaryDirectory() as tdir:
        tdir = Path(tdir)

        # Materialize uploads or use provided paths
        if up1 is not None:
            p1 = tdir / "hap1.fa"
            p1.write_bytes(up1.getbuffer())
            h1 = str(p1)
        else:
            h1 = hap1_path

        if up2 is not None:
            p2 = tdir / "hap2.fa"
            p2.write_bytes(up2.getbuffer())
            h2 = str(p2)
        else:
            h2 = hap2_path

        out_path = str(tdir / out_name)

        # RNG: reproducible if seed set; otherwise cryptographic randomness
        rng = random.Random(seed if seed is not None else int.from_bytes(os.urandom(8), "big"))

        try:
            with st.spinner("Simulating reads..."):
                simulate_reads(h1, h2, int(read_len), int(n_reads), float(mean_xovers), rng, out_path)
            st.success("Done! Download below.")
            with open(out_path, "rb") as fh:
                st.download_button("Download FASTQ.gz", fh, file_name=out_name, mime="application/gzip")
        except Exception as e:
            st.error(f"Simulation failed: {e}")
