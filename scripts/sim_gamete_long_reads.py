"""
Simulate long reads from a pool of gametes made by recombining two haplotypes.

Usage
-----
python scripts/sim_gamete_pool_long_reads.py \
  --hap1 alt_refs/A_thaliana_hap1.fa \
  --hap2 alt_refs/A_thaliana_hap2.fa \
  --read-len 100000 \
  --n-reads 1000 \
  --mean-xovers 0.1 \
  --seed 42 \
  --out reads/athal_gametes_hifi.fq.gz

Notes
-----
- Recombination is modeled as a Poisson(k; mean = --mean-xovers) number of
  crossovers *per read*, with crossover breakpoints placed uniformly at random
  along the read and segments alternating between haplotypes.
- Starts are sampled uniformly across contigs shared by both haplotypes,
  weighted by contig length. The segment is taken from the same coordinates
  in hap1/hap2 (assumes coordinate-compatible haplotypes).
- FASTQ qualities are constant Q60 ('I').
"""

import argparse
import gzip
import os
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

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
    """
    Return contigs present in BOTH haplotypes with usable length.
    """
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
    """
    Choose a contig name proportional to its length.
    """
    total = sum(L for _, L in contigs)
    r = rng.randrange(total)
    acc = 0
    for name, L in contigs:
        acc += L
        if r < acc:
            return name, L
    return contigs[-1]

# -------------------- recombination model --------------------

def breakpoints_for_read(read_len: int, mean_xovers: float, rng: random.Random) -> List[int]:
    """
    Draw k ~ Poisson(mean_xovers), then place k breakpoints uniformly in (1, read_len-1).
    Return sorted list of breakpoint positions (0-based, exclusive end of a segment).
    """
    # Python's random has no Poisson; simple thinning using numpy would be nice,
    # but we'll keep stdlib only. Use a Knuth Poisson sampler.
    L = pow(2.718281828459045, -mean_xovers)
    k, p = 0, 1.0
    while True:
        p *= rng.random()
        if p <= L:
            break
        k += 1
    if k == 0:
        return []
    # unique breakpoints
    bps = sorted(set(rng.randrange(1, read_len) for _ in range(k)))
    # remove any exact duplicates (set handled), ensure not at end
    bps = [b for b in bps if b < read_len]
    return bps

def mosaic_read(h1_seq: str, h2_seq: str, start: int, read_len: int,
                mean_xovers: float, rng: random.Random) -> str:
    """
    Build a recombinant read by alternating haplotypes at random breakpoints.
    """
    bps = breakpoints_for_read(read_len, mean_xovers, rng)
    # segment boundaries: [0, bp1, bp2, ..., read_len]
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

# -------------------- main --------------------

def simulate_reads(h1_fa: str, h2_fa: str, read_len: int, n_reads: int,
                   mean_xovers: float, seed: int, out_path: str):
    rng = random.Random(seed)

    # load haplotypes
    h1 = read_fasta(h1_fa)
    h2 = read_fasta(h2_fa)
    contigs = contig_index(h1, h2)

    # filter contigs that can accommodate read_len
    contigs = [(c, L) for c, L in contigs if L >= read_len]
    if not contigs:
        raise RuntimeError(f"No contigs >= read_len={read_len} in both haplotypes.")

    # writer (gz or plain)
    is_gz = out_path.endswith(".gz")
    os.makedirs(Path(out_path).parent, exist_ok=True)
    if is_gz:
        out_fh = gzip.open(out_path, "wt")
    else:
        out_fh = open(out_path, "w")

    qline = "I" * read_len  # Q40 for HiFi-like

    try:
        for i in range(1, n_reads + 1):
            cname, L = weighted_choice(contigs, rng)
            # choose start so the whole read fits
            start = rng.randrange(0, L - read_len + 1)
            seq = mosaic_read(h1[cname], h2[cname], start, read_len, mean_xovers, rng)
            # FASTQ record
            out_fh.write(f"@simRead_{i}_{cname}:{start}-{start+read_len}\n")
            out_fh.write(seq + "\n+\n")
            out_fh.write(qline + "\n")
            if i % 1000 == 0:
                # small progress hint for big jobs
                pass
    finally:
        out_fh.close()

def positive_int(x: str) -> int:
    v = int(x)
    if v <= 0:
        raise argparse.ArgumentTypeError("Must be positive")
    return v

def main():
    ap = argparse.ArgumentParser(description="Simulate HiFi-like long reads from recombining haplotypes.")
    ap.add_argument("--hap1", required=True, help="FASTA for haplotype 1")
    ap.add_argument("--hap2", required=True, help="FASTA for haplotype 2")
    ap.add_argument("--read-len", type=positive_int, default=100_000, help="Length of each simulated read (bp)")
    ap.add_argument("--n-reads", type=positive_int, default=10_000, help="Number of reads to simulate")
    ap.add_argument("--mean-xovers", type=float, default=1.0,
                    help="Mean number of crossovers per read (Poisson mean). Set to 0 for no recombination.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--out", required=True, help="Output FASTQ(.gz)")
    args = ap.parse_args()

    simulate_reads(args.hap1, args.hap2, args.read_len, args.n_reads,
                   args.mean_xovers, args.seed, args.out)

if __name__ == "__main__":
    main()

