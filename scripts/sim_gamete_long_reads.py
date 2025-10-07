# scripts/sim_gamete_long_reads.py
from pathlib import Path
import random
from typing import Dict, List, Tuple
import streamlit as st
from Bio import SeqIO

def _read_fasta_dict(path: str) -> Dict[str, str]:
    d = {}
    for rec in SeqIO.parse(path, "fasta"):
        d[rec.id] = str(rec.seq).upper()
    if not d: raise RuntimeError(f"No sequences read from {path}")
    return d

def _contig_index(h1: Dict[str,str], h2: Dict[str,str]) -> List[Tuple[str,int]]:
    shared = sorted(set(h1) & set(h2))
    out = []
    for c in shared:
        L = min(len(h1[c]), len(h2[c]))
        if L > 0: out.append((c, L))
    if not out: raise RuntimeError("No shared contigs between haplotypes.")
    return out

def _weighted_choice(contigs, rng):
    tot = sum(L for _,L in contigs); r = rng.randrange(tot); acc=0
    for name,L in contigs:
        acc += L
        if r < acc: return name, L
    return contigs[-1]

def _poisson_knuth(lam: float, rng: random.Random) -> int:
    if lam<=0: return 0
    L = (2.718281828459045)**(-lam); k,p = 0,1.0
    while True:
        p *= rng.random()
        if p <= L: return k
        k += 1

def _breakpoints(read_len: int, mean_xovers: float, rng: random.Random) -> List[int]:
    k = _poisson_knuth(mean_xovers, rng)
    if k==0: return []
    bps = sorted(set(rng.randrange(1, read_len) for _ in range(k)))
    return [b for b in bps if 0<b<read_len]

def _mosaic_read(h1_seq: str, h2_seq: str, start: int, read_len: int, mean_xovers: float, rng: random.Random) -> str:
    bps = _breakpoints(read_len, mean_xovers, rng)
    cuts = [0] + bps + [read_len]
    use_h1_first = rng.choice([True, False])
    out = []
    for i in range(len(cuts)-1):
        seg = cuts[i+1]-cuts[i]; s = start + cuts[i]
        out.append((h1_seq if (use_h1_first ^ (i%2==1)) else h2_seq)[s:s+seg])
    return "".join(out)

def _genome_size_from_hap(fa: str) -> int:
    return sum(len(rec.seq) for rec in SeqIO.parse(fa, "fasta"))

def run(h1_fa: str, h2_fa: str, read_len: int, n_reads: int, mean_xovers: float, seed: int, out_path: str):
    rng = random.Random(seed)
    h1 = _read_fasta_dict(h1_fa); h2 = _read_fasta_dict(h2_fa)
    contigs = [(c,L) for c,L in _contig_index(h1,h2) if L>=read_len]
    if not contigs: raise RuntimeError(f"No contigs >= {read_len} in both haplotypes.")
    out_p = Path(out_path); out_p.parent.mkdir(parents=True, exist_ok=True)
    qline = "I"*read_len
    with open(out_p, "w") as fh:
        for i in range(1, n_reads+1):
            cname, L = _weighted_choice(contigs, rng)
            start = rng.randrange(0, L - read_len + 1)
            seq = _mosaic_read(h1[cname], h2[cname], start, read_len, mean_xovers, rng)
            fh.write(f"@simRead_{i}_{cname}:{start}-{start+read_len}\n{seq}\n+\n{qline}\n")
    return str(out_p)

def streamlit_panel(state):
    st.header("2) Simulate HiFi long reads (recombining haplotypes)")
    colA, colB = st.columns(2)
    with colA:
        hap1_path = st.text_input("Haplotype 1 FASTA", value=str(state.alt_refs / "A_thaliana_hap1.fa"))
        out_name  = st.text_input("Output filename", value="thaliana_gametes_hifi.fq", key="simreads_outname")
    with colB:
        hap2_path = st.text_input("Haplotype 2 FASTA", value=str(state.alt_refs / "A_thaliana_hap2.fa"))
        species_cov = st.number_input("Target coverage (×)", 1.0, 300.0, float(state.target_cov), 1.0, key="simreads_cov")

    read_len    = int(state.read_len)
    mean_xovers = float(state.mean_xovers)

    # auto n_reads
    try:
        G = _genome_size_from_hap(hap1_path)
        n_reads = max(1, (G*species_cov + read_len - 1)//read_len)
        st.caption(f"Genome ≈ {G:,} bp → {species_cov:.0f}× at {read_len:,} bp ≈ {n_reads:,} reads")
    except Exception:
        n_reads = st.number_input("Number of reads", 1, 5_000_000, 10_000, 1_000, key="simreads_nreads")

    seed = st.number_input("Seed", 0, 2**31-1, 42, key="simreads_seed")

    if st.button("Simulate reads", key="simreads_go"):
        try:
            out_path = Path(state.long_reads) / out_name
            out = run(hap1_path, hap2_path, read_len, int(n_reads), mean_xovers, int(seed), str(out_path))
            st.success(f"Done → {out}")

            # Serve BYTES (no temp files / no stale handles)
            data_bytes = Path(out).read_bytes()
            st.download_button(
                "Download FASTQ",
                data=data_bytes,
                file_name=Path(out).name,
                mime="text/plain",
                key="simreads_dl",
            )
        except Exception as e:
            st.error(f"Simulation failed: {e}")
