#!/usr/bin/env python

from __future__ import annotations
from typing import List, Dict, Tuple, Iterator
import gzip
from pathlib import Path
from loguru import logger
import numpy as np
from .rate_map import make_sine_rate_map, RateMapInterferenceSampler

RNG = np.random.Generator


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Return dict of {header: seq} where '>' is removed from header."""

    xopen = gzip.open if path.name.endswith(".gz") else open
    with xopen(path, 'rt') as indat:
        header = None
        seq = []
        for line in indat:
            if line.startswith(">"):
                if seq:
                    yield (header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        yield (header, "".join(seq))        


def parse_fasta(path: Path) -> Dict[str, str]:
    """..."""
    return {name: seq for (name, seq) in iter_fasta(path)}


class FastaGzipWriter:
    """Stream FASTA records to a .gz with fixed-width wrapping (fast, low overhead)."""

    def __init__(
        self,
        path: Path,
        width: int = 80,
        buf_lines: int = 16_384,
        compresslevel: int = 5,
        mtime: int | None = 0,
        append: bool = False,
    ) -> None:
        """
        Args:
            path: Output .fa.gz path.
            width: Line width (e.g., 80 for FASTA).
            buf_lines: Lines per batch; higher = fewer writes.
            compresslevel: gzip compression level (1–9).
            mtime: Reproducible timestamp for gzip header (0 for deterministic;
                   None to use current time). Ignored if unsupported by Python.
            append: If True, open in 'ab' and append as a new gzip member.
        """
        self.path = Path(path)
        self.width = int(width)
        self.buf_lines = int(buf_lines)
        self._nl = b"\n"

        mode = "ab" if append else "wb"
        # Try mtime; fall back for Python versions where it's unsupported.
        try:
            self._fh = gzip.open(
                self.path, mode, compresslevel=compresslevel, mtime=mtime
            )
        except TypeError:
            self._fh = gzip.open(self.path, mode, compresslevel=compresslevel)

    def close(self) -> None:
        """Close the gzip file."""
        self._fh.close()

    def write_chrom(self, name: str, seq: str) -> None:
        """Write a single chromosome record (header + wrapped sequence)."""
        # Header
        self._fh.write(b">")
        self._fh.write(name.encode("ascii"))
        self._fh.write(self._nl)

        # Sequence: batched fixed-width wrapping
        bseq = memoryview(seq.encode("ascii"))
        step = self.width * self.buf_lines
        for start in range(0, len(bseq), step):
            block = bseq[start : start + step]
            lines = [
                block[i : i + self.width] for i in range(0, len(block), self.width)
            ]
            self._fh.write(self._nl.join(lines) + self._nl)

    def __enter__(self) -> "FastaGzipWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def mutate_fasta(reference: Path, outfile: Path, het: float, chroms: List[str], rng: RNG):
    """Write two mutated copies of the reference genome"""
    nts = ("A", "C", "G", "T")
    snp_rate = het / 2.

    with FastaGzipWriter(outfile, width=80, buf_lines=16384, compresslevel=5) as w:
        for header, seq in iter_fasta(reference):

            # skip other chroms
            if chroms:
                cname = header.split()[0]
                if cname not in chroms:
                    continue

            # samples N mutated positions at snp_rate
            logger.info(f"writing mutated haplotype of chromosome {header.split()[0]}")
            mutpos = np.where(np.random.binomial(n=1, p=snp_rate, size=len(seq)))[0]
            seq = list(seq)
            for pos in mutpos:
                base = seq[pos].upper()
                if base in nts:
                    seq[pos] = rng.choice([x for x in nts if x != base])
            w.write_chrom(header, "".join(seq))
    return outfile


def mosaic_read(hap1: str, hap2: str, start: int, read_len: int, crossover: float, rng: RNG) -> str:
    """Return a read sampled from one haplotype or a recombined chrom."""
    if rng.binomial(n=1, p=0.5):
        a, b = hap1, hap2
    else:
        a, b = hap2, hap1
    if crossover > 0:
        chrom = a[:crossover] + b[crossover:]
    else:
        chrom = a
    return chrom[start: start + read_len]


def run_simulate(
    reference: Path,
    outdir: Path,
    prefix: str,
    read_length: int,
    nreads: int,
    heterozygosity: float,
    crossover_rate: float,
    interference: float,
    chromosomes: List[str] | None,
    random_seed: int | None,
):  
    """

    """
    # init RNG
    rng = np.random.default_rng(random_seed)

    # out paths
    outdir.mkdir(exist_ok=True, parents=True)
    hap1 = outdir / f"{prefix}.hap1.fa.gz"
    hap2 = outdir / f"{prefix}.hap2.fa.gz"
    gametes = outdir / f"{prefix}.gametes.fastq.gz"

    # write mutated haplotypes
    hap1 = mutate_fasta(reference, hap1, heterozygosity, chromosomes, rng)
    hap2 = mutate_fasta(reference, hap2, heterozygosity, chromosomes, rng)    

    # parse fastas to {cname: seq}
    h1 = parse_fasta(path=hap1)
    h2 = parse_fasta(path=hap2)

    # set quality scores on reads to high
    qline = "I" * read_length

    # build a crossover probability density sampler for each chrom
    samplers = {}
    for cname, seq in h1.items():
        # make a sine rate map
        rmap = make_sine_rate_map(len(seq), bins=500, mean_rate=crossover_rate)
        # accepts any rate map
        samplers[cname] = RateMapInterferenceSampler(rmap)

    # generate recomb reads and write.
    recom_read_counter = 0
    read_id = 0
    with gzip.open(gametes, "w", compresslevel=5) as fh:
        reads = []
        while read_id < nreads:
            for chrom in h1:
                # sample crossovers from this scaffold
                sampler = samplers[chrom]
                length = sampler.rate_map.length
                crossovers = sampler.sample_meiosis_gamma(nu=interference, rng=rng)
                
                # sample a read start position on this scaffold
                start = rng.integers(0, length - read_length + 1)

                # check if read will cover a crossover in the scaffold
                cx = 0
                for crossover in crossovers:
                    if start < crossover < start + read_length:
                        cx = crossover
                        recom_read_counter += 1
                        if not recom_read_counter % 100:
                            logger.debug(f"recombinant read {recom_read_counter} / {read_id} reads")

                # if so, store the read
                if cx:
                    seq = mosaic_read(h1[chrom], h2[chrom], start, read_length, cx, rng)
                    cname = chrom.split()[0]
                    read = f"@simRead_{cname}:{start}-{start+read_length}:cx={cx}:rid={read_id}\n{seq}\n+\n{qline}\n"
                    reads.append(read)
                    
                    # flush occasionally
                    if len(reads) >= 5_000:
                        fh.write("".join(reads).encode("utf-8"))
                        reads = []
                read_id += 1
        # write reads to file
        if reads:
            logger.debug(f"recombinant read {recom_read_counter} / {read_id} reads")            
            fh.write("".join(reads).encode("utf-8"))
    return gametes




if __name__ == "__main__":
    pass
