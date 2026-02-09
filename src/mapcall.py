#!/usr/bin/env python
"""
Map long reads to a reference haplotype, call variants in parallel per contig,
and phase variants.

Parallelization strategy:
  - One bcftools mpileup+call per contig (multiprocessing)
  - One CPU per contig worker
  - Final bcftools concat + index

Notes:
  - bcftools mpileup is single-threaded per region; parallelism is per contig
  - Depth cap is configurable  by user
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
import subprocess as sp
from multiprocessing import Pool
from loguru import logger


# -------------------------
# Tool paths
# -------------------------

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_WHA = str(BIN / "whatshap")
BIN_MIN = str(BIN / "minimap2")


# -------------------------
# Helpers
# -------------------------

def contigs_from_fai(reference: Path) -> list[str]:
    """Get contigs from reference .fai (create if missing)."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        sp.run([BIN_SAM, "faidx", str(reference)], check=True)

    contigs: list[str] = []
    with open(fai) as f:
        for line in f:
            cid = line.split("\t", 1)[0]
            if cid:
                contigs.append(cid)
    return contigs


# -------------------------
# Mapping
# -------------------------

def map_reads_to_bam(reference: Path, reads: Path, base: Path, threads: int) -> Path:
    """Map reads, sort, and index BAM."""
    logger.info("aligning reads to reference")
    threads = max(1, int(threads))

    bam_file = base.with_suffix(".sorted.bam")

    cmd_map = [
        BIN_MIN, "-ax", "map-hifi",
        "-R", f"@RG\\tID:{base.name}\\tSM:{base.name}",
        "-t", str(threads),
        str(reference),
        str(reads),
    ]
    cmd_sort = [
        BIN_SAM, "sort",
        "-@", str(threads),
        "-O", "bam",
        "-o", str(bam_file),
        "-",
    ]

    p1 = sp.Popen(cmd_map, stdout=sp.PIPE)
    p2 = sp.Popen(cmd_sort, stdin=p1.stdout)
    p1.stdout.close()
    rc2 = p2.wait()
    rc1 = p1.wait()

    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_map)
    if rc2 != 0:
        raise sp.CalledProcessError(rc2, cmd_sort)

    sp.run([BIN_SAM, "index", "-@", str(threads), str(bam_file)], check=True)
    return bam_file


# -------------------------
# Variant calling (worker)
# -------------------------

def _mpileup_call_worker(args: tuple[str, str, str, int, int, int | None]) -> str:
    """
    Call variants for a single contig.
    One process == one contig.
    """
    contig, reference, bam, min_map_q, min_base_q, max_depth = args

    out_vcf = Path(f"{bam}.{contig}.vcf.gz")
    out_tmp = Path(str(out_vcf) + ".tmp")

    cmd_mpileup = [
        BIN_BCF, "mpileup",
        "-r", contig,
        "-Ou",
        "-f", reference,
        "-q", str(min_map_q),
        "-Q", str(min_base_q),
    ]

    if max_depth is not None:
        cmd_mpileup.extend(["-d", str(max_depth)])

    cmd_mpileup.append(bam)

    cmd_call = [
        BIN_BCF, "call",
        "-mv",
        "-Oz",
        "-o", str(out_tmp),
    ]

    p1 = sp.Popen(cmd_mpileup, stdout=sp.PIPE)
    p2 = sp.Popen(cmd_call, stdin=p1.stdout)
    p1.stdout.close()
    p2.wait()
    p1.wait()

    sp.run(["gzip", "-t", str(out_tmp)], check=True)
    os.replace(out_tmp, out_vcf)
    sp.run([BIN_BCF, "index", "-f", str(out_vcf)], check=True)

    return str(out_vcf)


# -------------------------
# Parallel variant calling
# -------------------------

def call_variants_parallel(
    reference: Path,
    bam_file: Path,
    base: Path,
    min_map_q: int,
    min_base_q: int,
    threads: int,
    max_depth: int | None,
) -> Path:
    """Call variants per contig in parallel, then concat."""
    contigs = contigs_from_fai(reference)
    if not contigs:
        raise RuntimeError("No contigs found in reference")

    logger.info(f"Calling variants on {len(contigs)} contigs with {threads} workers")

    work = [
        (c, str(reference), str(bam_file), min_map_q, min_base_q, max_depth)
        for c in contigs
    ]

    with Pool(processes=threads) as pool:
        pool.map(_mpileup_call_worker, work)

    final_vcf = base.with_suffix(".vcf.gz")
    tmp_vcf = Path(str(final_vcf) + ".tmp")

    part_files = [f"{bam_file}.{c}.vcf.gz" for c in contigs]

    sp.run(
        [BIN_BCF, "concat", "-Oz", "-o", str(tmp_vcf), *part_files],
        check=True,
    )

    sp.run(["gzip", "-t", str(tmp_vcf)], check=True)
    os.replace(tmp_vcf, final_vcf)
    sp.run([BIN_BCF, "index", "-f", str(final_vcf)], check=True)

    return final_vcf


# -------------------------
# Phasing
# -------------------------

def _whatshap_thread_flag() -> str | None:
    try:
        txt = sp.check_output([BIN_WHA, "phase", "-h"], text=True)
    except Exception:
        return None
    if "--threads" in txt:
        return "--threads"
    if "--jobs" in txt:
        return "--jobs"
    return None


def phase_vcf(reference: Path, bam: Path, vcf_gz: Path, base: Path, threads: int) -> Path:
    logger.info("phasing VCF")
    threads = max(1, int(threads))

    final_gz = base.with_suffix(".phased.vcf.gz")
    tmp_vcf = Path(str(final_gz) + ".tmp.vcf")
    tmp_gz = Path(str(final_gz) + ".tmp.gz")

    cmd = [BIN_WHA, "phase"]
    tf = _whatshap_thread_flag()
    if tf:
        cmd.extend([tf, str(threads)])

    cmd.extend([
        "--reference", str(reference),
        "-o", str(tmp_vcf),
        str(vcf_gz),
        str(bam),
    ])

    sp.run(cmd, check=True)

    sp.run([BIN_BCF, "view", "-Oz", "-o", str(tmp_gz), str(tmp_vcf)], check=True)
    sp.run(["gzip", "-t", str(tmp_gz)], check=True)
    os.replace(tmp_gz, final_gz)
    sp.run([BIN_BCF, "index", "-f", str(final_gz)], check=True)

    try:
        tmp_vcf.unlink()
    except FileNotFoundError:
        pass

    return final_gz


# -------------------------
# Main pipeline entry
# -------------------------

def run_mapcall(
    reference: Path,
    reads: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    min_map_q: int,
    min_base_q: int,
    max_depth: int | None = 200,
) -> Path:
    """Map → parallel call → phase."""
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    reference = reference.expanduser().absolute()
    reads = reads.expanduser().absolute()

    base = outdir / (prefix if prefix else reads.stem)

    bam = map_reads_to_bam(reference, reads, base, threads)

    vcf = call_variants_parallel(
        reference,
        bam,
        base,
        min_map_q,
        min_base_q,
        threads,
        max_depth,
    )

    phased = phase_vcf(reference, bam, vcf, base, threads)

    logger.info(f"BAM: {bam}")
    logger.info(f"VCF: {vcf}")
    logger.info(f"PHASED VCF: {phased}")

    return phased

