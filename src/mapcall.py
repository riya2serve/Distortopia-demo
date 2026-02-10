#!/usr/bin/env python3
"""
mapcall.py

Map long reads to a reference haplotype, call variants in parallel per contig,
and phase variants in parallel per contig.

Parallelization strategy:
  - Variant calling: one bcftools mpileup+call per contig (multiprocessing)
  - Phasing: one whatshap phase per contig (multiprocessing)
  - Concatenate per-contig VCFs into a final genome-wide VCF (bcftools concat)

Important note on CPU usage:
  - This script uses multiprocessing pools sized by --threads.
  - Each worker runs a single tool invocation; tools are invoked as single-threaded
    per worker (except mapping/sorting, which uses --threads).
  - That means "true parallelization" comes from multiple contigs running at once.

Requirements in PATH (or in the active conda env):
  - minimap2, samtools, bcftools, whatshap
"""

from __future__ import annotations

import argparse
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

def _run(cmd: list[str], **kwargs) -> None:
    """Run a command with logging."""
    logger.debug("CMD: {}", " ".join(cmd))
    sp.run(cmd, check=True, **kwargs)


def contigs_from_fai(reference: Path) -> list[str]:
    """Get contigs from reference .fai (create if missing)."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        _run([BIN_SAM, "faidx", str(reference)])

    contigs: list[str] = []
    with open(fai) as f:
        for line in f:
            cid = line.split("\t", 1)[0].strip()
            if cid:
                contigs.append(cid)

    return contigs


def _whatshap_thread_flag() -> str | None:
    """Detect which threading flag whatshap supports (varies by version/build)."""
    try:
        txt = sp.check_output([BIN_WHA, "phase", "-h"], text=True)
    except Exception:
        return None
    if "--threads" in txt:
        return "--threads"
    if "--jobs" in txt:
        return "--jobs"
    return None


# -------------------------
# Mapping
# -------------------------

def map_reads_to_bam(reference: Path, reads: Path, base: Path, threads: int) -> Path:
    """Map reads, sort, and index BAM."""
    logger.info("Aligning reads to reference and sorting BAM")
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
    assert p1.stdout is not None
    p1.stdout.close()

    rc2 = p2.wait()
    rc1 = p1.wait()

    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_map)
    if rc2 != 0:
        raise sp.CalledProcessError(rc2, cmd_sort)

    _run([BIN_SAM, "index", "-@", str(threads), str(bam_file)])
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
    out_tmp = Path(str(out_vcf) + ".tmp.gz")

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
    assert p1.stdout is not None
    p1.stdout.close()

    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_mpileup)
    if rc2 != 0:
        raise sp.CalledProcessError(rc2, cmd_call)

    _run(["gzip", "-t", str(out_tmp)])
    os.replace(out_tmp, out_vcf)
    _run([BIN_BCF, "index", "-f", str(out_vcf)])

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
    workers: int,
    max_depth: int | None,
) -> tuple[Path, list[str], list[str]]:
    """
    Call variants per contig in parallel, then concatenate.

    Returns:
      final_vcf, contigs, per_contig_vcfs
    """
    contigs = contigs_from_fai(reference)
    if not contigs:
        raise RuntimeError("No contigs found in reference")

    workers = max(1, int(workers))
    logger.info(f"Calling variants on {len(contigs)} contigs with {workers} workers")

    work = [(c, str(reference), str(bam_file), min_map_q, min_base_q, max_depth) for c in contigs]

    with Pool(processes=workers) as pool:
        per_vcfs = pool.map(_mpileup_call_worker, work)

    final_vcf = base.with_suffix(".vcf.gz")
    tmp_vcf = Path(str(final_vcf) + ".tmp.gz")

    # IMPORTANT: concat in contig order (not completion order)
    part_files_ordered = [f"{bam_file}.{c}.vcf.gz" for c in contigs]

    _run([BIN_BCF, "concat", "-Oz", "-o", str(tmp_vcf), *part_files_ordered])
    _run(["gzip", "-t", str(tmp_vcf)])
    os.replace(tmp_vcf, final_vcf)
    _run([BIN_BCF, "index", "-f", str(final_vcf)])

    return final_vcf, contigs, per_vcfs


# -------------------------
# Phasing (worker)
# -------------------------

def _phase_worker(args: tuple[str, str, str, str, int]) -> str:
    """
    Phase one contig VCF with whatshap.
    One process == one contig.
    """
    contig, reference, bam, contig_vcf_gz, whatshap_threads = args

    reference = str(reference)
    bam = str(bam)
    contig_vcf_gz = str(contig_vcf_gz)

    out_gz = Path(f"{contig_vcf_gz}.phased.vcf.gz")
    tmp_vcf = Path(str(out_gz) + f".tmp.{os.getpid()}.vcf")
    tmp_gz = Path(str(out_gz) + f".tmp.{os.getpid()}.gz")

    cmd = [BIN_WHA, "phase"]
    tf = _whatshap_thread_flag()
    if tf and whatshap_threads > 1:
        cmd.extend([tf, str(whatshap_threads)])

    cmd.extend([
        "--reference", reference,
        "-o", str(tmp_vcf),
        contig_vcf_gz,
        bam,
    ])

    _run(cmd)

    # compress to bgzip with bcftools view (more reliable than gzip on VCF text)
    _run([BIN_BCF, "view", "-Oz", "-o", str(tmp_gz), str(tmp_vcf)])
    _run(["gzip", "-t", str(tmp_gz)])
    os.replace(tmp_gz, out_gz)
    _run([BIN_BCF, "index", "-f", str(out_gz)])

    try:
        tmp_vcf.unlink()
    except FileNotFoundError:
        pass

    return str(out_gz)


# -------------------------
# Parallel phasing
# -------------------------

def phase_variants_parallel(
    reference: Path,
    bam_file: Path,
    contigs: list[str],
    per_contig_vcfs: list[str],
    base: Path,
    workers: int,
    whatshap_threads_per_worker: int = 1,
) -> Path:
    """
    Phase per-contig VCFs in parallel, then concatenate into a final phased VCF.

    - workers: number of concurrent contigs to phase
    - whatshap_threads_per_worker: threads passed to each whatshap invocation.
      Keep this at 1 if you want "true parallelization" across contigs using a fixed CPU budget.
      (If you set >1, you may oversubscribe CPUs unless you reduce workers accordingly.)
    """
    workers = max(1, int(workers))
    whatshap_threads_per_worker = max(1, int(whatshap_threads_per_worker))

    logger.info(
        f"Phasing {len(contigs)} contigs with {workers} workers "
        f"(whatshap threads per worker: {whatshap_threads_per_worker})"
    )

    # Ensure ordered mapping: contigs -> expected per-contig vcf name
    # We ignore per_contig_vcfs list order and build canonical list in contig order.
    ordered_inputs = [f"{bam_file}.{c}.vcf.gz" for c in contigs]

    work = [
        (c, str(reference), str(bam_file), ordered_inputs[i], whatshap_threads_per_worker)
        for i, c in enumerate(contigs)
    ]

    with Pool(processes=workers) as pool:
        pool.map(_phase_worker, work)

    final_gz = base.with_suffix(".phased.vcf.gz")
    tmp_gz = Path(str(final_gz) + ".tmp.gz")

    phased_parts = [f"{inp}.phased.vcf.gz" for inp in ordered_inputs]

    _run([BIN_BCF, "concat", "-Oz", "-o", str(tmp_gz), *phased_parts])
    _run(["gzip", "-t", str(tmp_gz)])
    os.replace(tmp_gz, final_gz)
    _run([BIN_BCF, "index", "-f", str(final_gz)])

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
    whatshap_threads_per_worker: int = 1,
) -> Path:
    """Map → parallel call → parallel phase."""
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    reference = reference.expanduser().absolute()
    reads = reads.expanduser().absolute()

    base = outdir / (prefix if prefix else reads.stem)

    bam = map_reads_to_bam(reference, reads, base, threads)

    vcf, contigs, per_vcfs = call_variants_parallel(
        reference=reference,
        bam_file=bam,
        base=base,
        min_map_q=min_map_q,
        min_base_q=min_base_q,
        workers=threads,
        max_depth=max_depth,
    )

    phased = phase_variants_parallel(
        reference=reference,
        bam_file=bam,
        contigs=contigs,
        per_contig_vcfs=per_vcfs,
        base=base,
        workers=threads,
        whatshap_threads_per_worker=whatshap_threads_per_worker,
    )

    logger.info(f"BAM: {bam}")
    logger.info(f"VCF: {vcf}")
    logger.info(f"PHASED VCF: {phased}")

    return phased

