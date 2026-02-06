#!/usr/bin/env python
"""
Map long reads to a reference haplotype and call variants.

Pipeline:
1) minimap2 -> samtools sort -> BAM (+ index)
2) bcftools mpileup/call -> VCF.gz (+ index)
3) whatshap phase -> phased VCF.gz (+ index)

Robustness:
- WhatsHap thread flag differs across installs (--threads vs --jobs). We detect.
- WhatsHap may write UNCOMPRESSED VCF even if output filename ends with .gz.
  We therefore always write a tmp .vcf, then bgzip it via bcftools to .vcf.gz.
"""

from __future__ import annotations

import sys
from pathlib import Path
import subprocess as sp
from loguru import logger
import os

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_WHA = str(BIN / "whatshap")
BIN_MIN = str(BIN / "minimap2")


def map_reads_to_bam(reference: Path, gametes: Path, base: Path, threads: int) -> Path:
    """Map HiFi reads, sort/index BAM. Returns path to sorted BAM."""
    logger.info("aligning reads to reference")
    threads = max(1, int(threads))
    bam_file = base.with_suffix(".sorted.bam")

    # Ensure reference index exists for mpileup.
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        sp.run([BIN_SAM, "faidx", str(reference)], check=True)

    cmd_map = [
        BIN_MIN, "-ax", "map-hifi",
        "-R", f"@RG\\tID:{base.name}\\tSM:{base.name}",
        "-t", str(threads),
        str(reference),
        str(gametes),
    ]
    cmd_sort = [
        BIN_SAM, "sort",
        "-@", str(threads),
        "-O", "bam",
        "-o", str(bam_file),
        "-",  # read SAM from stdin
    ]

    # Stream minimap2 -> samtools sort
    p1 = sp.Popen(cmd_map, stdout=sp.PIPE, stderr=None)
    p2 = sp.Popen(cmd_sort, stdin=p1.stdout, stderr=None)
    assert p1.stdout is not None
    p1.stdout.close()

    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_map)
    if rc2 != 0:
        raise sp.CalledProcessError(rc2, cmd_sort)

    # BAM index
    sp.run([BIN_SAM, "index", "-@", str(threads), str(bam_file)], check=True)
    return bam_file


def call_variants_bcftools(
    reference: Path,
    bam_file: Path,
    base: Path,
    min_map_q: int,
    min_base_q: int,
    threads: int,
) -> Path:
    """Call variants with bcftools. Returns bgzipped VCF (.vcf.gz)."""
    logger.info("calling variants")
    threads = max(1, int(threads))

    final_vcf_gz = base.with_suffix(".vcf.gz")
    tmp_vcf_gz = Path(str(final_vcf_gz) + ".tmp")

    cmd_mpileup = [
        BIN_BCF, "mpileup",
        "--threads", str(threads),
        "-f", str(reference),
        "-q", str(min_map_q),
        "-Q", str(min_base_q),
        "-Ou",
        "-a", "AD,DP,SP",
        str(bam_file),
    ]
    cmd_call = [
        BIN_BCF, "call",
        "--threads", str(threads),
        "-m",
        "-v",
        "--ploidy", "2",
        "-Oz",
        "-o", str(tmp_vcf_gz),   # write TMP, not final
    ]

    p1 = sp.Popen(cmd_mpileup, stdout=sp.PIPE, stderr=None)
    p2 = sp.Popen(cmd_call, stdin=p1.stdout, stderr=None)
    assert p1.stdout is not None
    p1.stdout.close()

    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_mpileup)
    if rc2 != 0:
        raise sp.CalledProcessError(rc2, cmd_call)

    # Validate compressed stream is complete (works for bgzip too)
    sp.run(["gzip", "-t", str(tmp_vcf_gz)], check=True)

    # Atomic publish
    os.replace(tmp_vcf_gz, final_vcf_gz)

    # Index the published VCF
    sp.run([BIN_BCF, "index", "-f", str(final_vcf_gz)], check=True)
    return final_vcf_gz


def _whatshap_thread_flag() -> str | None:
    """Detect whether whatshap supports --threads or --jobs (or neither)."""
    try:
        help_txt = sp.check_output([BIN_WHA, "phase", "-h"], text=True)
    except Exception:
        return None

    if "--threads" in help_txt:
        return "--threads"
    if "--jobs" in help_txt:
        return "--jobs"
    return None


def phase_vcf(reference: Path, bam_path: Path, vcf_gz: Path, base: Path, threads: int) -> Path:
    """
    Phase VCF file with whatshap.

    Always outputs bgzipped .phased.vcf.gz, regardless of whether WhatsHap itself
    writes compressed or uncompressed output.
    """
    logger.info("phasing VCF")
    threads = max(1, int(threads))

    final_phased_gz = base.with_suffix(".phased.vcf.gz")

    # Always write an UNCOMPRESSED tmp VCF first
    tmp_vcf = Path(str(final_phased_gz) + ".tmp.vcf")
    # Then compress it ourselves to tmp.gz, validate, and publish atomically
    tmp_gz = Path(str(final_phased_gz) + ".tmp.gz")

    thread_flag = _whatshap_thread_flag()

    cmd = [BIN_WHA, "phase"]
    if thread_flag is not None:
        cmd.extend([thread_flag, str(threads)])
    cmd.extend([
        "--reference", str(reference),
        "-o", str(tmp_vcf),
        str(vcf_gz),
        str(bam_path),
    ])

    sp.run(cmd, check=True)

    # bgzip the VCF -> .vcf.gz using bcftools (already in your env)
    sp.run([BIN_BCF, "view", "-Oz", "-o", str(tmp_gz), str(tmp_vcf)], check=True)

    # Validate compressed stream is complete
    sp.run(["gzip", "-t", str(tmp_gz)], check=True)

    # Atomic publish
    os.replace(tmp_gz, final_phased_gz)

    # Index published VCF.gz
    sp.run([BIN_BCF, "index", "-f", str(final_phased_gz)], check=True)

    # Cleanup
    try:
        tmp_vcf.unlink()
    except FileNotFoundError:
        pass

    return final_phased_gz


def run_mapcall(
    reference: Path,
    gametes: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    min_map_q: int,
    min_base_q: int,
) -> Path:
    """Map reads, call variants, then phase. Returns phased VCF.gz."""
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    reference = reference.expanduser().absolute()
    gametes = gametes.expanduser().absolute()

    base = (outdir / prefix) if prefix else (outdir / gametes.stem)

    # 1) Map reads -> BAM
    bam_file = map_reads_to_bam(reference, gametes, base, threads)

    # 2) Call variants -> VCF.gz
    vcf_gz = call_variants_bcftools(reference, bam_file, base, min_map_q, min_base_q, threads)

    logger.info(f"BAM: {bam_file}")
    logger.info(f"VCF: {vcf_gz}")

    # 3) Phase VCF -> phased VCF.gz
    phased_vcf_gz = phase_vcf(reference, bam_file, vcf_gz, base, threads)
    logger.info(f"PHASED VCF: {phased_vcf_gz}")

    return phased_vcf_gz

