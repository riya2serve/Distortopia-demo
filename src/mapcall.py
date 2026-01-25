#!/usr/bin/env python

"""Map long reads to a reference haplotype and call variants.

TODO
----
    - add back in more filtering options on the variant calls.
    - ...

# Consider alternative implementation. But can it work for phasing?
freebayes \
  -f ref.fa \
  --pooled-continuous \
  --min-alternate-fraction 0.01 \
  --min-alternate-count 3 \
  --use-best-n-alleles 2 \
  pooled.bam > calls_freebayes.vcf

"""

import sys
from pathlib import Path
import subprocess as sp
from loguru import logger

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_WHA = str(BIN / "whatshap")
BIN_MIN = str(BIN / "minimap2")


# def map_reads_to_bam(reference: Path, gametes: Path, base: Path, threads: int) -> Path:
#     threads = max(1, threads // 2)
#     bam_file = base.with_suffix(".sorted.bam")
#     # Ensure reference index exists for mpileup.
#     fai = reference.with_suffix(reference.suffix + ".fai")
#     if not fai.exists():
#         sp.run(["samtools", "faidx", str(reference)], check=True)
#     # --- map & sort ---
#     cmd_map = [
#         BIN_MIN, "-ax", "map-hifi",
#         "-R", f"@RG\\tID:{base}\\tSM:{base}",  # store sample name.
#         "-t", str(threads),
#         str(reference),
#         str(gametes),
#     ]
#     cmd_sort = [
#         BIN_SAM, "sort",
#         "-@", str(threads),
#         "-O", "bam",
#         "-o", str(bam_file),
#     ]
#     p1 = sp.Popen(cmd_map, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
#     p2 = sp.Popen(cmd_sort, stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
#     if p1.stdout:
#         p1.stdout.close()
#     sort_stdout, sort_stderr = p2.communicate()
#     map_stderr = p1.stderr.read() if p1.stderr else b""
#     rc1 = p1.wait()
#     if rc1 != 0:
#         raise sp.CalledProcessError(rc1, cmd_map, None, map_stderr)
#     if p2.returncode != 0:
#         raise sp.CalledProcessError(p2.returncode, cmd_sort, sort_stdout, sort_stderr)
#     # BAM index (not strictly required for sequential mpileup, but recommended).
#     sp.run([BIN_SAM, "index", str(bam_file)], check=True)
#     return bam_file


def map_reads_to_bam(reference: Path, gametes: Path, base: Path, threads: int) -> Path:
    """Map HiFi reads and return sorted indexed BAM.
    """
    logger.info("aligning reads to reference")
    threads = max(1, threads // 2)
    bam_file = base.with_suffix(".sorted.bam")

    fai = reference.with_suffix(reference.suffix + ".fai")
    if not fai.exists():
        sp.run([BIN_SAM, "faidx", str(reference)], check=True)

    cmd_map = [
        BIN_MIN, "-ax", "map-hifi",
        "-R", f"@RG\\tID:{base.name}\\tSM:{base.name}",
        "-t", str(threads),
        # "--secondary", "no",  # drop secondaries right away
        str(reference),
        str(gametes),
    ]
    cmd_sort = [
        BIN_SAM, "sort",
        "-@", str(threads),
        "-O", "bam",
        "-o", str(bam_file),
        "-T", str(base.with_suffix(".sorttmp")),
    ]

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

    sp.run([BIN_SAM, "index", "-@", str(threads), str(bam_file)], check=True)
    return bam_file




def call_variants_bcftools(reference: Path, bam_file: Path, base: Path, min_map_q: int, min_base_q: int):
    """Return VCF file with variants called by bcftools."""
    logger.info("calling variants")
    vcf_gz = base.with_suffix(".vcf.gz")

    # --- call variants ---
    cmd_mpileup = [
        BIN_BCF, "mpileup",
        "-f", str(reference),
        "-q", str(min_map_q),
        "-Q", str(min_base_q),
        "-Ou",
        "-a", "AD,DP,SP",
        str(bam_file),
    ]
    cmd_call = [
        BIN_BCF, "call",
        "-m",  # multiallelic caller
        "-v",  # variants only
        "--ploidy", "2", # compare to freebayes
        "-Oz",
        "-o",
        str(vcf_gz),
    ]

    p1 = sp.Popen(cmd_mpileup, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
    p2 = sp.Popen(cmd_call, stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
    if p1.stdout:
        p1.stdout.close()
    call_stdout, call_stderr = p2.communicate()
    mpileup_stderr = p1.stderr.read() if p1.stderr else b""
    rc1 = p1.wait()
    if rc1 != 0:
        raise sp.CalledProcessError(rc1, cmd_mpileup, None, mpileup_stderr)
    if p2.returncode != 0:
        raise sp.CalledProcessError(p2.returncode, cmd_call, call_stdout, call_stderr)

    # Index the bgzipped VCF.
    sp.run([BIN_BCF, "index", "-f", str(vcf_gz)], check=True)
    return vcf_gz


def phase_vcf(reference: Path, bam_path: Path, vcf_gz: Path, base: Path):
    """Phase VCF file in whatshap"""
    logger.info("phasing VCF")
    phased_vcf_gz = base.with_suffix(".phased.vcf.gz")
    cmd1 = [
        BIN_WHA, "phase",
        "--reference", str(reference),
        "-o", str(phased_vcf_gz),
        str(vcf_gz),
        str(bam_path),
    ]
    p = sp.Popen(cmd1, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
    out, err = p.communicate()
    if p.returncode != 0:
        raise sp.CalledProcessError(p.returncode, cmd1, out, err)
    sp.run([BIN_BCF, "index", "-f", str(phased_vcf_gz)], check=True)
    return phased_vcf_gz


def run_mapcall(
    reference: Path,
    gametes: Path,
    outdir: Path,
    prefix: str,
    threads: int,
    min_map_q: int,
    min_base_q: int,
):
    """Map reads, call variants, then phase. """
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)
    reference = reference.expanduser().absolute()
    gametes = gametes.expanduser().absolute()

    if prefix:
        base = outdir / f"{prefix}"
    else:
        base = outdir / gametes.stem

    #1. Map reads --> BAM
    bam_file = map_reads_to_bam(reference, gametes, base, threads)

    #2. Call variants --> VCF.gz
    vcf_gz = call_variants_bcftools(reference, bam_file, base, min_map_q, min_base_q)
    bam_file = base.with_suffix(".sorted.bam")
    vcf_gz = base.with_suffix(".vcf.gz")
    logger.info(bam_file)
    logger.info(vcf_gz)

    #3. Phase VCF with WhatsHap
    phased_vcf_gz = phase_vcf(reference, bam_file, vcf_gz, base)
    return phased_vcf_gz

