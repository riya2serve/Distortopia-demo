# scripts/haplotag_helper.py
from pathlib import Path
import streamlit as st, subprocess, shutil
from Bio import SeqIO

def have(t): return shutil.which(t) is not None
def sh(cmd):
    st.code("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=True, capture_output=True, text=True)

def make_simple_phased_vcf(h1_fa: Path, h2_fa: Path, out_vcf: Path):
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with open(out_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        H1 = {r.id:str(r.seq).upper() for r in SeqIO.parse(str(h1_fa), "fasta")}
        H2 = {r.id:str(r.seq).upper() for r in SeqIO.parse(str(h2_fa), "fasta")}
        for chrom in sorted(set(H1)&set(H2)):
            s1, s2 = H1[chrom], H2[chrom]
            for i,(a,b) in enumerate(zip(s1, s2), start=1):
                if a!=b and a in "ACGT" and b in "ACGT":
                    f.write(f"{chrom}\t{i}\t.\t{a}\t{b}\t.\tPASS\t.\tGT\t0|1\n")

def streamlit_panel(state):
    st.header("4) Haplotag reads with whatshap")

    col1, col2 = st.columns(2)
    with col1:
        ref_fa = Path(st.text_input(
            "Reference FASTA",
            value=str(state.raw_data / "A_thaliana.fna"),
            key="haplotag_ref_fa",
        ))
        bam = Path(st.text_input(
            "Aligned BAM",
            value=str(state.alignments / "Athaliana_gametes_hifi__to__A_thaliana.bam"),
            key="haplotag_bam",
        ))
    with col2:
        hap1 = Path(st.text_input(
            "Haplotype 1 FASTA",
            value=str(state.alt_refs / "A_thaliana_hap1.fa"),
            key="haplotag_hap1",
        ))
        hap2 = Path(st.text_input(
            "Haplotype 2 FASTA",
            value=str(state.alt_refs / "A_thaliana_hap2.fa"),
            key="haplotag_hap2",
        ))

    out_bam = Path(st.text_input(
        "Output haplotagged BAM",
        value=str(state.alignments / (bam.stem + ".haplotagged.bam")),
        key="haplotag_out_bam",
    ))
    phased_vcf = Path(st.text_input(
        "Phased VCF (0|1 SNPs)",
        value=str(state.data / "A_thaliana_hap1_vs_hap2.phased.vcf"),
        key="haplotag_phased_vcf",
    ))
    make_vcf = st.checkbox(
        "If VCF missing, build from hap FASTAs",
        value=True,
        key="haplotag_make_vcf",
    )

    if st.button("Run haplotag", key="haplotag_run_button"):
        for t in ["whatshap", "samtools"]:
            if not have(t):
                st.error(f"`{t}` not found on PATH"); return
        try:
            if make_vcf and not phased_vcf.exists():
                make_simple_phased_vcf(hap1, hap2, phased_vcf)
                st.success(f"Made phased VCF: {phased_vcf.name}")

            # ensure a SAMPLE is present in the BAM
            sm_bam = out_bam.with_suffix(".SM.bam")
            sh(["samtools","addreplacerg","-r","ID:sim\tSM:SAMPLE\tPL:HiFi","-o",str(sm_bam),str(bam)])
            sh(["samtools","index",str(sm_bam)])

            # haplotag
            sh(["whatshap","haplotag","--reference",str(ref_fa),
                "-o",str(out_bam), str(phased_vcf), str(sm_bam)])
            sh(["samtools","index",str(out_bam)])
            st.success(f"Haplotag complete → {out_bam.name}")
        except subprocess.CalledProcessError as e:
            st.error((e.stderr or e.stdout or str(e)).strip())

