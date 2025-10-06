python - <<'PY'
from Bio import SeqIO
from pathlib import Path

hap1 = Path("alt_refs/A_thaliana_hap1.fa")
hap2 = Path("alt_refs/A_thaliana_hap2.fa")
out  = Path("data/A_thaliana_hap1_vs_hap2.phased.vcf")
out.parent.mkdir(exist_ok=True)

h1 = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(hap1), "fasta")}
h2 = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(hap2), "fasta")}

with open(out, "w") as f:
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for chrom in sorted(set(h1) & set(h2)):
        seq1, seq2 = h1[chrom], h2[chrom]
        for i, (a, b) in enumerate(zip(seq1, seq2), start=1):
            if a != b and a in "ACGT" and b in "ACGT":
                f.write(f"{chrom}\t{i}\t.\t{a}\t{b}\t.\tPASS\t.\tGT\t0|1\n")

