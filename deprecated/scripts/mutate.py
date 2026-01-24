import random
from Bio import SeqIO
from Bio.Seq import Seq

def mutate_fasta(input_fasta, output_fasta, snp_rate=0.001, indel_rate=0.0, seed: random.Random=123):
    rng = random.Random(seed)
    nts = ("A", "C", "G", "T")

    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            seq = list(str(rec.seq))
            i = 0
            while i < len(seq):
                base = seq[i].upper()

                # SNP mutation
                if base in nts and rng.random() < snp_rate:
                    seq[i] = rng.choice([x for x in nts if x != base])

                # INDEL mutation
                if indel_rate > 0 and rng.random() < indel_rate:
                    if rng.random() < 0.5:  # deletion
                        del seq[i]
                        continue
                    else:                   # insertion
                        seq.insert(i+1, rng.choice(nts))
                        i += 1

                i += 1

            rec.seq = Seq("".join(seq))   # FIX: wrap in Seq object
            SeqIO.write(rec, out, "fasta")

