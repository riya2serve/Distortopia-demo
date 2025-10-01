
# Vision

```bash
tool -b BAM -r REF -o OUT
```

## Dependencies
- minimap2
- bcftools
- whatshap
- matplotlib
- pandas
- numpy

## Steps
1. Call variant positions in BAM file using bcftools => .vcf
2. Filter VCF to keep only high quality hetero sites using bcftools => .filtered.vcf
3. Phase variants in whatshap => .phased.filtered.vcf
4. Write TSV of read phases => .tsv
5. Analyse TSV and plot histogram..

## Step 2
- Exclude INDELS
- Exclude MNPs (variants should be bi-allelic)
- Exclude Low Quality variants (Q>30)
- Require variant call to be heterozygous.

## Step 3
- Write Python function(VCF, BAM)
- Iterate over reads in bam and get the following information:
  - [Scaff, start, end, nsnps, phased-snps (e.g., 00000), crossover-left, crossover-right]
  - e.g., [Chr1, 5000, 55000, 10, 0000000000, NA, NA]
  - e.g., [Chr1, 8000, 62000, 8, 00001111, 25000, 35000]
- SUBSTEPS:
  - Parse read chrom, start, end from BAM
  - Get number of variants in this region from VCF (e.g., 10)
  - Get positions of each variant (e.g., [10, 20, 200, 3000, 6000, 7000, 10000, 20000, ...]
  - Get phase vector of these positions (e.g., [0000001111]
  - If one cross-over occurred store position on either side of it: (7000, 10000)
  - If no cross-over then store (NA, NA)
  - If >1 cross-over then exclude read as bad variant calls.
  - Write results to TSV
 
  
