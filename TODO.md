
## TODO
### Tests to validate our tool:
- test1: run simulate with larger dataset (50M reads x 100Kb per chrom)  (best)
- test2: run simulate with larger dataset (1M reads x 100Kb per chrom)
- test3: run simulate with larger dataset (50M reads x 15Kb per chrom)
- test4: run simulate with larger dataset (1M reads x 15Kb per chrom)    (worst)
- test5: run simulate with true recomb map of Chr? of At. Show that we can hit the true result.
### Tests on empirical Amaranthus data:
- test6: run on At pollen pool data. Plot histogram that also shows centromeres, telomeres, and SDR.
- test7: run on Ap pollen pool data. Plot histogram that also shows centromeres, telomeres, and SDR.


## Installation

Set up and activate a new conda environment
```bash
conda create -f src/environment.yml -n disto
conda activate disto
```

Install the local `disto` package into the conda env
```bash
cd Distortopia-demo/
pip install -e . --no-deps
```

Call the `disto` program top-level CLI
```bash
disto -h
```

## Running an example

Simulate 10M gamete reads for chr 1 of a REF genome.
```bash
disto simulate -r REF.fa -c 1 -n 10_000_000 -s 123 -p test
```

Map reads, call and phase variants.
```bash
disto mapcall -r REF.fa -g test.gametes.fastq.gz
```

Infer crossovers
```bash
disto infer -r REF.fa -b test.sorted.bam -v test.phased.vcf.gz
```

Plot crossovers
```bash
disto plot -t test.tsv
```

