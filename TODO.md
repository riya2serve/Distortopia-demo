


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

