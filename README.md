# Distortopia

**Distortopia: Simulation of long-read sequences, high-resolution mapping of meiotic crossovers, and identification of segregation distorters in *Arabidopsis thaliana* and *A. lyrata*.**

# Getting Started

First, clone this repository:

```bash 
git clone https://github.com/your-username/Distortopia-demo.git
cd Distortopia-demo
```
## Activate Conda Environment 
### Prerequisites
We recommend using the following conda environment to manage package dependencies. The following tools and libraries are required:

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://www.htslib.org/)
- [htslib](https://www.htslib.org/)
- [bcftools](https://samtools.github.io/bcftools/)
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Streamlit](https://streamlit.io/)
- [haptools](https://github.com/ajmazurie/haptools)

Required versions are specified in the [`environment.yml`](environment.yml). To create environment using provided YAML: 

```bash
conda env create --name distortopia --file environment.yml
```
## Generate Alternate Haploid Genomes 
1. Download TAIR9 and hardmasked reference genome sequences from [JGI Phytozome](https://phytozome-next.jgi.doe.gov). 
2. Create index file 
```bash

```
### Explanation of files in variants directory
[`Athaliana.fna`](Athaliana.fna), [`Alyrata.fna`](Alyrata.fna): parental genome assemblies in FASTA format

[`Athaliana.mmi`](Athaliana.mmi), [`Alyrata.mmi`](A.lyrata): minimap2 index files

### Input data
Distortopia requires long reads (e.g. Nanopore, Illumina, PacBio) in `.fq.gz`, `.fna`, (or `.fa`) formats. Uploaded files can be **up to** 2GB in size. 

Test input files are in the `raw_data/` directory. 

## Generate Long Reads From Gametes of Artificial Hybrids 

For the purposes of our simulation, these hybrids have uniform recombination rates. 

## Map Long Reads to Referene Genomes

## Call Variants 

## Run Crossover Detection Using [COMapper](https://github.com/KyuhaChoi-Lab/COmapper)








