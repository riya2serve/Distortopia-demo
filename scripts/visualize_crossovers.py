#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def plot_crossovers(tsv_path, species_name, outdir="plots"):
    df = pd.read_csv(tsv_path, sep="\t")
    Path(outdir).mkdir(exist_ok=True)

    # Convert crossover positions to numeric (NaN for 'NA')
    df["crossover_left"] = pd.to_numeric(df["crossover_left"], errors="coerce")
    df["crossover_right"] = pd.to_numeric(df["crossover_right"], errors="coerce")

    # 1️⃣ Histogram of crossover counts per scaffold
    counts = df.groupby("scaff")["crossover_left"].apply(lambda x: x.notna().sum()).reset_index()
    plt.figure(figsize=(8, 5))
    sns.barplot(x="scaff", y="crossover_left", data=counts, color="skyblue")
    plt.title(f"{species_name}: Number of crossovers per scaffold")
    plt.xlabel("Scaffold")
    plt.ylabel("Crossover count")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{outdir}/{species_name}_crossovers_per_scaffold.png", dpi=300)
    plt.close()

    # 2️⃣ Genome-wide crossover density plot
    plt.figure(figsize=(10, 5))
    sns.histplot(df["crossover_left"].dropna(), bins=50, color="coral", alpha=0.7)
    plt.title(f"{species_name}: Distribution of crossover breakpoints")
    plt.xlabel("Genomic position (bp)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"{outdir}/{species_name}_crossover_density.png", dpi=300)
    plt.close()

    # 3️⃣ Read-level crossover summary
    df["has_CO"] = df["crossover_left"].notna()
    total_reads = len(df)
    n_co = df["has_CO"].sum()
    print(f"{species_name}: {n_co}/{total_reads} reads ({n_co/total_reads:.2%}) show a crossover")

if __name__ == "__main__":
    plot_crossovers("data/A_thaliana_gamete_crossovers.tsv", "A. thaliana")
    plot_crossovers("data/A_lyrata_gamete_crossovers.tsv", "A. lyrata")

