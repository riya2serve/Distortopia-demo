import pandas as pd
import numpy as np
import ast

# --- PATHS: edited as needed ------------------------------------
# COmapper output
co_file = "/Users/riyarampalli/hacks/COmapper/test1_COmapper.csv"

# Distortopia infer output TSV (from disto infer)
dist_file = "/Users/riyarampalli/hacks/Distortopia-demo/crossovers/test1.tsv"
# ---------------------------------------------------------------------


# 1. Load Distortopia crossovers and keep only real Chr3 CO windows
dist = pd.read_csv(dist_file, sep="\t")

# Distortopia columns: scaff start end nsnps phased_snps crossover_left crossover_right read
dist_chr3 = dist[(dist["scaff"] == "Chr3") &
                 (dist["crossover_left"] != "NA") &
                 (dist["crossover_right"] != "NA")].copy()

# Making sure left/right are numeric
dist_chr3["crossover_left"] = dist_chr3["crossover_left"].astype(float)
dist_chr3["crossover_right"] = dist_chr3["crossover_right"].astype(float)

print("Distortopia Chr3 crossovers:", len(dist_chr3))


# 2. Loading COmapper output and extracting Chr3 COs
co = pd.read_csv(co_file)

# COmapper columns: qname, chr_num, pos, length, num_of_snp, snp_pos, bases,
# types, SNP_genotype, haplotype, CO_info
# Keep only Chr3 and non-parental haplotypes
co_chr3 = co[(co["chr_num"] == 3) & (co["haplotype"] != 0)].copy()
print("COmapper Chr3 non-parental reads:", len(co_chr3))


# 3. Parsing CO_info to get left/right/mid/width
def parse_co_info(x):
    """
    CO_info is like [chr_num, lb, rb, co, width] or 0.
    Return a series of (lb, rb, mid, width).
    """
    if pd.isna(x) or x == 0 or x == "0":
        return pd.Series([np.nan, np.nan, np.nan, np.nan],
                         index=["co_left", "co_right", "co_mid", "co_width"])
    try:
        arr = ast.literal_eval(str(x))
        if not isinstance(arr, (list, tuple)) or len(arr) < 5:
            return pd.Series([np.nan, np.nan, np.nan, np.nan],
                             index=["co_left", "co_right", "co_mid", "co_width"])
        _, lb, rb, mid, width = arr
        return pd.Series([float(lb), float(rb), float(mid), float(width)],
                         index=["co_left", "co_right", "co_mid", "co_width"])
    except Exception:
        return pd.Series([np.nan, np.nan, np.nan, np.nan],
                         index=["co_left", "co_right", "co_mid", "co_width"])


co_chr3 = co_chr3.join(co_chr3["CO_info"].apply(parse_co_info))

# Keep only reads with a valid CO interval
co_chr3_valid = co_chr3[~co_chr3["co_mid"].isna()].copy()
print("COmapper Chr3 reads with valid CO_info:", len(co_chr3_valid))


# 4. Define an overlap rule
# Let a COmapper CO "match" Distortopia if its midpoint
# falls inside any Distortopia [crossover_left, crossover_right] window.
def midpoint_in_any_dist_window(mid, dist_df):
    if np.isnan(mid):
        return False
    return ((dist_df["crossover_left"] <= mid) &
            (dist_df["crossover_right"] >= mid)).any()


co_chr3_valid["overlaps_distortopia"] = co_chr3_valid["co_mid"].apply(
    lambda m: midpoint_in_any_dist_window(m, dist_chr3)
)

n_overlap = co_chr3_valid["overlaps_distortopia"].sum()
print(f"COmapper Chr3 COs whose midpoint falls in a Distortopia window: {n_overlap}")
print(f"Fraction of COmapper Chr3 COs matching Distortopia: {n_overlap / len(co_chr3_valid):.3f}")


# 5. Optionally also measure the reverse: how many Distortopia windows are hit
def dist_window_hit(row, co_df):
    # Check if any COmapper midpoint lies in this Distortopia window
    return ((co_df["co_mid"] >= row["crossover_left"]) &
            (co_df["co_mid"] <= row["crossover_right"])).any()


dist_chr3["hit_by_COmapper"] = dist_chr3.apply(
    dist_window_hit, axis=1, co_df=co_chr3_valid
)

n_dist_hit = dist_chr3["hit_by_COmapper"].sum()
print(f"Distortopia Chr3 windows hit by at least one COmapper CO midpoint: {n_dist_hit}")
print(f"Fraction of Distortopia Chr3 CO windows matched: {n_dist_hit / len(dist_chr3):.3f}")

