import pandas as pd
import matplotlib.pyplot as plt

def classify_te_count(df):
    # Get genotype columns
    genotype_data = df.iloc[:, 6:]
    
    # Count number of "1/1" per TE
    count_1_1 = (genotype_data == "1/1").sum(axis=1)
    
    # Classify TEs into count-based categories
    df["count_category"] = pd.cut(
        count_1_1,
        bins=[1, 2, 3, 5, 7, float("inf")],  
        labels=["less than 2", "2 genomes", "At least 3", "At least 5", "More than 7"],
        include_lowest=True,  
        right=False  
    )

    # more than 7 TEs
    print(df[df["count_category"] == "More than 7"])
    # more than 5
    print(df[df["count_category"] == "At least 5"])
    # more than 3
    print(df[df["count_category"] == "At least 3"])
    # more than 2
    print(df[df["count_category"] == "2 genomes"])

    return df

def count_te_by_category(df):
    # Get genotype columns
    genotype_data = df.iloc[:, 6:-1]
    
    # Create a dictionary to store counts
    genome_te_counts = {"less than 2": {}, "2 genomes": {}, "At least 3": {}, "At least 5": {}, "More than 7": {}}
    
    # Iterate over each genome 
    for genome in genotype_data.columns:
        # Count "1/1" occurrences for each count category
        genome_te_counts["less than 2"][genome] = ((df[genome] == "1/1") & (df["count_category"] == "less than 2")).sum()
        genome_te_counts["2 genomes"][genome] = ((df[genome] == "1/1") & (df["count_category"] == "2 genomes")).sum()
        genome_te_counts["At least 3"][genome] = ((df[genome] == "1/1") & (df["count_category"] == "At least 3")).sum()
        genome_te_counts["At least 5"][genome] = ((df[genome] == "1/1") & (df["count_category"] == "At least 5")).sum()
        genome_te_counts["More than 7"][genome] = ((df[genome] == "1/1") & (df["count_category"] == "More than 7")).sum()
    
    return pd.DataFrame(genome_te_counts)

# Load data
nonreference_path = "/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes/merged_nonreference/VCF_nonreference_post_processing.csv"
df = pd.read_csv(nonreference_path)

# Classify TEs into count-based categories
df = classify_te_count(df)

# Count TEs per genome for each category
te_counts_by_genome = count_te_by_category(df)

# Plot stacked bar chart
te_counts_by_genome.plot(kind="bar", stacked=True, figsize=(40, 35), color=["blue", "orange", "red", "green", "purple"])
plt.xlabel("Genome")
plt.ylabel("TE Count")
plt.title("TE Counts per Genome by Presence in Genomes", fontsize=40)
plt.legend(title="TE Presence", fontsize=20)
plt.xticks(fontsize=8, rotation=90)  # Rotate x-axis labels for better readability
plt.savefig("genome_TE_counts_by_presence.png")
