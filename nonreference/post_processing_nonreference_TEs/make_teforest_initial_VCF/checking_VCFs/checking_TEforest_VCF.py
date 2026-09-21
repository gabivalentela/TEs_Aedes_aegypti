import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../merging_VCFs/TEforest_VCF_nonreference_TEs.csv")
base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'score',
       'strand']
genotype_columns = [col for col in df.columns if col not in base_columns]
print('the number of genome columns:')
print(len(genotype_columns))

# count the number of TEs in total
print('the number of TEs in total:')
print(df.shape[0])

# number of TEs per genome - count the number of 1/1 in each column
df_genotypes = df[genotype_columns]
te_counts_per_genome = (df_genotypes == '1/1').sum()
print('the number of TEs per genome:')
te_counts_per_genome.to_csv('TE_counts_per_genome.csv', header=True)
plt.figure(figsize=(15, 6))
te_counts_per_genome.plot(kind='bar')
plt.xlabel('Genome')
plt.ylabel('Number of TEs')
plt.title('Number of TEs per Genome')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('TE_counts_per_genome.png')
plt.close()
