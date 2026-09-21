import os
import pandas as pd

list_files = os.listdir("../")
combined_df = []
for file in list_files:
    if file.endswith('.vcf'):
        df = pd.read_csv(os.path.join("../", file), sep='\t')
        print(f'{file}: {len(df)}')
        combined_df.append(df)
combined_df = pd.concat(combined_df)
# print total number of rows in the combined dataframe
print(f'Total number of TEs in combined dataframe: {len(combined_df)}')
combined_df.to_csv('combined_reference_tes.vcf', sep='\t', index=False)