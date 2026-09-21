import pandas as pd
import os

def calculate_frequency(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    df['count_11'] = df[genotype_columns].apply(lambda row: sum(row == '1/1'), axis=1)
    df['count_00'] = df[genotype_columns].apply(lambda row: sum(row == '0/0'), axis=1)
    # keep only rows where count_11 is 122 and count_00 is 0
    df_filtered = df[(df['count_11'] == 122) & (df['count_00'] == 0)]
    print(f"File: {file}, Count of rows with count_11 = 122 and count_00 = 0: {len(df_filtered)}")
    

list_files = os.listdir('../../post_processing_overlap_df/')
for file in list_files:
    if file.endswith('.csv'):
        df = pd.read_csv('../../post_processing_overlap_df/' + file)
        calculate_frequency(df)