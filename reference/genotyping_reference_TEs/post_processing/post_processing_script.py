import sys
import pandas as pd

def remove_small_TEs(df):
    threshold = 150
    df['size'] = df['position_end'] - df['position_start']
    df = df[df['size'] > threshold]
    # drop the 'size' column after filtering
    df.drop(columns=['size'], inplace=True)
    return df

def remove_high_freq_TEs_or_all_missing(df):
    frequency_threshold = 0.95
    base_columns = ['chromosome', 'position_start', 'position_end', 'strand', 'TE_name', 'score']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    print(f'the number of genotype columns is: {len(genotype_columns)}')
    df['count_11'] = (df[genotype_columns] == '1/1').sum(axis=1)
    df['count_01'] = (df[genotype_columns] == '0/1').sum(axis=1)
    df['count_00'] = (df[genotype_columns] == '0/0').sum(axis=1)
    df['count_missing'] = (df[genotype_columns] == './.').sum(axis=1)
    # remove if all genotypes are missing
    df = df[df['count_missing'] < len(genotype_columns)]
    # calculate frequency of TEs
    df['total_count'] = df['count_11'] + df['count_01'] + df['count_00']
    df['frequency_1'] = (df['count_11']*2 + df['count_01']) / (df['total_count'] * 2)
    df['frequency_0'] = (df['count_00']*2 + df['count_01']) / (df['total_count'] * 2)
    # remove TEs that are fixed for absence (all genotypes are 0/0, frequency_0 == 1)
    df = df[df['frequency_0'] != 1]
    # remove TEs that have frequency_1 > frequency_threshold
    df = df[df['frequency_1'] <= frequency_threshold]
    # drop the temporary count and frequency columns
    df.drop(columns=['count_11', 'count_01', 'count_00', 'count_missing', 'total_count', 'frequency_1', 'frequency_0'], inplace=True)
    return df

file_path = sys.argv[1]
file_name = file_path.split('/')[-1]
df = pd.read_csv(file_path)
print(f'original number of rows in the dataframe: {len(df)}')
# remove small TEs
df_no_small_TEs = remove_small_TEs(df)
print(f'number of rows after removing small TEs: {len(df_no_small_TEs)}')
# remove TEs with frequency above 0.95
df_no_small_TEs = remove_high_freq_TEs_or_all_missing(df_no_small_TEs)
print(f'number of rows after removing high frequency TEs: {len(df_no_small_TEs)}')
# make a csv file with the filtered dataframe
output_file_path = file_name.replace('.csv', '_filtered.csv')
df_no_small_TEs.to_csv(output_file_path, index=False)
print(f'filtered dataframe saved to: {output_file_path}')

# now make removing_high_frequency keeping small TEs
df = remove_high_freq_TEs_or_all_missing(df)
print(f'number of rows after removing high frequency TEs (keeping small TEs): {len(df)}')
output_file_path_small_TEs = file_name.replace('.csv', '_filtered_keep_small_TEs.csv')
df.to_csv(output_file_path_small_TEs, index=False)
print(f'filtered dataframe (keeping small TEs) saved to: {output_file_path_small_TEs}')

