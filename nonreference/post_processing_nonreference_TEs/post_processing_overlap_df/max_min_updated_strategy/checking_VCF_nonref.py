import pandas as pd
import sys

def remove_singletons(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    # count number of 1/1 
    df['count_1/1'] = df[genotype_columns].apply(lambda x: x.str.count('1/1').sum(), axis=1)
    # remove rows with count_1/1 = 1
    df = df[df['count_1/1'] != 1]
    # drop the count_1/1 column
    df = df.drop(columns=['count_1/1'])
    return df

def remove_fixed_almost_fixed_df(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    print(f'the number of genomes is {len(genotype_columns)}')
    # count number of 1/1 
    df['count_1/1'] = df[genotype_columns].apply(lambda x: x.str.count('1/1').sum(), axis=1)
    df['count_0/0'] = df[genotype_columns].apply(lambda x: x.str.count('0/0').sum(), axis=1)
    # remove rows with count_1/1 = 122 and 121
    df = df[df['count_1/1'] != 122]
    df = df[df['count_1/1'] != 121]
    df = df[df['count_0/0'] != 122]
    # drop the count_1/1 column
    df = df.drop(columns=['count_1/1', 'count_0/0'])
    return df

#path_VCF = '../combine_retroseq_TEforest/matching_positions_retroseq_teforest.csv'
path_VCF = sys.argv[1]
variation = sys.argv[2]
df_VCF = pd.read_csv(path_VCF)
# remove these columns 
#columns_to_remove = ['Unnamed: 0', 'Unnamed: 5',
#       'position_start_forest.2', 'position_start_forest.3',
#       'position_end_forest.2', 'position_end_forest.3', 'score_forest', 'position_start_forest', 'position_start_forest.1',
#       'position_end_forest', 'position_end_forest.1']
#df_VCF = df_VCF.drop(columns=columns_to_remove)
# sort the TEs based on chromosome, position start and position end
df_VCF = df_VCF.sort_values(by=['chromosome', 'position_start', 'position_end'])
# remove singletons
no_singletons_df = remove_singletons(df_VCF)
# remove fixed and almost fixed
final_df = remove_fixed_almost_fixed_df(no_singletons_df)
# save the final dataframe to a new csv file
final_df.to_csv(f'final_df_overlap_retroseq_teforest_nonreference_{variation}.csv', index=False)