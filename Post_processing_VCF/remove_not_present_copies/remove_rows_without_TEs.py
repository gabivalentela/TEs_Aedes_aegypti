import pandas as pd

def remove_rows_without_TEs(df):
    mask = (df.iloc[:, 5:] == '1/1').any(axis=1)
    df = df[mask]
    return df

non_reference_VCF = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/nonreference/VCF_nonreference_all_data.csv'
reference_VCF = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/reference/subset/VCF_TEs_reference_subset_100000_strand.csv'

df_nonref = pd.read_csv(non_reference_VCF)
df_ref = pd.read_csv(reference_VCF)
# ignore first column
df_ref = df_ref.iloc[:, 1:]

df_nonref_updated = remove_rows_without_TEs(df_nonref)
# print results
print('Number of rows in original nonreference VCF:', df_nonref)
print('Number of rows in non-reference VCF after removing rows without TEs:', df_nonref_updated)

df_ref_updated = remove_rows_without_TEs(df_ref)
# print results
print('Number of rows in original reference VCF:', df_ref)
print('Number of rows in reference VCF after removing rows without TEs:', df_ref_updated)