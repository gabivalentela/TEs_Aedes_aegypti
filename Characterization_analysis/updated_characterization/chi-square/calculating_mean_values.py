import pandas as pd

reference_file = './chi_square_results_reference.csv'
nonreference_file = './chi_square_results_nonreference.csv'

df_ref = pd.read_csv(reference_file)
df_nonref = pd.read_csv(nonreference_file)

mean_values_ref = df_ref['p-value'].mean()
mean_values_nonref = df_nonref['p-value'].mean()

# print the mean values with all digits
print(f'Mean p-value for reference: {mean_values_ref}')
print(f'Mean p-value for non-reference: {mean_values_nonref}')
