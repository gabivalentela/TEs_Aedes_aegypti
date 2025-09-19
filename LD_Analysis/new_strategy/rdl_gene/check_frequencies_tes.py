import pandas as pd
import ast  # to safely parse the list-like strings

ref_df = './genotype_info_Reference.csv'
non_ref_df = './genotype_info_Non-Reference.csv'

df_ref = pd.read_csv(ref_df)
print(df_ref)
df_nonref = pd.read_csv(non_ref_df)

# Function to safely convert string representation of lists to actual lists
def safe_literal_eval(val):
    if isinstance(val, str):
        try:
            return ast.literal_eval(val)
        except (ValueError, SyntaxError) as e:
            print(f"Error parsing: {val[:50]}... Error: {e}")
            return []
    else:
        return val

# Convert the string representation of lists into actual lists
df_ref['Genotypes'] = df_ref['Genotypes'].apply(safe_literal_eval)
df_nonref['Genotypes'] = df_nonref['Genotypes'].apply(safe_literal_eval)

# Verify all genotypes are now lists and check lengths
print("Reference genotype lengths:")
print(df_ref['Genotypes'].apply(len))
print("Non-reference genotype lengths:")
print(df_nonref['Genotypes'].apply(len))

# ---- Option 2: Calculate frequency of "1/1" per TE ----
df_ref['freq_1_1'] = df_ref['Genotypes'].apply(lambda g: g.count('1/1') / len(g) if len(g) > 0 else 0)
print(df_ref[['TE', 'freq_1_1']].head())

df_nonref['freq_1_1'] = df_nonref['Genotypes'].apply(lambda g: g.count('1/1') / len(g) if len(g) > 0 else 0)
print(df_nonref[['TE', 'freq_1_1']].head())