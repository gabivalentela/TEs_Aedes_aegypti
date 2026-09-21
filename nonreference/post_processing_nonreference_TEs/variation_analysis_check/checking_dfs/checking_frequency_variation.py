import pandas as pd
import os
import matplotlib.pyplot as plt

def calculate_frequency(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    df['count_11'] = df[genotype_columns].apply(lambda row: sum(row == '1/1'), axis=1)
    df['count_00'] = df[genotype_columns].apply(lambda row: sum(row == '0/0'), axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, col, color in [
        (axes[0], 'count_11', 'steelblue'),
        (axes[1], 'count_00', 'orange'),
    ]:
        freq = df[col].value_counts().sort_index()

        ax.bar(freq.index, freq.values, color=color, alpha=0.7, edgecolor='none')

        # Annotate min and max count values
        min_count, max_count = freq.index.min(), freq.index.max()
        for count_val, annot_color, offset in [
            (min_count, 'red',   (10, 5)),
            (max_count, 'green', (10, 5)),
        ]:
            ax.annotate(f'count={count_val}\nn={freq[count_val]}',
                        xy=(count_val, freq[count_val]),
                        xytext=offset, textcoords='offset points',
                        fontsize=8, color=annot_color,
                        arrowprops=dict(arrowstyle='->', color=annot_color, lw=1))

        ax.set_xlabel(f'{col} value')
        ax.set_ylabel('number of variants')
        ax.set_title(f'How often each {col} value appears')

    plt.tight_layout()
    plt.savefig('frequency_distribution.png', dpi=150)
    plt.close()

list_dfs = []
list_files = os.listdir('../../post_processing_overlap_df/')
for file in list_files:
    if file.endswith('.csv'):
        variation = file.split('_')[-1]
        df = pd.read_csv('../../post_processing_overlap_df/' + file)
        df['variation'] = variation
        list_dfs.append(df)

merged_df = pd.concat(list_dfs, ignore_index=True)
merged_df = merged_df.drop_duplicates(subset=['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], keep=False)
merged_df_1000 = merged_df[merged_df['variation'] == '1000.csv']
merged_df_1000 = merged_df_1000.drop(columns=['variation'])
calculate_frequency(merged_df_1000)