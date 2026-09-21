import pandas as pd
import os

list_files = os.listdir('../')
for file in list_files:
    if file.endswith('.csv'):
        variation = file.split('_')[-1].split('.')[0]
        df = pd.read_csv(f'../{file}')
        # drop the frequency_diff
        df = df.drop(columns=['frequency_diff'])
        # remove duplicates - TE_information
        df = df.drop_duplicates(subset=['TE_information'], keep='first')
        # make csv with final df - without duplicates
        df.to_csv(f'final_df_nonreference_{variation}.csv', index=False)
        # make correlation matrix from this df
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Define country groups
        group1 = ['frequency_Kenya', 'frequency_Senegal', 'frequency_Gabon']
        group2 = ['frequency_USA', 'frequency_Brazil', 'frequency_Colombia']

        # Calculate mean frequency for each country in both groups
        group1_means = df[group1].mean()
        group2_means = df[group2].mean()


        # Remove 'frequency_' from the labels for plotting
        group1_labels = [country.replace('frequency_', '') for country in group1]
        group2_labels = [country.replace('frequency_', '') for country in group2]

        # Plot setup for mean frequencies
        plt.figure(figsize=(10, 6))

        # Plot mean frequencies for Group 1
        plt.bar(group1_labels, group1_means.values, color=['blue', 'green', 'orange'], label='Group 1 (Kenya, Senegal, Gabon)', alpha=0.7)

        # Plot mean frequencies for Group 2
        plt.bar(group2_labels, group2_means.values, color=['brown', 'purple', 'red'], label='Group 2 (USA, Brazil, Colombia)', alpha=0.7)

        # Labels and title
        plt.xlabel('Country', fontsize=12)
        plt.ylabel('Mean Frequency', fontsize=12)
        plt.title('Mean Reference TE Insertion Frequencies by Country', fontsize=14)
        plt.legend(title="Country Group")

        # Show plot
        plt.tight_layout()
        plt.savefig(f'mean_frequencies_nonreference_{variation}.svg')

        # Correlation calculation between all countries
        correlation_matrix = df[group1 + group2].corr()

        # Rename columns and index in the correlation matrix for plotting
        correlation_matrix.index = [x.replace('frequency_', '') for x in correlation_matrix.index]
        correlation_matrix.columns = [x.replace('frequency_', '') for x in correlation_matrix.columns]

        # Plot the correlation matrix as a heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', cbar=True, linewidths=0.5)
        plt.title('Correlation Matrix of reference TE Frequencies Across Countries', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'correlation_matrix_nonreference_{variation}.svg')

        # Print the correlation matrix for reference
        print("Correlation Matrix for All Countries:")
        print(correlation_matrix)
