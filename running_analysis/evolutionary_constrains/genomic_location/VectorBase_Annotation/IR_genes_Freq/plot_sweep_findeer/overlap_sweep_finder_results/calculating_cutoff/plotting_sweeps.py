import pandas as pd
import matplotlib.pyplot as plt
import re
import os

general_path = '/overflow/dschridelab/users/rrlove/aedes/results/sweepfinder/outputs/merged/'
list_files = os.listdir(general_path)

information_sweeps = {}
for file in list_files:
    if 'no_relatives' in file and file.startswith('AaegL5'):
        chromosome = re.search(r'(AaegL5_\d)_no_relatives_(.+)_raw', file).group(1)
        country = re.search(r'(AaegL5_\d)_no_relatives_(.+)_raw', file).group(2)
        if 'trimmed' in file:
            full_path = os.path.join(general_path, file)
            information_sweeps[(chromosome, country, 'trimmed')] = full_path
        else:
            full_path = os.path.join(general_path, file)
            information_sweeps[(chromosome, country, 'not_trimmed')] = full_path

# make df with chromosomes - AaegL5_1, AaegL5_2, AaegL5_3 as columns and countries as rows
df = pd.DataFrame(
    columns=['AaegL5_1', 'AaegL5_2', 'AaegL5_3'],
    index = ['Brazil', 'Colombia', 'USA', 'Senegal', 'Gabon', 'Kenya']
)

# Define colors for each country
country_colors = {
    'Brazil': '#FF6B6B',
    'Colombia': '#4ECDC4',
    'USA': '#45B7D1',
    'Senegal': '#FFA07A',
    'Gabon': '#98D8C8',
    'Kenya': '#C7A4FF'
}

# Organize data by chromosome
chromosomes = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
data_by_chromosome = {chrom: {} for chrom in chromosomes}
position_limits = {'AaegL5_1': (12993807, 310821866), 'AaegL5_2': (495149, 467999272), 'AaegL5_3': (1749994, 393999083)}

# Collect all data first
for key, values in information_sweeps.items():
    chromosome, country, file_type = key
    if file_type == 'trimmed':
        trimmed = pd.read_csv(values, sep='\t')
        not_trimmed = pd.read_csv(information_sweeps[(chromosome, country, 'not_trimmed')], sep='\t')
        
        #first_position = trimmed.iloc[0]['location']
        #last_position = trimmed.iloc[-1]['location']
        first_position = position_limits[chromosome][0]
        last_position = position_limits[chromosome][1]
        
        data_by_chromosome[chromosome][country] = {
            'trimmed': trimmed,
            'not_trimmed': not_trimmed,
            'first_position': first_position,
            'last_position': last_position
        }
        
        df.loc[country, chromosome] = f'{first_position}-{last_position}'

# Create plots for each chromosome - all countries together
for chromosome in chromosomes:
    if data_by_chromosome[chromosome]:
        # Full chromosome plot
        plt.figure(figsize=(14, 6))
        for country, data in data_by_chromosome[chromosome].items():
            plt.plot(data['not_trimmed']['location'], data['not_trimmed']['LR'], 
                    color=country_colors[country], alpha=0.4, linewidth=1.5, label=country)
            plt.axvline(x=data['first_position'], color=country_colors[country], 
                       linestyle='--', alpha=0.6, linewidth=1)
            plt.axvline(x=data['last_position'], color=country_colors[country], 
                       linestyle='--', alpha=0.6, linewidth=1)
        
        plt.title(f'SweepFinder2 Results - {chromosome} (All Countries)')
        plt.xlabel('Genomic Position')
        plt.ylabel('Likelihood Ratio (LR)')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f'sweep_plots_new/{chromosome}_all_countries.png', dpi=300)
        plt.close()
        
        # Beginning of chromosome (0 to first trim + 10000)
        plt.figure(figsize=(14, 6))
        max_start = max([data['first_position'] for data in data_by_chromosome[chromosome].values()])
        
        for country, data in data_by_chromosome[chromosome].items():
            mask = data['not_trimmed']['location'] <= (max_start + 10000)
            subset = data['not_trimmed'][mask]
            plt.plot(subset['location'], subset['LR'], 
                    color=country_colors[country], alpha=0.5, linewidth=2, label=country)
            plt.axvline(x=data['first_position'], color=country_colors[country], 
                       linestyle='--', alpha=0.7, linewidth=1.5)
        
        plt.title(f'SweepFinder2 Results - {chromosome} (Chromosome Start)')
        plt.xlabel('Genomic Position')
        plt.ylabel('Likelihood Ratio (LR)')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f'sweep_plots_new/{chromosome}_all_countries_start.png', dpi=300)
        plt.close()
        
        # End of chromosome (last trim - 10000 to end)
        plt.figure(figsize=(14, 6))
        min_end = min([data['last_position'] for data in data_by_chromosome[chromosome].values()])
        
        for country, data in data_by_chromosome[chromosome].items():
            mask = data['not_trimmed']['location'] >= (min_end - 10000)
            subset = data['not_trimmed'][mask]
            plt.plot(subset['location'], subset['LR'], 
                    color=country_colors[country], alpha=0.5, linewidth=2, label=country)
            plt.axvline(x=data['last_position'], color=country_colors[country], 
                       linestyle='--', alpha=0.7, linewidth=1.5)
        
        plt.title(f'SweepFinder2 Results - {chromosome} (Chromosome End)')
        plt.xlabel('Genomic Position')
        plt.ylabel('Likelihood Ratio (LR)')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(f'sweep_plots_new/{chromosome}_all_countries_end.png', dpi=300)
        plt.close()

print('cutoff for each chromosome and country:')
print(df)