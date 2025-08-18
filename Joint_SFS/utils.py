import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def create_joint_sfs_subplot(ax, country_1, country_2, combined_freq_df):
    # Find the maximum frequency in each country separately
    max_freq_country1 = max(combined_freq_df[f'count_1/1_{country_1}'])
    max_freq_country2 = max(combined_freq_df[f'count_1/1_{country_2}'])
    
    # Initialize an empty grid using the correct max frequencies
    grid = np.zeros((max_freq_country1 + 1, max_freq_country2 + 1))
    
    # Fill the grid with counts of TEs for each pair of allele counts
    for i in range(max_freq_country1 + 1):
        for j in range(max_freq_country2 + 1): 
            grid[i, j] = np.sum(
                (combined_freq_df[f'count_1/1_{country_1}'] == i) & 
                (combined_freq_df[f'count_1/1_{country_2}'] == j)
            )
    
    # Create a logarithmic color normalization
    vmin = max(1, grid.min()) 
    vmax = grid.max()
    norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot on the given subplot axis
    cax = ax.imshow(grid, cmap='Spectral', origin='lower', norm=norm)
    ax.set_xlabel(f'{country_2} allele count', fontsize=15)
    ax.set_ylabel(f'{country_1} allele count', fontsize=15)
    ax.set_title(f'{country_1} vs {country_2}', fontsize=15)
    
    return cax
