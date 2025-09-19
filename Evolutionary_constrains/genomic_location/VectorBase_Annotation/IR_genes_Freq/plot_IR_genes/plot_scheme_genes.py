import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

def exon_list_retrieve(annotation_vector_base, gene_code):
    # open gff file as pandas df
    gff_df = pd.read_csv(annotation_vector_base, sep='\t', comment='#', header=None)
    # get the rows where the column 2 has 'exon' 
    exon_df = gff_df[gff_df[2] == 'exon']
    # get the rows where the column 8 has the gene name == AAEL023266
    gene_exon_df = exon_df[exon_df[8].str.contains(gene_code)]
    # get information of exon name from last colum after ID and before the first semicolon
    gene_exon_df['exon_info'] = gene_exon_df[8].str.extract(r'ID=([^;]+)')[0].tolist()
    # make a list of dictionaries with the exon information
    exon_list = []
    for index, row in gene_exon_df.iterrows():
        exon_info = {
            'exon_name': row['exon_info'],
            'start_position': row[3],
            'end_position': row[4]
        }
        exon_list.append(exon_info)

    return exon_list

#annotation_vector_base = '../../../../../Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
#annotation_vector_base = '../../../../../../final_scripts/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
annotation_vector_base = '../../../../../Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
# Gene and TE information
gene_sizes = {}
gene_sizes['CYP6P12'] = {'code_gene':'AAEL014891', 'start':271403322, 'end':271406949, 'te_name': 'rnd_1_family_337_DNA_Sola', 'te_start': 271402372, 'te_end': 271402953}
gene_sizes['GSTD11'] = {'code_gene':'AAEL010582', 'start':301237372, 'end':301238287, 'te_name': 'rnd_1_family_844_Unknown', 'te_start': 301236516, 'te_end': 301236900}
gene_sizes['GSTZ1'] = {'code_gene':'AAEL011934', 'start':286326225, 'end':286447671, 'te_name': 'rnd_1_family_696_Unknown', 'te_start': 286419545, 'te_end': 286420010}

# Margin information
margin_size = 1000

# Colors for genes and TEs
gene_color = 'lightpink'  # Green for genes
te_colors = {'rnd_1_family_337_DNA_Sola': '#FF5722', 
             'rnd_1_family_844_Unknown': '#2196F3', 
             'rnd_1_family_696_Unknown': '#FF9800'}

# Create individual plots for each gene-TE combination
for gene_name, gene_info in gene_sizes.items():
    # Create a new figure for each gene
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Calculate plot boundaries with margin
    plot_start = min(gene_info['start'], gene_info['te_start']) - margin_size
    plot_end = max(gene_info['end'], gene_info['te_end']) + margin_size
    
    # Y position for this single gene plot
    y_pos = 0.5
    
    # Draw margin area (light gray background)
    margin_rect = patches.Rectangle(
        (plot_start, y_pos - 0.3), 
        plot_end - plot_start, 
        0.6, 
        linewidth=1, 
        edgecolor='lightgray', 
        facecolor='lightgray',
        alpha=0.3
    )
    ax.add_patch(margin_rect)
    
    # Draw gene rectangle
    gene_rect = patches.Rectangle(
        (gene_info['start'], y_pos - 0.2), 
        gene_info['end'] - gene_info['start'], 
        0.4, 
        linewidth=2, 
        edgecolor='black', 
        facecolor=gene_color,
        alpha=0.8
    )
    ax.add_patch(gene_rect)
    
    # Draw TE rectangle
    te_color = te_colors.get(gene_info['te_name'], '#9C27B0') 
    te_rect = patches.Rectangle(
        (gene_info['te_start'], y_pos - 0.15), 
        gene_info['te_end'] - gene_info['te_start'], 
        0.3, 
        linewidth=2, 
        edgecolor=te_color, 
        facecolor=te_color,
        alpha=0.8
    )
    ax.add_patch(te_rect)
    
    # Retrieve exon list for the current gene
    exon_list = exon_list_retrieve(annotation_vector_base, gene_info['code_gene'])
    print(f"Exon list for {gene_name}: {exon_list}")
    # Draw each exon as a black rectangle
    for exon in exon_list:
        exon_rect = patches.Rectangle(
            (exon['start_position'], y_pos - 0.1), 
            exon['end_position'] - exon['start_position'],
            0.2,  
            linewidth=1,
            edgecolor='black',
            facecolor='black',
            alpha=1
        )
        ax.add_patch(exon_rect)

    # Add gene name label
    ax.text(gene_info['start'], y_pos + 0.35, f"Gene: {gene_name}", 
            fontsize=14, fontweight='bold', ha='left')
    
    ax.text(
        gene_info['te_end'],       # anchor at TE end
        y_pos - 0.45, 
        f"TE: {gene_info['te_name']}", 
        fontsize=12, 
        ha='right',   # text ends at this x-coordinate
        style='italic'
    )
    
    # Add coordinates as text
    #gene_coords = f"Gene coordinates: {gene_info['start']:,} - {gene_info['end']:,} ({gene_info['end'] - gene_info['start']:,} bp)"
    #te_coords = f"TE coordinates: {gene_info['te_start']:,} - {gene_info['te_end']:,} ({gene_info['te_end'] - gene_info['te_start']:,} bp)"
    
    # Position coordinate text at the bottom
    #ax.text(plot_start + (plot_end - plot_start) * 0.02, y_pos - 0.65, gene_coords, 
    #        fontsize=10, ha='left')
    #ax.text(plot_start + (plot_end - plot_start) * 0.02, y_pos - 0.75, te_coords, 
    #       fontsize=10, ha='left')

    # Set plot properties
    ax.set_ylim(-0.9, 1.0)
    ax.set_xlim(plot_start, plot_end)
    ax.set_xlabel('Genomic Position (bp)', fontsize=12)
    ax.set_title(f'{gene_name} and Associated Transposable Element', fontsize=14, fontweight='bold')
    
    # Remove y-axis ticks and labels since they're not meaningful
    ax.set_yticks([])
    
    # Format x-axis to show coordinates in a readable format
    ax.ticklabel_format(style='plain', axis='x')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Add legend
    legend_elements = [
        patches.Patch(color=gene_color, label=f'Gene ({gene_name})'),
        patches.Patch(color=te_color, label=f'TE ({gene_info["te_name"]})'),
        patches.Patch(color='black', label='Exons'),
        patches.Patch(color='lightgray', alpha=0.3, label=f'intergenic region ({margin_size:,} bp)')
    ]

    # make the legend be outside the plot
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), fontsize=10)    
    
    # Add grid and adjust layout
    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    
    # Save the plot
    filename = f"{gene_name}_TE_plot.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {filename}")
    # Close the figure to free memory
    plt.close()

# Print summary information
print("Gene and TE Summary:")
print("-" * 50)
for gene_name, info in gene_sizes.items():
    gene_length = info['end'] - info['start']
    te_length = info['te_end'] - info['te_start']
    print(f"{gene_name}:")
    print(f"  Gene length: {gene_length:,} bp")
    print(f"  TE length: {te_length:,} bp")
    print(f"  TE name: {info['te_name']}")
    print()