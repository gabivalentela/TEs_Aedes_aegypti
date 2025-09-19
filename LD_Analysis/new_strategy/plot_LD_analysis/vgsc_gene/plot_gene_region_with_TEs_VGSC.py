import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import re
import numpy as np
import pandas as pd

def plot_elements_and_exons(list_of_elements, exon_info, gene, output_file):
    fig, ax = plt.subplots(figsize=(20, 10))
    
    start = gene['start_position']
    end = gene['end_position']

    width = end - start
    height = 1
    gene_rect = patches.Rectangle((start, 1), width, height, linewidth=1, edgecolor='lightblue', facecolor='lightblue')
    ax.add_patch(gene_rect)

    exon_color = 'black'

    for idx, exon in enumerate(exon_list):
        exon_start = exon['start_position']
        exon_end = exon['end_position']
        exon_width = exon_end - exon_start
        exon_rect = patches.Rectangle((exon_start, 1), exon_width, height, linewidth=1, edgecolor='black', facecolor=exon_color)
        ax.add_patch(exon_rect)
        ax.text(exon_start + exon_width / 2, 2.5, f"exon_{idx+1}", ha='center', va='bottom', fontsize=8, rotation=90)

    colors = plt.cm.tab20b(np.linspace(0, 1, len(list_of_elements)))

    num_heights = 6  
    for idx, TE in enumerate(list_of_elements):
        te_start = TE['start_position']
        te_end = TE['end_position']
        te_width = te_end - te_start
        te_color = colors[idx]
        TE_rectangle = patches.Rectangle((te_start, 1), te_width, height, linewidth=1, edgecolor=te_color, facecolor=te_color)
        ax.add_patch(TE_rectangle)
        
        arrow_height = 0.3 + (idx % num_heights) * 0.11  
        ax.annotate(TE['TE_name'], xy=(te_start + te_width / 2, arrow_height), xytext=(te_start + te_width / 2, arrow_height - 0.3),
                    arrowprops=dict(facecolor=te_color, shrink=0.05), ha='center', va='top', fontsize=8)
        
        if 'rsquared' in TE:
            rsquared_text = f"R²={TE['rsquared']:.2f}"
            ax.text(te_start + te_width / 2, arrow_height + 0.1, rsquared_text, ha='center', va='bottom', fontsize=8, color=te_color)

    plt.xlim(start - 130500, end + 100000)
    plt.ylim(0, num_heights * 0.5)
    plt.xlabel('Position')
    plt.ylabel(' ')
    plt.title('Gene Region and Exons Plot')

    handles = [patches.Patch(color=colors[idx], label=TE['TE_name']) for idx, TE in enumerate(list_of_elements)]
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.15, 1))

    plt.savefig(output_file)
    plt.close()

def exon_list_retrieve(annotation_vector_base):
    gff_df = pd.read_csv(annotation_vector_base, sep='\t', comment='#', header=None)
    exon_df = gff_df[gff_df[2] == 'exon']
    gene_exon_df = exon_df[exon_df[8].str.contains('AAEL023266')]
    gene_exon_df['exon_info'] = gene_exon_df[8].str.extract(r'ID=([^;]+)')[0].tolist()
    exon_list = []
    for index, row in gene_exon_df.iterrows():
        exon_info = {
            'exon_name': row['exon_info'],
            'start_position': row[3],
            'end_position': row[4]
        }
        exon_list.append(exon_info)

    return exon_list

annotation_vector_base = '../../../../Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
exon_list = exon_list_retrieve(annotation_vector_base)

DNA_CMC_Chapaev_ref = {'TE_name': 'rnd_1_family_56_DNA_CMC_Chapaev', 'start_position': 316003069, 'end_position': 316003804, 'rsquared': 0.7124999999999999}

V1016I = {'TE_name': 'V1016I', 'start_position': 315983763, 'end_position': 315983764} 
I915K = {'TE_name': 'I915K', 'start_position': 315999297, 'end_position': 315999298} 
S723T = {'TE_name': 'S723T', 'start_position': 316014588, 'end_position': 316014589} 
V410L = {'TE_name': 'V410L', 'start_position': 316080722, 'end_position': 316080723} 

gene_region_VGSC = {'gene_name': 'VGSC', 'start_position': 315926360, 'end_position': 316405639}

list_of_elements = [DNA_CMC_Chapaev_ref, V1016I, I915K, S723T, V410L]

output_file = 'vgsc_plot_exons_TEs.svg'
plot_elements_and_exons(list_of_elements, exon_list, gene_region_VGSC, output_file)
