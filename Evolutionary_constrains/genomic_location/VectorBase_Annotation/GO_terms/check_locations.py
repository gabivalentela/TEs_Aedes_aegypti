import pandas as pd
import mygene
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns

# Load nonref and ref datasets
nonref = '../genes_overlap_tes_nonref_no_duplicates.csv'
ref = '../genes_overlap_tes_ref_no_duplicates.csv'

df_nonref = pd.read_csv(nonref)
df_ref = pd.read_csv(ref)

# remove ID= from the gene_symbol column
df_nonref['gene_symbol'] = df_nonref['gene_symbol'].str.replace('ID=', '', regex=False)
df_ref['gene_symbol'] = df_ref['gene_symbol'].str.replace('ID=', '', regex=False)

# remove first column
df_nonref = df_nonref.iloc[:, 1:]
df_ref = df_ref.iloc[:, 1:]

# Get gene symbols from each DataFrame, dropping NaN values
gene_symbols_nonref = df_nonref['gene_symbol'].dropna().unique().tolist()
gene_symbols_ref = df_ref['gene_symbol'].dropna().unique().tolist()
print(gene_symbols_nonref)

# Initialize mygene
mg = mygene.MyGeneInfo()

def extract_go_terms(gene_list):
    go_dict = {}
    for gene in gene_list:
        try:
            # Query using AAEL IDs for Aedes aegypti (NCBI taxid: 7159)
            gene_info = mg.query(gene, scopes='symbol', fields='go', species=7159)
            if gene_info and 'hits' in gene_info and len(gene_info['hits']) > 0:
                hit = gene_info['hits'][0]
                go_terms = []
                if 'go' in hit and isinstance(hit['go'], dict):
                    for go_type in ['BP', 'CC', 'MF']:
                        if go_type in hit['go']:
                            terms = hit['go'][go_type]
                            if isinstance(terms, list):
                                go_terms.extend([term['term'] for term in terms])
                            elif isinstance(terms, dict):
                                go_terms.append(terms['term'])
                go_dict[gene] = "; ".join(go_terms) if go_terms else None
            else:
                go_dict[gene] = None
        except Exception as e:
            print(f"Error retrieving GO annotations for {gene}: {e}")
            go_dict[gene] = None
    return go_dict

# Extract GO terms for nonref and ref gene symbols
go_terms_dict_nonref = extract_go_terms(gene_symbols_nonref)
go_terms_dict_ref = extract_go_terms(gene_symbols_ref)

# Map GO terms back to DataFrames
df_nonref['GO_terms'] = df_nonref['gene_symbol'].map(go_terms_dict_nonref)
df_ref['GO_terms'] = df_ref['gene_symbol'].map(go_terms_dict_ref)

# Save the updated DataFrames
df_nonref.to_csv('./genes_overlap_tes_nonref_with_GO_noduplicate.csv', index=False)
df_ref.to_csv('./genes_overlap_tes_ref_with_GO_noduplicate.csv', index=False)

'''def plot_go_bar(go_terms_series, title, filename):
    # Ensure go_terms_series is a Pandas Series
    go_terms_series = pd.Series(go_terms_series)

    # Flatten GO terms: split semicolon-separated terms and filter out None/NaN values
    go_terms_list = [
        term.strip() for terms in go_terms_series.dropna() 
        if isinstance(terms, str)  # Ensure terms are strings before splitting
        for term in terms.split(";") if term
    ]

    # Count GO term occurrences
    go_count = Counter(go_terms_list)
    
    if not go_count:  # If no GO terms, exit early
        print(f"No GO terms available for {title}. Skipping plot.")
        return
    
    # Get the top 10 GO terms
    top_10 = go_count.most_common(10)
    labels, counts = zip(*top_10)
    
    # Compute proportions
    total = sum(counts)
    proportions = [count / total for count in counts]
    
    plt.figure(figsize=(30, 15))  # Adjusted for better visibility
    plt.bar(labels, proportions, color="darkblue")
    plt.yticks(fontsize=20)
    plt.xticks(rotation=45, fontsize=20)
    plt.ylabel('Proportion')
    plt.title(title)
    
    # Save plot to a file
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()  # Close the plot to avoid displaying it in interactive mode

# Extract GO terms for nonref and ref gene symbols
go_terms_nonref = extract_go_terms(gene_symbols_nonref)
go_terms_ref = extract_go_terms(gene_symbols_ref)

# Save GO term distributions to files
plot_go_bar(go_terms_nonref, 'Proportion of GO Terms for Nonref Genes', 'go_terms_nonref_bar_top10.png')
plot_go_bar(go_terms_ref, 'Proportion of GO Terms for Ref Genes', 'go_terms_ref_bar_top10.png')'''

def plot_go_table_grouped(go_terms_series, title, filename, color):
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    # Ensure go_terms_series is a Pandas Series
    go_terms_series = pd.Series(go_terms_series)

    # Flatten GO terms: split semicolon-separated terms and filter out None/NaN values
    go_terms_list = [
        term.strip() for terms in go_terms_series.dropna()
        if isinstance(terms, str)
        for term in terms.split(";") if term
    ]

    # Count GO term occurrences
    go_count = Counter(go_terms_list)

    if not go_count:
        print(f"No GO terms available for {title}. Skipping table.")
        return

    # Filter and group GO terms with at least 2 occurrences
    count_to_terms = {}
    for term, count in go_count.items():
        if count >= 2:
            count_to_terms.setdefault(count, []).append(term)

    if not count_to_terms:
        print(f"No GO terms found at least twice for {title}. Skipping table.")
        return

    # Sort by descending count
    grouped_data = sorted(count_to_terms.items(), key=lambda x: x[0], reverse=True)

    # Expand into DataFrame: one GO term per row, shared count per group
    expanded_rows = []
    for count, terms in grouped_data:
        terms = sorted(terms)
        for i, term in enumerate(terms):
            expanded_rows.append((count if i == 0 else "", term))

    df = pd.DataFrame(expanded_rows, columns=["Count", "GO Term"])

    # Create figure with better proportions for publication
    fig, ax = plt.subplots(figsize=(16, max(10, len(df) * 0.4)))
    ax.axis('tight')
    ax.axis('off')

    # Create table with adjusted column widths
    table = ax.table(cellText=df.values, colLabels=df.columns,
                     cellLoc='left', loc='center',
                     colWidths=[0.12, 0.88])

    # Enhanced styling for publication quality
    table.auto_set_font_size(False)
    table.set_fontsize(14)  # Larger font size
    table.scale(1.0, 2.0)   # Increased row height for better readability

    # Style header with professional look
    for i in range(len(df.columns)):
        table[(0, i)].set_facecolor(color)
        table[(0, i)].set_text_props(weight='bold', color='white', fontsize=16)
        table[(0, i)].set_height(0.08)

    # Remove all cell borders and apply subtle styling
    current_count = None
    group_color = 'white'  # Very light gray for grouping
    
    for i in range(1, len(df) + 1):
        row_count = df.iloc[i-1]['Count']
        
        # Determine if this is a new group
        if row_count != "":
            current_count = row_count
            use_group_color = True
        else:
            use_group_color = True
            
        for j in range(len(df.columns)):
            cell = table[(i, j)]
            
            # Remove all borders
            cell.set_linewidth(0)
            
            # Apply background color for grouping
            if use_group_color and current_count and int(current_count) % 2 == 0:
                cell.set_facecolor(group_color)
            else:
                cell.set_facecolor('white')
            
            # Increase font size for data cells
            cell.set_text_props(fontsize=14)
            
            # Align count column to center, GO terms to left
            if j == 0:  # Count column
                cell.set_text_props(ha='center', fontweight='bold' if row_count != "" else 'normal')
            else:  # GO Term column
                cell.set_text_props(ha='left')

    # Remove header cell borders
    for i in range(len(df.columns)):
        table[(0, i)].set_linewidth(0)

    # Professional title styling
    #plt.title(title, fontsize=18, fontweight='bold', pad=25, color='#2C3E50')

    # Adjust layout for publication
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    
    # Save with high DPI for publication quality
    plt.savefig(filename, bbox_inches='tight', dpi=300, facecolor='white', 
                edgecolor='none', pad_inches=0.2)
    plt.close()

    # Save to CSV with better formatting
    csv_filename = filename.replace('.png', '.csv')
    
    # Create a cleaner CSV version
    csv_df = df.copy()
    # Fill empty count cells for CSV clarity
    current_count_csv = None
    for i in range(len(csv_df)):
        if csv_df.iloc[i]['Count'] != "":
            current_count_csv = csv_df.iloc[i]['Count']
        else:
            csv_df.iloc[i, csv_df.columns.get_loc('Count')] = current_count_csv
    
    csv_df.to_csv(csv_filename, index=False)
    print(f"Publication-ready table saved as {filename} and data saved as {csv_filename}")

# Extract GO terms for nonref and ref gene symbols
go_terms_nonref = extract_go_terms(gene_symbols_nonref)
go_terms_ref = extract_go_terms(gene_symbols_ref)

plot_go_table_grouped(go_terms_nonref, 'GO Terms for Non-reference Genes (≥2 occurrences)', 
                     'go_terms_nonref_table.png', '#E74C3C')  
plot_go_table_grouped(go_terms_ref, 'GO Terms for Reference Genes (≥2 occurrences)', 
                     'go_terms_ref_table.png', '#3498DB')   

'''# Function to plot bar plot for the distribution of all unique values in the 'Name' column
# Function to plot bar plot for the top 5 most frequent 'gene_name' values
def plot_name_distribution_bar(df, title, filename):
    name_counts = df['gene_name'].value_counts()
    top_5 = name_counts.nlargest(10)  # Get the top 5 most frequent names
    labels = top_5.index
    counts = top_5.values
    total = sum(counts)
    proportions = [count / total for count in counts]

    plt.figure(figsize=(30, 25))
    plt.bar(labels, proportions, color="blue")
    plt.xticks(rotation=45, ha="right", fontsize=35)
    plt.yticks(fontsize=35)
    plt.ylabel('Proportion', fontsize=35)
    plt.title(title)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


# Plot and save distribution of the unique values in the 'Name' column
plot_name_distribution_bar(df_nonref, 'Proportion of Name Distribution for Nonref Genes', 'name_distribution_nonref_bar_top10.png')
plot_name_distribution_bar(df_ref, 'Proportion of Name Distribution for Ref Genes', 'name_distribution_ref_bar_top10.png')'''
