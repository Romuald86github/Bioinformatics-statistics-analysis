import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Load the wasp gut data with larvae origins and metadata
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
data_df = pd.read_csv(data_path, sep='\t')

# Load metadata containing nests and locations (larvae IDs are set as index)
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'
metadata_df = pd.read_csv(metadata_path, sep='\t', index_col='Larvae')

# Remove the ASV and taxonomy columns (keeping only the larvae abundance data)
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
sample_cols = [col for col in data_df.columns if col not in ['ASV'] + taxonomy_cols]

# Transpose the larvae abundance data so that larvae IDs are rows
abundance_df = data_df[sample_cols].T  # Transpose
abundance_df.index.name = 'Larvae'  # Set index name to match metadata

# Merge the transposed abundance data with the metadata using larvae IDs
merged_df = pd.merge(abundance_df, metadata_df, left_index=True, right_index=True)








# Abundance matrix and metadata (merged_df) are assumed to be ready
abundance_matrix = merged_df.iloc[:, :2180].values  # abundance data without metadata
nests = merged_df['Nest'].values  # nests for each sample
locations = merged_df['Location'].values  # locations for each sample

# Define Hill number function for q=1
def hill_number_q1(abundances):
    abundances = np.array(abundances)
    abundances = abundances[abundances > 0]  # Exclude zeros
    proportions = abundances / abundances.sum()
    return np.exp(-np.sum(proportions * np.log(proportions)))

# Accumulate diversity by location for q=1
def accumulate_diversity_by_location_q1(abundance_matrix, nests, locations):
    results = {}
    unique_locations = np.unique(locations)
    
    for location in unique_locations:
        loc_indices = np.where(locations == location)[0]
        loc_abundances = abundance_matrix[loc_indices, :]
        
        # Initialize diversity accumulation for the location
        hill_accum = []
        for i in range(1, loc_abundances.shape[0] + 1):
            combined_abundance = loc_abundances[:i, :].sum(axis=0)
            diversity = hill_number_q1(combined_abundance)
            hill_accum.append(diversity)
        
        results[location] = hill_accum
    
    return results

# Calculate Hill number for q=1
hill_1 = accumulate_diversity_by_location_q1(abundance_matrix, nests, locations)

# Plot the diversity accumulation curve for q=1
def plot_diversity_accumulation_q1(hill_1, save_dir):
    plt.figure(figsize=(6, 6))
    
    for location, values in hill_1.items():
        plt.plot(range(1, len(values) + 1), values, label=location)
    
    plt.title("Diversity Accumulation (q=1)")
    plt.xlabel("Number of Samples")
    plt.ylabel("Exponent of Shannon Entropy (Hill 1)")
    plt.legend()
    
    # Save the plot
    plt.savefig(f'{save_dir}/diversity_accumulation_q1.png')
    plt.show()

# Define the save directory path
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'

# Plot and save the curve for q=1
plot_diversity_accumulation_q1(hill_1, save_dir)
