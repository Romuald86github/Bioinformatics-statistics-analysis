import pandas as pd
from skbio.diversity import alpha_diversity
import matplotlib.pyplot as plt
import seaborn as sns

# Load data (species counts and metadata)
species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'

# Load the origin_df (species counts and taxonomy)
origin_df = pd.read_csv(species_counts_path, sep='\t')

# Drop the ASV column as instructed
origin_df = origin_df.drop('ASV', axis=1)

# Separate taxonomy columns from sample columns
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
sample_cols = [col for col in origin_df.columns if col not in taxonomy_cols]

# Extract species counts
species_counts = origin_df[sample_cols].fillna(0)

# Extract the metadata (nests, larvae, locations)
metadata = pd.DataFrame({
    'Larvae': sample_cols,
    'Nest': [larvae[0] for larvae in sample_cols],  # First letter represents the nest
})

# Assign locations based on the nest information
def assign_location(nest):
    if nest in list('ABCDEFGHIJ'):
        return 'Apir'
    elif nest in list('KLMNOPQRST'):
        return 'Ihala'
    elif nest in list('WX'):
        return 'Nsuka'
    # Exclude 'University of Agriculture'
    return 'Unknown'

metadata['Location'] = metadata['Nest'].apply(assign_location)

# (2) Alpha diversity (Shannon) calculation for each larvae
shannon_diversity = alpha_diversity('shannon', species_counts.T.values, ids=species_counts.columns)
shannon_diversity_df = pd.DataFrame({'Larvae': species_counts.columns, 'Shannon': shannon_diversity})
shannon_diversity_df.set_index('Larvae', inplace=True)

# Merge alpha diversity with metadata
alpha_diversity_metadata = shannon_diversity_df.merge(metadata, left_index=True, right_on='Larvae')

# Remove "Unknown" locations (this will exclude only 'University of Agriculture')
alpha_diversity_metadata = alpha_diversity_metadata[alpha_diversity_metadata['Location'] != 'Unknown']

# (3) Plotting Shannon diversity across locations
plt.figure(figsize=(10, 6))
sns.violinplot(x='Location', y='Shannon', data=alpha_diversity_metadata)
plt.title("Shannon Diversity by Location (Excluding University of Agriculture)")

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/violinplot1.png'
plt.savefig(plot_output_path)


plt.show()
