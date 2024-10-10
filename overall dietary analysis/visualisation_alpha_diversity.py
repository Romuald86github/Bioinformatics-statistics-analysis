import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal

# (1) Load data (species counts and metadata)
species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'

# Load the origin_df (species counts and taxonomy)
origin_df = pd.read_csv(species_counts_path, sep='\t')

# Drop the ASV column
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
    elif nest in list('U'):
        return 'University of Agriculture'
    elif nest in list('WX'):
        return 'Nsukka'
    return 'Unknown'

metadata['Location'] = metadata['Nest'].apply(assign_location)

# Save metadata
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'
metadata.to_csv(metadata_path, sep='\t', index=False)

# (2) Alpha diversity (Shannon) calculation for each larvae
shannon_diversity = alpha_diversity('shannon', species_counts.T.values, ids=species_counts.columns)
shannon_diversity_df = pd.DataFrame({'Larvae': species_counts.columns, 'Shannon': shannon_diversity})
shannon_diversity_df.set_index('Larvae', inplace=True)

# Merge alpha diversity with metadata
alpha_diversity_metadata = shannon_diversity_df.merge(metadata, left_index=True, right_on='Larvae')

# (4) Visualize alpha diversity (overall larvae)
plt.figure(figsize=(12, 6))
sns.scatterplot(x='Larvae', y='Shannon', data=alpha_diversity_metadata)
plt.xticks(rotation=90)
plt.title("Shannon Diversity for All Larvae")
plt.xlabel("Larvae")
plt.ylabel("Shannon Diversity")
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/alpha_diversity_all_larvae.png')
plt.show()