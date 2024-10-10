import pandas as pd
from skbio.diversity import alpha_diversity
from scipy.stats import kruskal

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

# View the resulting dataset
print(alpha_diversity_metadata.head())

# Remove rows with missing Shannon diversity values
alpha_diversity_metadata_clean = alpha_diversity_metadata.dropna(subset=['Shannon'])

# Function to perform Kruskal-Wallis test for each location
def test_nests_within_location(location):
    location_data = alpha_diversity_metadata_clean[alpha_diversity_metadata_clean['Location'] == location]
    if location_data['Nest'].nunique() < 2:
        print(f"Not enough nests to compare in location: {location}")
        return
    groups = [location_data.loc[location_data['Nest'] == nest, 'Shannon'] 
              for nest in location_data['Nest'].unique()]
    H_stat, p_value = kruskal(*groups)
    print(f"Location: {location}, H-statistic = {H_stat}, p-value = {p_value}")

# Perform Kruskal-Wallis test for nests within each location
for location in alpha_diversity_metadata_clean['Location'].unique():
    test_nests_within_location(location)
