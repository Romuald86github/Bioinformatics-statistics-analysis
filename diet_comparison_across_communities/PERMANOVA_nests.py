import pandas as pd
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova

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

# Calculate Bray-Curtis distance matrix
bray_curtis_distances = beta_diversity("braycurtis", species_counts.T)

# Create metadata DataFrame from sample columns
metadata = pd.DataFrame()

# Extract nest IDs (first character of each sample column)
metadata['Nest'] = [sample[0] for sample in sample_cols]

# Extract larvae IDs (full column names)
metadata['Larvae'] = sample_cols

# Define locations based on nest IDs
def assign_location(nest_id):
    if nest_id in 'ABCDEFGHIJ':  # Nests A-J sampled in Apir
        return 'Apir'
    elif nest_id in 'KLMNOPQRST':  # Nests K-T sampled in Ihala
        return 'Ihala'
    elif nest_id == 'U':  # Nests U sampled at University of Agriculture
        return 'University of Agriculture'
    elif nest_id in 'WX':  # Nests W-X sampled in Nsukka
        return 'Nsukka'
    else:
        return 'Unknown'

metadata['Location'] = metadata['Nest'].apply(assign_location)

# Set 'Larvae' as the index
metadata.set_index('Larvae', inplace=True)

# Remove NaN values from the Bray-Curtis distance matrix
bray_curtis_distances_no_nan = bray_curtis_distances.copy()
bray_curtis_distances_no_nan.data[np.isnan(bray_curtis_distances_no_nan.data)] = 0

# Perform PERMANOVA for beta diversity comparison between nests
# Check if all nest values are unique
if len(metadata['Nest'].unique()) == len(metadata['Nest']):
    print("All nests contain unique groups, PERMANOVA cannot be performed.")
else:
    # Run PERMANOVA to test for significant differences in beta diversity between nests
    permanova_results = permanova(bray_curtis_distances_no_nan, metadata['Nest'], permutations=99999)
    
    print("PERMANOVA results for nests:")
    print(permanova_results)