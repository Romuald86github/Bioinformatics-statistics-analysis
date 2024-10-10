import pandas as pd
from skbio.diversity import alpha_diversity
import re

# Load the wasp gut data
wasp_data_path = "/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results.tsv"
wasp_data = pd.read_csv(wasp_data_path, sep="\t")

# (1) Replace each ASV with the associated species name
wasp_data['species'] = wasp_data['species'].fillna(wasp_data['ASV'])
wasp_data.drop('ASV', axis=1, inplace=True)

# (2) Name the new data as species_asv.tsv and save it back to the same path
species_asv_path = "/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/species_asv.tsv"
wasp_data.to_csv(species_asv_path, sep="\t", index=False)

# Load the species ASV data
species_asv_path = "/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/species_asv.tsv"
species_asv = pd.read_csv(species_asv_path, sep="\t")

# Load the sample-origin mapping DataFrame
sample_origin_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/sample_origin_mapping_cleaned.tsv'
sample_origin_df = pd.read_csv(sample_origin_file_path, sep='\t')

# Create a dictionary to map each sample to its origin (larvae ID)
sample_to_origin = dict(zip(sample_origin_df['sample'], sample_origin_df['origin']))

# Rename the sample columns in species_asv to their corresponding larvae IDs
species_asv.rename(columns=sample_to_origin, inplace=True)

# Drop the 'Unnamed: 0' column if it exists
species_asv = species_asv.drop('Unnamed: 0', axis=1, errors='ignore')

# Separate the taxonomy columns from the larvae columns
taxonomy_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
larvae_cols = [col for col in species_asv.columns if col not in taxonomy_cols]

# Convert the larvae data to a NumPy array
asv_matrix = species_asv[larvae_cols].to_numpy()

# Compute alpha diversity for each larvae
alpha_diversity_results = []
for larvae in larvae_cols:
    sample_counts = asv_matrix[:, species_asv.columns.get_loc(larvae)]
    alpha_div = alpha_diversity('shannon', sample_counts)
    alpha_diversity_results.append({'Larvae': larvae, 'Alpha_Diversity': alpha_div})

# Create a DataFrame from the alpha diversity results
alpha_diversity_df = pd.DataFrame(alpha_diversity_results)

# Remove the 'dtype: float64' and whitespaces from the DataFrame before saving
alpha_diversity_df = alpha_diversity_df.applymap(lambda x: re.sub(r'\s+dtype: float64', '', str(x)))


# Remove any leading or trailing whitespace (spaces) from the 'Alpha_Diversity' column
# and replace any leading '0' with an empty string using a regular expression
# Convert the 'Alpha_Diversity' column to a float data type
alpha_diversity_df['Alpha_Diversity'] = alpha_diversity_df['Alpha_Diversity'].str.replace(r'^\s*0\s*', '', regex=True).astype(float)

# Disable the pandas display of data types before saving the output
with pd.option_context('display.float_format', None):
    # Save the alpha diversity results
    alpha_diversity_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/alpha_diversity_with_larvae.tsv'
    alpha_diversity_df.to_csv(alpha_diversity_output_path, sep='\t', index=False)

print(f"Alpha diversity results with larvae IDs saved to: {alpha_diversity_output_path}")