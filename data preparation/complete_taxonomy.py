import pandas as pd

# Load the merged DataFrame (wasp gut results)
merged_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results.tsv'
merged_df = pd.read_csv(merged_file_path, sep='\t')

# Load the sample-origin mapping DataFrame
sample_origin_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/sample_origin_mapping_cleaned.tsv'
sample_origin_df = pd.read_csv(sample_origin_file_path, sep='\t')

# Create a dictionary to map each sample to its origin
sample_to_origin = dict(zip(sample_origin_df['sample'], sample_origin_df['origin']))

# Replace sample columns in merged_df with their corresponding origins
merged_df.rename(columns=sample_to_origin, inplace=True)

merged_df = merged_df.drop('Unnamed: 0', axis =1)

# Save the updated merged DataFrame
updated_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
merged_df.to_csv(updated_file_path, sep='\t', index=False)

print(f"Updated merged table with origins saved as: {updated_file_path}")
