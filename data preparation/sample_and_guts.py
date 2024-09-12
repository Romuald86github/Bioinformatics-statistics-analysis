import os
import pandas as pd
import re

# Define the path to the cutadapt folder
cutadapt_folder_path = '/Users/romualdchristialtcheutchoua/Downloads/IK_data/cutadapt'

# Initialize lists to store samples and their origins
samples = []
origins = []

# Regular expression pattern to extract sample (e.g., S1) and origin (e.g., A5)
pattern = re.compile(r'([A-Z]\d+)_S(\d+)_L001_R1_001.fastq.gz')

# Iterate over each file in the cutadapt folder
for filename in os.listdir(cutadapt_folder_path):
    # Search for the pattern in the filename
    match = pattern.search(filename)
    if match:
        origin = match.group(1)  # Extract the specific origin (e.g., A5)
        sample = f"S{match.group(2)}"  # Extract sample (e.g., S1)

        samples.append(sample)
        origins.append(origin)

# Create a DataFrame with the extracted data
sample_origin_df = pd.DataFrame({
    'sample': samples,
    'origin': origins
})

# Sort the DataFrame by sample name for clarity
sample_origin_df.sort_values('sample', inplace=True)

# Print duplicate values in the 'sample' column along with their corresponding 'origin' values
print("Duplicate samples and their origins:")
duplicates = sample_origin_df[sample_origin_df.duplicated('origin', keep=False)]
for index, row in duplicates.iterrows():
    print(f"Sample: {row['sample']}, Origin: {row['origin']}")

# Drop duplicate rows based on the 'sample' column
sample_origin_df.drop_duplicates('origin', keep='first', inplace=True)

# Display the resulting DataFrame
print("\nCleaned DataFrame:")
print(sample_origin_df)

# Save the DataFrame to a CSV or TSV file
output_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/sample_origin_mapping_cleaned.tsv'
sample_origin_df.to_csv(output_file_path, sep='\t', index=False)

print(f"\nCleaned sample-origin mapping saved as: {output_file_path}")