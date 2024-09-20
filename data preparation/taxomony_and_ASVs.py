import pandas as pd

# Define file paths
asv_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/table_ASVs.tsv'
taxonomy_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/TAXONOMY_processed.97.tsv'
output_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results.tsv'


# Load the ASV count table
asv_df = pd.read_csv(asv_file_path, sep='\t')

# Rename the first column to 'ASV'
asv_df.rename(columns={asv_df.columns[0]: 'ASV'}, inplace=True)

# Load the taxonomy table, skipping the first row (original header)
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t', header=None, skiprows=1, engine='python')

# Create meaningful column names in reversed order including subspecies
taxonomy_df.columns = [
    'ASV', 'nt_database_id', 'similarity_percentage', 'alignment_length',
    'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end',
    'e_value', 'bit_score', 'sequence_length', 'kingdom', 'phylum', 'class', 
    'order', 'family', 'genus', 'species', 'subspecies'
]

# Reduce to only the relevant columns and reorder them
taxonomy_df = taxonomy_df[['ASV', 'similarity_percentage', 'nt_database_id', 
                           'kingdom', 'phylum', 'class', 'order', 
                           'family', 'genus', 'species', 'subspecies']]

# Filter asv_df to include only ASVs present in the taxonomy_df
asv_df_filtered = asv_df[asv_df['ASV'].isin(taxonomy_df['ASV'])]

# Merge the filtered ASV count table with the taxonomy table on the 'ASV' column
merged_df = pd.merge(asv_df_filtered, taxonomy_df, on='ASV', how='left')

# Drop the 'similarity_percentage', 'nt_database_id' and duplicate samples columns to create merged_df
merged_df = merged_df.drop(columns=['S177', 'S94', 'S29', 'S93', 'S39', 'S82', 'similarity_percentage', 'nt_database_id'])


# Save the resulting dataframe to a new TSV file
merged_df.to_csv(output_file_path, sep='\t', index=True)

print(f"Merged table without 'similarity_percentage' and 'nt_database_id' saved as: {output_file_path}")
