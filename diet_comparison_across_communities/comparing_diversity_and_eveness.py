import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# (1) Load data (species counts and metadata)
species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'

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
    elif nest in list('U'):
        return 'University of Agriculture'
    elif nest in list('WX'):
        return 'Nsukka'
    return 'Unknown'

metadata['Location'] = metadata['Nest'].apply(assign_location)

# Save metadata
metadata.to_csv(metadata_path, sep='\t', index=False)

# (2) Alpha diversity (Shannon) calculation for each larvae
from skbio.diversity import alpha_diversity
shannon_diversity = alpha_diversity('shannon', species_counts.T.values, ids=species_counts.columns)
shannon_diversity_df = pd.DataFrame({'Larvae': species_counts.columns, 'Shannon': shannon_diversity})
shannon_diversity_df.set_index('Larvae', inplace=True)

# Merge alpha diversity with metadata
alpha_diversity_metadata = shannon_diversity_df.merge(metadata, left_index=True, right_on='Larvae')

# Set the colors for Belonogaster (S and T) and Ropalidia (all other nests)
belonogaster_color = 'red'
ropalidia_color = 'blue'

# Create a custom color palette as a dictionary based on the nest types
palette = {nest: belonogaster_color if nest in ['S', 'T'] else ropalidia_color for nest in alpha_diversity_metadata['Nest'].unique()}

# (1) Compare and visualize alpha diversity between larvae in each nest overall using violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Nest', y='Shannon', data=alpha_diversity_metadata, palette=palette, inner='quartile')

# Add title and legend
plt.title("Shannon Diversity Between Larvae in Each Nest (Violin Plot)")
handles = [
    plt.Line2D([0], [0], color=belonogaster_color, lw=4, label='Belonogaster (Nests S and T)'),
    plt.Line2D([0], [0], color=ropalidia_color, lw=4, label='Ropalidia (Other Nests)')
]
plt.legend(handles=handles, title="Wasp Genus")

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/violinplot.png'
plt.savefig(plot_output_path)


plt.show()

