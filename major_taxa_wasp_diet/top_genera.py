
import pandas as pd
import matplotlib.pyplot as plt

species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'

# Load the origin_df (species counts and taxonomy)
origin_df = pd.read_csv(species_counts_path, sep='\t')


# (1) Count the number of species per genus
genus_counts = origin_df['genus'].value_counts()

# (2) Select the top 10 genera
top_10_genera = genus_counts.nlargest(10)

# (3) Prepare the bar plot
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(top_10_genera.index, top_10_genera.values)

# (4) Customize the plot
ax.set_xlabel('Genus')
ax.set_ylabel('Number of Species')
ax.set_title('Top 10 Genera by Number of Species')
plt.xticks(rotation=90)

plt.tight_layout()

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/top_10_genera.png'
plt.savefig(plot_output_path)

plt.show()

print(f"Bar plot saved to: {plot_output_path}")