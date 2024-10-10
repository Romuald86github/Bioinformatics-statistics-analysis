import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the alpha diversity results
alpha_diversity_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/alpha_diversity_with_larvae.tsv'
alpha_diversity_df = pd.read_csv(alpha_diversity_path, sep='\t')

# Plot the distribution of alpha diversity indices
plt.figure(figsize=(10, 6))
sns.histplot(alpha_diversity_df['Alpha_Diversity'], bins=20, kde=True)
plt.title('Distribution of Alpha Diversity (Shannon Index) for Wasp Larvae')
plt.xlabel('Alpha Diversity (Shannon Index)')
plt.ylabel('Frequency')
plt.grid(True)
plt.tight_layout()

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/alpha_diversity_distribution.png'
plt.savefig(plot_output_path)

plt.show()

print(f"Alpha diversity distribution plot saved to: {plot_output_path}")
