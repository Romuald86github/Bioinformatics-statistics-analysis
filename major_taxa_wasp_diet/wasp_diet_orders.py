import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

species_counts_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'

# Load the origin_df (species counts and taxonomy)
origin_df = pd.read_csv(species_counts_path, sep='\t')


# (1) Count the number of species per insect order
order_counts = origin_df['order'].value_counts()

# (2) Prepare the bar plot
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(order_counts.index, order_counts.values)

# (3) Customize the plot
ax.set_xlabel('Insect Order')
ax.set_ylabel('Number of Species')
ax.set_title('Number of Species per Insect Order')
plt.xticks(rotation=90)

plt.tight_layout()

plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/order')
plt.show()
