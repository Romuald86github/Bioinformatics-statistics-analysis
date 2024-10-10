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


# (1) Filter the origin_df for rows where species is "Spodoptera frugiperda"
spodoptera_df = origin_df[origin_df['species'].str.contains('Spodoptera frugiperda', case=False, na=False)]

# (2) Extract all the larvae columns (ignoring non-larvae columns like taxonomy)
larvae_columns = origin_df.columns[:-8]  # Assuming last 8 columns are taxonomic columns

# (3) Create a dictionary to store the count of ASVs for each larva
spodoptera_counts_by_larvae = {}

# (4) Loop through each larva and count the non-zero ASVs
for larva in larvae_columns:
    spodoptera_counts_by_larvae[larva] = (spodoptera_df[larva] > 0).sum()

# List of output lines (example)
output_lines = [
    f"Count of Spodoptera frugiperda ASVs in larva {larva}: {count}"
    for larva, count in spodoptera_counts_by_larvae.items()
]


# Loop through each line and extract the larva ID and count
for line in output_lines:
    # Split the line to isolate the larva and count
    larva_info = line.replace("Count of Spodoptera frugiperda ASVs in larva ", "")
    larva, count = larva_info.split(": ")
    spodoptera_counts_by_larvae[larva] = int(count)

# Convert the dictionary to a DataFrame
df = pd.DataFrame(list(spodoptera_counts_by_larvae.items()), columns=['Larva', 'ASV_Count'])

# Print the DataFrame to verify
print(df)







# List of output lines (example)
output_lines = [
    f"Count of Spodoptera frugiperda ASVs in larva {larva}: {count}"
    for larva, count in spodoptera_counts_by_larvae.items()
]

# Initialize an empty dictionary
spodoptera_counts_by_larvae = {}

# Loop through each line and extract the larva ID and count
for line in output_lines:
    # Split the line to isolate the larva and count
    larva_info = line.replace("Count of Spodoptera frugiperda ASVs in larva ", "")
    larva, count = larva_info.split(": ")
    spodoptera_counts_by_larvae[larva] = int(count)

# Convert the dictionary to a DataFrame
df = pd.DataFrame(list(spodoptera_counts_by_larvae.items()), columns=['Larva', 'ASV_Count'])

# Extract the nest and location information from the Larva column
df['Nest'] = df['Larva'].str[0]
df['Location'] = df['Nest'].apply(lambda x: 'Apir' if x in list('ABCDEFGHIJ')
                                 else 'Ihala' if x in list('KLMNOPQRST')
                                 else 'University of Agriculture' if x == 'U'
                                 else 'Nsukka' if x in list('WX')
                                 else 'Unknown')

# Print the updated DataFrame
print(df)



# Create a new column for FAW presence
df['FAW_Presence'] = df['ASV_Count'].apply(lambda x: 'presence FAW' if x > 0 else 'no FAW')


# Proportion of larvae with 'presence FAW' within each nest
nests_proportion = df.groupby('Nest')['FAW_Presence'].apply(lambda x: (x == 'presence FAW').mean())


# Proportion of larvae with 'presence FAW' within each location
locations_proportion = df.groupby('Location')['FAW_Presence'].apply(lambda x: (x == 'presence FAW').mean())


# Set colors for nests S and T (Belonogaster) vs other nests (Ropalidia)
belonogaster_color = 'red'
ropalidia_color = 'blue'

# Create a color palette based on whether the nest is S or T (Belonogaster) or others (Ropalidia)
nest_colors = [belonogaster_color if nest in ['S', 'T'] else ropalidia_color for nest in nests_proportion.index]

# Visualization for nests (bar plot)
plt.figure(figsize=(10, 6))
sns.barplot(x=nests_proportion.index, y=nests_proportion.values, palette=nest_colors)
plt.title("Proportion of Larvae with 'presence FAW' per Nest")
plt.ylabel('Proportion')
plt.xlabel('Nest')

# Add custom legend
handles = [
    plt.Line2D([0], [0], color=belonogaster_color, lw=4, label='Belonogaster (Nests S and T)'),
    plt.Line2D([0], [0], color=ropalidia_color, lw=4, label='Ropalidia (Other Nests)')
]
plt.legend(handles=handles, title="Wasp Genus")

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/faw_presence_per_nest2.png'
plt.savefig(plot_output_path)

plt.show()

print(f"Bar plot for nests saved to: {plot_output_path}")

