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

# Calculate the proportion of larvae with 'presence FAW'
faw_presence_proportion = (df['FAW_Presence'] == 'presence FAW').mean()
print(f"Overall proportion of larvae with 'presence FAW': {faw_presence_proportion:.2%}")

# Visualization (Bar plot)
plt.figure(figsize=(6, 4))
sns.barplot(x=['presence FAW', 'no FAW'], y=[faw_presence_proportion, 1 - faw_presence_proportion])
plt.title('Proportion of Larvae with FAW Presence')
plt.ylabel('Proportion')
plt.show()







# Create a new column for FAW presence
df['FAW_Presence'] = df['ASV_Count'].apply(lambda x: 'presence FAW' if x > 0 else 'no FAW')

# Calculate the proportion of larvae with 'presence FAW'
faw_presence_proportion = (df['FAW_Presence'] == 'presence FAW').mean()
print(f"Overall proportion of larvae with 'presence FAW': {faw_presence_proportion:.2%}")

# Visualization (Pie chart)
plt.figure(figsize=(6, 4))
labels = ['presence FAW', 'no FAW']
sizes = [faw_presence_proportion, 1 - faw_presence_proportion]
plt.pie(sizes, labels=labels, autopct='%1.1f%%')
plt.title('Proportion of Larvae with FAW Presence')
plt.axis('equal')  # Equal aspect ratio ensures that pie is circular.
plt.tight_layout()

# Save the plot
plot_output_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/faw_presence_pie.png'
plt.savefig(plot_output_path)

plt.show()

print(f"Pie chart saved to: {plot_output_path}")