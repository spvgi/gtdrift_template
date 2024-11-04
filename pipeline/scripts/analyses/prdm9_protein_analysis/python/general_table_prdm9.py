import pandas as pd
import argparse

# Configure the parser for input and output arguments
parser = argparse.ArgumentParser(description='Generates a summary table of domain presence/absence in candidate genes.')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input directory with CSV files')
parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory for the resulting CSV')

args = parser.parse_args()

# Get input and output directories from arguments
input_dir = args.input_dir
output_dir = args.output_dir

# Load the data
dataframes = []
for file in input_dir.split(','):
    file = file.strip()
    df = pd.read_csv(file, sep=';')
    dataframes.append(df)

# Combine all dataframes into one
data = pd.concat(dataframes, ignore_index=True)

# Normalize column names (remove whitespace)
data.columns = data.columns.str.strip()

# Initialize the list for the output DataFrame
summary_data = []

# Group by Species_name and taxid (one row per species)
grouped = data.groupby(['species_name', 'taxid'])

# Iterate over each group
for name, group in grouped:
    kxsz_proteins_list = []
    kxs_proteins_list = []
    other_proteins_list = []
    other_domain_content = []

    # Classify proteins within the group
    for _, row in group.iterrows():
        domains = []
        domain_positions = []

        # Add domain and position if present
        if row['Nb KRAB domains'] > 0:
            domains.append('K')
            domain_positions.append(row['KRAB domain start'])
        if row['Nb SSXRD domains'] > 0:
            domains.append('X')
            domain_positions.append(row['SSXRD domain start'])
        if row['Nb SET domains'] > 0:
            domains.append('S')
            domain_positions.append(row['SET domain start'])
        if row['Nb ZF domains'] > 0:
            domains.append('Z')
            domain_positions.append(row['ZF domain start'])

        # Sort domains by positions to verify order
        sorted_domains = [dom for _, dom in sorted(zip(domain_positions, domains))]

        # Classification based on the domain order
        if sorted_domains == ['K', 'X', 'S', 'Z']:
            kxsz_proteins_list.append(row['SeqID'])
        elif sorted_domains == ['K', 'X', 'S']:
            kxs_proteins_list.append(row['SeqID'])
        else:
            other_proteins_list.append(row['SeqID'])
            other_domain_content.append(''.join(sorted_domains))

    # Add classification results to the summary
    summary_data.append({
        'Species_name': name[0],
        'taxid': name[1],
        'KXSZ nb': len(kxsz_proteins_list),
        'KXSZ list': ', '.join(kxsz_proteins_list),
        'KXS nb': len(kxs_proteins_list),
        'KXS list': ', '.join(kxs_proteins_list),
        'Other nb': len(other_proteins_list),
        'Other list': ', '.join(other_proteins_list),
        'Other domain content': ', '.join(other_domain_content)
    })

# Convert summary_data to DataFrame
summary_df = pd.DataFrame(summary_data)

# Save to CSV
summary_df.to_csv(output_dir, sep=';', index=False)
print(f"Output saved to {output_dir}")
