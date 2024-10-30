import pandas as pd
import argparse
import os
from ete3 import NCBITaxa

'''
This script takes the "summary_table_(Genome_Accesion).csv" file as input and filters only
the PRDM9 candidates whose "Best Match" is PRDM9. The output is a .csv file with the input information 
for these filtered sequences, along with new columns for "species_name", "taxid", and "Assembly".
'''

def filter_prdm9(input_file, output_file):
    # Initialize NCBITaxa to retrieve the species name
    ncbi = NCBITaxa()
    
    # Extract the Genome_assembly code from three levels up in the directory path
    genome_assembly = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(input_file))))
    
    try:
        # Read the original CSV file using a semicolon as separator
        df = pd.read_csv(input_file, sep=';')
    except pd.errors.EmptyDataError:
        # If the input file is empty, create an empty output file
        pd.DataFrame(columns=['species_name', 'taxid', 'Assembly']).to_csv(output_file, sep=';', index=False)
        return

    # Check if the DataFrame is empty after reading
    if df.empty:
        # Create an empty CSV file if the input has no data
        pd.DataFrame(columns=['species_name', 'taxid', 'Assembly']).to_csv(output_file, sep=';', index=False)
    else:
        # Filter rows where the 'Best Match' column is equal to 'PRDM9'
        df_filtered = df[df['Best Match'] == 'PRDM9']
        
        # If the filtered DataFrame has no rows, create an empty output file
        if df_filtered.empty:
            pd.DataFrame(columns=['species_name', 'taxid', 'Assembly']).to_csv(output_file, sep=';', index=False)
            return
        
        # Get the taxid from the first row of the filtered rows
        taxid = int(df_filtered['Taxid'].iloc[0])
        
        # Retrieve the species name using the taxid
        species_name = ncbi.get_taxid_translator([taxid]).get(taxid, 'Unknown species')
        
        # Add new columns for "species_name", "taxid", and "Assembly"
        df_filtered.insert(0, 'species_name', species_name)
        df_filtered.insert(1, 'taxid', taxid)
        df_filtered.insert(2, 'Assembly', genome_assembly)
        
        # Save the result to a new CSV file
        df_filtered.to_csv(output_file, sep=';', index=False)

if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Filter sequences where the Best Match is PRDM9 and add species information.")
    
    # Add arguments for the input file and output file
    parser.add_argument('input_file', type=str, help="Path to the input CSV file")
    parser.add_argument('output_file', type=str, help="Path to the output CSV file (filtered result)")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to perform the filtering
    filter_prdm9(args.input_file, args.output_file)

