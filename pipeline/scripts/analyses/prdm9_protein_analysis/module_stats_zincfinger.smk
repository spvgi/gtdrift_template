import json
import os

configfile: "config.json"
ACCESSNB = config["assembly_list"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

rule generate_prdm9_candidates:
    input:
        prdm9_prot_dir = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/SET_sequences"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates.txt"
    shell:
        """
        ls {input.prdm9_prot_dir} | sed 's/\.[^.]*$//' > {output.candidate_list}
        """

rule run_seqkit_extract:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """

rule zincfinger_analysis:
    """
    Run the zinc finger analysis on each protein sequence using R.
    """
    input:
        protein_seq = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        zincfinger_out = pathGTDriftGlobalResults + "analyses_summaries/zinc_finger/ZFD_{accession}.csv"
    
    run:
        # Verify if the input file exists
        if not os.path.exists(input.protein_seq):
            raise FileNotFoundError(f"Input file does not exist: {input.protein_seq}")
        
        # Run the R script
        command = f"Rscript --vanilla ./zincfinger_analysis.R {input.protein_seq} {output.zincfinger_out}"
        shell(command)

rule combine_zinc_finger:
    """
    Combine all ZFD_{accession}.csv files into one zinc_finger.csv.
    """
    input:
        zfd_files = expand(pathGTDriftGlobalResults + "analyses_summaries/zinc_finger/ZFD_{accession}.csv", accession=ACCESSNB)
    output:
        zinc_finger_combined = pathGTDriftGlobalResults + "analyses_summaries/table_results/zinc_finger.csv"
    run:
        # Check if there are input files before trying to combine
        if not input.zfd_files:
            raise FileNotFoundError("No input files found for combining.")
        
        # Combine all ZFD_{accession}.csv files into one zinc_finger.csv
        with open(output.zinc_finger_combined, 'w') as outfile:
            for fname in input.zfd_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")  # Add a newline between files
