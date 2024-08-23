import pandas as pd
import os

# This script is used to count the number of proteins with 5 or more Zinc finger domains for every organism

accession = [elt for elt in os.listdir('results/') if elt.startswith('GC')]
full_data = []

for accession_number in accession:
    with open(f"results/{accession_number}/Step4_Hmm/domtbl/ZF_domains_summary") as reader:
        lines = reader.readlines()[1:]
        i = 0
        prot_count = 0
        while i < len(lines):
            domains_nb = int(lines[i].split('\t')[10])
            if domains_nb >= 5:
                prot_count += 1
            i += domains_nb
        full_data.append([accession_number, prot_count])

zf_data = pd.DataFrame(full_data, columns=['Accession', '5+ ZF'])
zf_data.to_csv('summarized_results/zf_count.csv', sep = ';', index= False)
