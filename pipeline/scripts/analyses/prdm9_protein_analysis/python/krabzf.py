import pandas as pd
import os
import json


f = open('config.json')
data = json.load(f)
accession = data['assembly_list']
f.close()
# This script takes the previously created krab_data dataframe, that contains the list of every protein presenting a krab domain for every organism,
# and checks for every entries in all of the lists whether these proteins also carry a zinc finger domain or not.

#accession = [elt for elt in os.listdir('results/') if elt.startswith('GC')]
full_data = []

krab_data = pd.read_csv('/beegfs/banque/gtdrift/data/analyses_summaries/table_results/krab_data.csv', sep=';')
for accession_number in accession:
    with open(f"/beegfs/banque/gtdrift/data/genome_assembly/{accession_number}/analyses/prdm9_prot/hmm_search/tbl/ZF_processed") as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name in krab_data.loc[krab_data['Accession'] == accession_number, 'Protein List'].iloc[0]: 
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'KRAB+ZF protein list', 'KRAB+ZF nb'])
zf_data.to_csv('/beegfs/banque/gtdrift/data/analyses_summaries/table_results/krabzf_data.csv', sep= ';', index = False)
