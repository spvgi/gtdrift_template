import os

import json

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

include: pathGTDriftScripts + "analyses/collecting_genome_annotation/Snakefile"


#pathBanque = "/beegfs/banque/gtdrift/data/"
#pathScript = "/beegfs/banque/gtdrift/pipeline/"

configfile: "config.json"
#ACCESSNB = [elt for elt in os.listdir('data/ncbi/') if elt.startswith('GC') == True]
ACCESSNB  = config["assembly_list"]
DOMAIN = ['KRAB', 'SET', 'SSXRD', 'ZF']

rule all:
    """
    Get the prdm9 stats
    """
    input:
        stats_prdm9 = expand(pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv", accession=ACCESSNB),
        krab = pathBanque +"analyses_summaries/table_results/krab_data.csv",
        krabzf = pathBanque + "analyses_summaries/table_results/krabzf_data.csv",
        zf = pathBanque + "analyses_summaries/table_results/zf_count.csv",
        table = pathBanque + "analyses_summaries/table_results/table_prdm9.csv"

#include: "module_get_faa.smk"
include: "module_stats_prdm9.smk"
