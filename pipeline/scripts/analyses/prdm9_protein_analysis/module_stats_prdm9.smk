configfile: "config.json"
ACCESSNB  = config["assembly_list"]
DOMAIN = ['KRAB', 'SET', 'SSXRD', 'ZF']
pathBanque = "/beegfs/banque/gtdrift/data/"
pathScript = "/beegfs/banque/gtdrift/pipeline/scripts/"

if config["mode"] == "guix":
    RUNCMD="guix shell hmmer -- "
else:
    RUNCMD=""

rule get_blastdb:
    """
    Genere a blast db
    WARNING: formatdb is used and not makeblastdb. Suffixes and options are different with makeblastdb
    """
    input:
        fasta = pathBanque + "genome_assembly/{accession}/annotation/protein.faa"
    output:    
        psq= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
        psi= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psi",
        psd= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psd",
        pin= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        phr= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr"        
    shell:
        "formatdb -i {input.fasta} -t protdb -n "+pathBanque+"genome_assembly/{wildcards.accession}/analyses/prdm9_prot/protdb -p T -o T"

         
rule hmm_build:
    """
    HMM creation.
    """
    input:
        pathBanque + "files_for_analyses/ref_align/Prdm9_Metazoa_Reference_alignment/Domain_{domain}_ReferenceAlignment.fa"
    output:
       	pathBanque + "files_for_analyses/hmm_build/{domain}.hmm"
    shell:
        "{RUNCMD} hmmbuild {output} {input}"
 
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        model=pathBanque + "files_for_analyses/hmm_build/{domain}.hmm",
        protein=pathBanque + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        table = pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}",
        domains = pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains"
    shell:
        "{RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}"


rule tbl_processing:
    """
    Result file processing for a later use.
    """
    input:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}"
    output:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed"
    shell:
        "python3 "+pathScript+"analyses/prdm9_protein_analysis/python/hmmsearch_parser.py -i {input} -o {output}"

rule domain_processing:
    """
    Result file processing for a later use.
    """
    input:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed",
        domain_data=pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains"
    output:
        processed=pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_processed",
        summary=pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary"
    shell:
       "python3  "+pathScript+"analyses/prdm9_protein_analysis/python/domain_parser.py -i {input.domain_data} -o {output.processed} -s {output.summary}"


def domain_done(wildcards):
    return expand(pathBanque + "genome_assembly/"+ wildcards.accession + "/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary" , domain=DOMAIN)
    #return expand("results/" + wildcards.accession + "/hmm_search/domtbl/{domain}_domains_summary" , domain=DOMAIN)

rule table_editing:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        domain_done,
    output:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv"
    shell:
        "python3 "+pathScript+"analyses/prdm9_protein_analysis/python/table_builder.py -a {wildcards.accession} -o {output}"


rule read_table:
    """
    Reads each summary table and runs a blastp analysis on every candidate
    """
    input:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv",
        psq= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
        psi= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psi",
        psd= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psd",
        pin= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        phr= pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr"   
        
    output:
        pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
    shell:
        "python3 "+pathScript+"analyses/prdm9_protein_analysis/python/blastp_analysis.py "+ pathBanque + "genome_assembly/{wildcards.accession}/analyses/prdm9_prot/summary_table_prdm9_{wildcards.accession}.csv {wildcards.accession}"

rule summary:
    """
    Concatenation of each proteome blastp results.
    """
    input: 
        expand(pathBanque + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt", accession=ACCESSNB)
    output:
        pathBanque + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    shell:
        """
        cat {input} > {output}
        """

rule blastp_results:
    """
    Writing a table from the concatenation
    """
    input:
        pathBanque + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    output:
        pathBanque + "analyses_summaries/BLASTP_results/blastp_results.csv",
    shell:
        "python3 "+pathScript+"analyses/prdm9_protein_analysis//python/blastp_table.py"

rule taxonomy:
    """
    Creation of a table associating a genome accession number to its complete taxonomy
    """
    input:
        pathBanque + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    output:
        pathBanque + "resources/sorted_taxonomy.csv"
    shell:
        "python3 "+pathScript+"analyses/prdm9_protein_analysis/python/taxonomy.py"

rule create_table:
    """
    Creation of multiple result table using blastp results and hmm search results
    """
    input:
        pathBanque + "analyses_summaries/BLASTP_results/blastp_results.csv",
        pathBanque + "resources/sorted_taxonomy.csv"
    output:
        pathBanque + "analyses_summaries/table_results/krab_data.csv",
        pathBanque + "analyses_summaries/table_results/krabzf_data.csv",
        pathBanque + "analyses_summaries/table_results/zf_count.csv",
        pathBanque + "analyses_summaries/table_results/table_prdm9.csv"
    shell:
        """
        python3 /beegfs/banque/gtdrift/pipeline/scripts/analyses/prdm9_protein_analysis/python/krab.py\
        && python3 /beegfs/banque/gtdrift/pipeline/scripts/analyses/prdm9_protein_analysis/python/krabzf.py\
        && python3 /beegfs/banque/gtdrift/pipeline/scripts/analyses/prdm9_protein_analysis/python/zf_analysis.py\
        && python3 /beegfs/banque/gtdrift/pipeline/scripts/analyses/prdm9_protein_analysis/python/table_prdm9.py
        """
