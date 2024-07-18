
Importer le depot  gtdrift_template

```
git clone https://github.com/simonpenel/gtdrift_template.git
```

Aller dans le repertoire  gtdrift_template/pipeline/scripts/analyses


```
cd gtdrift_template/pipeline/scripts/analyses
```


Moifier le fichier environment_path.json : remplacer 'my_directory' par le repertoire dans lequel se trouve gtdrift_template



Aller dans le repertoire collecting_genome_annotation et lancer snakemake


```
cd collecting_genome_annotation 
snakemake collect_everything --configfile config.json  --cores 1
```


Aller dans le repertoire frêre prdm9_protein_analysi et lancer snakemake

```
cd ../prdm9_protein_analysis/

snakemake -s process_stats_prdm9.smk --cores 1

```
