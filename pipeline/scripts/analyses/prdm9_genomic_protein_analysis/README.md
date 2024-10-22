Récupere l'indfo sur les assemblages avec la commande:

```
snakemake -s fetch_data.smk --cores 1
```

Generer le fichier de conf pour telecharger les assemblages

```
python3 generate_conf_for_collecting_genome_annotation.py data/resources/organisms_data config_prdm9_genomic.json
```


Lancer l'analyse avec la commande:

```
snakemake -s process_stats_prdm9.smk --cores 1
```

