Generer la liste de GC en json a partir des donnees collectees
export  gc=`grep ">"   ../../../../data_results_per_assembly/genome_assembly/*/annotation/protein.faa |cut -f7 -d"/"|sort -u `
for acc in $gc; do echo "\"$acc\","; done > acc.json

Lancer l'analyse avec la commande:

```
snakemake -s process_stats_prdm9.smk --cores 1
```

