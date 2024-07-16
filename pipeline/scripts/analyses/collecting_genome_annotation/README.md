# Collecting genome annotation from NCBI

Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage, dans `/beegfs/banque/gtdrift/data/genome_assembly/XXassemblyXX/genome_seq` et `/beegfs/banque/gtdrift/data/genome_assembly/``XXassemblyXX``/annotation` respectivement. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers bug.

Une commande pour tester ce pipeline est :

``` bash
snakemake collect_everything --configfile /beegfs/banque/gtdrift/pipeline/scripts/analyses/collecting_genome_annotation/config.json  --use-singularity --singularity-args "--bind /beegfs/:/beegfs/" --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --rerun-incomplete --rerun-triggers mtime -j 100 -n --forceall
```

To create a dag file :

``` bash
snakemake collect_everything --configfile /beegfs/banque/gtdrift/pipeline/scripts/analyses/collecting_genome_annotation/config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```
