# Collecting genome annotation from NCBI

Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage, dans `/beegfs/banque/gtdrift/data/genome_assembly/XXassemblyXX/genome_seq` et `/beegfs/banque/gtdrift/data/genome_assembly/``XXassemblyXX``/annotation` respectivement. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers bug.

Une commande pour lancer ce pipeline :

``` bash
snakemake collect_everything --configfile config.json  --cores 1
```

Un autre commande avec singularity

``` bash
snakemake collect_everything --configfile config.json  --use-singularity --singularity-args "--bind /beegfs/:/beegfs/" --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --rerun-incomplete --rerun-triggers mtime -j 100 -n --forceall
```

Créer un  dag file :

``` bash
snakemake collect_everything --configfile config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```
