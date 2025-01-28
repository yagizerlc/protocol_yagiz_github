```bash
#!/bin/bash
#SBATCH --job-name=TASK
#SBATCH --output=TASK.out
#SBATCH --error=TASK.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --partition=base
#SBATCH --time=5:00:00
#SBATCH --reservation=biol217

#load necessary modules
module load gcc12-env/12.1.0
module load gcc12-env/12.1.0
module load micromamba/1.3.1

cd $WORK

micromamba activate .micromamba/envs/00_anvio/
```


Prewius day we performed bin refinment and chimera detection by using `GUNC`. And today, we will continue with assigning taxnomy to our metagenome assembly. This part will be the last steps of our metagenome work and then we will do some genome assembly next day. 

**Now let's start.**

## Taxonomic Assignment

Here we will add taxonomic annotations to out MAGs.Taxonomic assignment of MAGs is crucial to:
1. Identifying the organisms in the sample: It allows us to determine the taxonomic identity of each reconstructed genome, providing insights into the microbial diversity in the metagenome.
2. Linking genome function to taxonomy: Knowing the taxonomy enables the association of metabolic and functional potential with specific taxa.
3. Enabling ecological and evolutionary studies: Taxonomic classification helps study microbial interactions, niche specialization, and evolutionary relationships.

For this, first we need to upload taxonomy information to our contig.db database, because the tool we are going to use associates the single-copy core genes (SCGs) identified in the contigs database with known taxonomy information.

```bash
# Uploading Taxonomy Information to contig.db
anvi-run-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -T 20 -P 2
```

I used the `anvi-run-scg-taxonomy` command to assign taxonomy to our MAGs.

```bash
# Estimate Coverage of Dataset (Abundance of rRNA)
anvi-estimate-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -p /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy

# Saving the Output
anvi-estimate-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -p /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt
```

then we will finalise our MAGs with taxonomical information.

```bash
anvi-summarize -p /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof/PROFILE.db -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/SUMMARY_METABAT2_FINAL -C METABAT2
```

## Questions
1- Did you get a species assignment to archea bins previously identified?
- Yes, previosly assigned archeal bins assigned with fallowing species:
    
    1- METABAT__10 - *Methanoculleus sp012797575*
    
    2- METABAT__36 - *Methanoculleus thermohydrogenotrophicum*

    3- METABAT__8 - *Methanosarcina flavescens*


2- Does the HIGH-QUALITY assignment of the bin needed revision?
- Yes. none of the Archeal bins got HIGH-QUALITY assignment.