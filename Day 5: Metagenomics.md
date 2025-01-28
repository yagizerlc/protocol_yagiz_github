```bash
module load gcc12-env/12.1.0
module load micromamba/1.3.1
micromamba activate 00_anvio
```

Now we will add taxonomic annotations to our MAGs. for this we will first upload taxonomy information to our contig database.

```bash
anvi-run-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -T 20 -P 2

anvi-estimate-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -p /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy

anvi-estimate-scg-taxonomy -c /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -p /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt
```

then we will finalise our MAGs with taxonomical information

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
- Yes. only one of the archeal bin is assigned as high qualiy. The other tqo are in low quality.