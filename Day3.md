# Assembly Visualisation

Using megahit, metagenomes optimesed for coassembly

```bash
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
```
Then final.contigs.fa file converted to fastg file to visualise in Bandage by megahit toolkit

```bash
megahit_toolkit contig2fastg 99 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fastg
```
then using Bandage, contig graph can be visualised
![image](./resources/graph.png)


with metaquast, quality oc the genome assembly is can be evaluated

```bash
metaquast -t 6 -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/3_metaquast -m 1000 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa
```
# Genome Binning

firstly fasta sequence IDs needed to be reformatted for binning step to work  properly.

```bash
anvi-script-reformat-fasta /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa --min-len 1000 --simplify-names --report-file nameconversion.txt
```
Then using bowtie2, the clean reads needed to be mapped onto assembled contigs

```bash
bowtie2-build /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index

# for sample 1
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam

# for sample 2
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam

# for sample 3
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam
```

the output of the bowtie2 is in .sam (Sequence Mapping File). After obtaining our .sam file, we will convert it to .bam file for further processes by using samstool.

```bash
# sample 1
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.bam

# sample 2
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.bam

# sample 3
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.bam
```

# Contig Data Preperation

For contig data preperation steps we require anvi´o contigs-db database. This contigs database contains information about sequences in our samples.

```bash
anvi-gen-contigs-database -f /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -n 'biol217'
```

then Hidden Markov Model search should be run on contigs to search for specific genes.

```bash
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/
anvi-run-hmms -c contigs.db
```

# Binning with Anvi´O

```bash
# Sorting and Indexing bam files
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
for file in *.bam; do anvi-init-bam $file -o "$file".sorted.bam; done
```

For binning we need to create anvi´o profile to store sample specific information about contigs

```bash
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean

# Sample 1
anvi-profile -i sample1.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile1

# Sample 2
anvi-profile -i sample2.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile2

# Sample 3
anvi-profile -i sample3.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile3

# Then the created profiles for each sample we had needed to be merged into single profile
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-merge ./profile1/PROFILE.db ./profile2/PROFILE2.db ./profile3/PROFILE3.db -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof -c contigs.db --enforce-hierarchical-clustering
```

Now genome binning can be performed by using both Metabat2 or MaxBin2. In genome binning, contigs belonging to same genome are grouping together.

```bash
# Binning with Metabat2
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-cluster-contigs -p ./mergedprof/PROFILE.db -c contigs.db -C METABAT2 --driver metabat2 --just-do-it --log-file log-metabat2

anvi-summarize -p ./mergedprof/PROFILE.db -c contigs.db -o SUMMARY_METABAT2 -C METABAT2

# binning with MaxBin2
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-cluster-contigs -p ./mergedprof/PROFILE.db -c contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2

anvi-summarize -p ./mergedprof/PROFILE.db -c contigs.db -o SUMMARY_MAXBIN2 -C MAXBIN2
```

In the summary we can see that number of Archeal bin we got is 3 for Metabat2 and 0 for Maxbin2.

After binning we need to estimate our genome completnes and contamination levels

```bash
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-estimate-genome-completeness -c contigs.db -p ./mergedprof/PROFILE.db -C METABAT2
anvi-estimate-genome-completeness -p ./mergedprof/PROFILE.db -c contigs.db --list-collections
```
