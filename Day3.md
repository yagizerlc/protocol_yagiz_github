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

metaquast -t 6 -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/3_metaquast -m 1000 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa

anvi-script-reformat-fasta /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa --min-len 1000 --simplify-names --report-file nameconversion.txt

bowtie2-build /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index

bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam

samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.bam
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.bam
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.bam


# For the contig Data Preperation
anvi-gen-contigs-database -f /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -n 'biol217'



# HMM search 
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/
anvi-run-hmms -c contigs.db




# Sorting and indexing bam files
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
for file in *.bam; do anvi-init-bam $file -o "$file".sorted.bam; done


# creating anvio profile
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-profile -i sample1.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile1
anvi-profile -i sample2.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile2
anvi-profile -i sample3.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile3

# merging profiles
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-merge ./profile1/PROFILE.db ./profile2/PROFILE2.db ./profile3/PROFILE3.db -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof -c contigs.db --enforce-hierarchical-clustering
