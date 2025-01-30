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
module load micromamba
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba

cd $WORK

micromamba activate .micromamba/envs/00_anvio/
```

Previous day we prepared our metagenome assembly and now we will continue from then. We will start with checking quality of our assembly.

## Assembly Quality Estimation
I used the `grep` command to count the number of contigs in the `final.contigs.fa` file. This step provides an overview of the assembly output by determining the total number of contigs generated, which can indicate the fragmentation level of the assembly.

```bash
grep -c ">" final.contigs.fa
```

Then I used `MetaQUAST` to assess the quality of the metagenome assembly. With `MetaQUAST`, quality of the genome assembly is can be evaluated

```bash
metaquast -t 6 -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/3_metaquast -m 1000 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa
```

`MetaQUAST` evaluates various metrics like `N50 Value`, that is a measure of assembly contiguity. A higher N50 value indicates a more contiguous assembly, which is particularly important in metagenomics as it reflects better resolution of genomic regions and fewer fragmented contigs.

### Questions
1. What is your N50 value? Why is this value relevant?
2. How many contigs are assembled?
- 54812
3. What is the total length of the contigs?
- 140329960




## Genome Binning

Genome binning is the process of grouping contigs from metagenome assemblies into bins, each representing a draft genome of an individual species or strain in the sample.

In metagenomics, the assembly typically contains contigs from multiple species in the community. Binning enables us to reconstruct individual genomes, which is critical for:
- Identifying and characterizing species present in the sample.
- Studying the functional and metabolic potential of each organism.
- Understanding microbial interactions in the environment.

To ensure compatibility with Anvi’o, I reformatted the FASTA sequence IDs using `anvi-script-reformat-fasta` command. This step is required for downstream processing in Anvi’o.

```bash
# Preparing the Fasta File for Anvi'o
anvi-script-reformat-fasta /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa --min-len 1000 --simplify-names --report-file nameconversion.txt
```

### Mapping

I used `bowtie2-build` to create an index of the reference fasta file. Indexing allows `Bowtie2` to efficiently map reads to the reference sequences.

```bash
# Indexing the Referance Genome
bowtie2-build /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index
```

Next, I used `Bowtie2` to map the processed reads back to the assembled contigs, generating `.sam` (Sequence Alignment/Map) files. A SAM file is a text-based format that stores read alignments to the reference genome. It contains detailed information about each mapped read, such as:
- The read sequence
- The reference sequence it aligns to
- The position of the alignment
- Alignment quality scores

```bash
# Mapping with Bowtie2
# for sample 1
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam

# for sample 2
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam

# for sample 3
bowtie2 --very-fast -x /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs-anvio.fa.index -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz -S /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam
```

Then I converted the `.sam` files into `.bam` files using `samtools`. A BAM file is the binary version of a SAM file.
- Smaller in size, saving storage space.
- Faster to process due to its binary format.
- Enables downstream analysis with tools that require BAM files.

```bash
# Conversing SAM Files to BAM Files
# sample 1
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample1.bam

# sample 2
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample2.bam

# sample 3
samtools view -bS /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.sam > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/sample3.bam
```

For contig data preperation steps we require anvi´o `contigs-db` database. This contigs database contains information about sequences in our samples.

```bash
anvi-gen-contigs-database -f /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/contigs.anvio.fa -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/contigs.db -n 'biol217'
```

After, I ran anvi-run-hmms to perform a Hidden Markov Model (HMM) search on the contigs to search for specific genes.

```bash
# HMM Search
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/
anvi-run-hmms -c contigs.db
```

### Binning with Anvi´O

Then I performed genome binning using Anvi’o, which groups contigs into bins based on:
- K-mer frequencies
- Coverage data from the mapped reads
- Taxonomic signals from marker genes

In order to do that first, we will sort and index our BAM files.

```bash
# Sorting and Indexing BAM files
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
for file in *.bam; do anvi-init-bam $file -o "$file".sorted.bam; done
```

After indexing our BAM files, I generated Anvi’o profiles for each sample using `anvi-profile`, which stores sample-specific information for each contig. These profiles then were merged into a single profile using `anvi-merge`.

```bash
# Anvi' Profile Generation
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean

# Sample 1
anvi-profile -i sample1.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile1

# Sample 2
anvi-profile -i sample2.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile2

# Sample 3
anvi-profile -i sample3.bam.sorted.bam -c contigs.db --output-dir  /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/profile3
```

```bash
# Merging Anvi'o Profiles
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-merge ./profile1/PROFILE.db ./profile2/PROFILE2.db ./profile3/PROFILE3.db -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/mergedprof -c contigs.db --enforce-hierarchical-clustering
```

To complement Anvi’o binning, I used `MetaBAT2` and `MaxBin2`, two widely used automated binning tools. Each tool uses different algorithms. `MetaBAT2` groups contigs based on coverage and sequence composition across multiple samples, while `MaxBin2` uses probabilistic methods to assign contigs to bins.

```bash
# Binning with Metabat2
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-cluster-contigs -p ./mergedprof/PROFILE.db -c contigs.db -C METABAT2 --driver metabat2 --just-do-it --log-file log-metabat2

anvi-summarize -p ./mergedprof/PROFILE.db -c contigs.db -o SUMMARY_METABAT2 -C METABAT2
```

### !!!PUT METABAT RESULT!!!

```bash
# binning with MaxBin2
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-cluster-contigs -p ./mergedprof/PROFILE.db -c contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2

anvi-summarize -p ./mergedprof/PROFILE.db -c contigs.db -o SUMMARY_MAXBIN2 -C MAXBIN2
```


## Questions
1- Number of ins you got from MetaBAT2 and MaxBin2?
- In the summary we can see that number of Archeal bin we got is 3 for Metabat2 and 0 for Maxbin2.

## MAGs Quality Estimation

MAGs (Metagenome-Assembled Genomes) are genomes reconstructed from metagenomic data. Instead of sequencing individual organisms in a microbial community, metagenomics captures all the DNA in a sample. MAGs are created by assembling and binning this mixed genomic data, where each MAG ideally represents the genome of a single species or strain.

Importance of MAGs
- Enable the reconstruction of genomes for uncultured or hard-to-culture microbes directly from environment.
- Allow us to study the metabolic potential, ecological roles, and interactions of individual organisms within the microbial community.
- Provide insights into microbial diversity, evolution, and functional pathways in complex environments.

While MAGs are powerful, they can sometimes be incomplete or contaminated (contain sequences from multiple organisms), which is why assessing their quality is crucial.

Here I used `anvi-estimate-genome-completeness` to evaluate the quality of the MAGs. This step estimates:
- **Genome Completeness**: The proportion of conserved single-copy marker genes present in the genome. Higher completeness indicates a more complete representation of the genome.
- **Contamination**: The presence of duplicate or non-target marker genes, which suggests that the bin may contain sequences from multiple organisms. Minimizing contamination ensures the MAG represents a single organism.

This step is critical because accurate genome reconstructions rely on maximizing completeness while minimizing contamination.

```bash
# Quality Assesment of MAGs
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
anvi-estimate-genome-completeness -c contigs.db -p ./mergedprof/PROFILE.db -C METABAT2
anvi-estimate-genome-completeness -p ./mergedprof/PROFILE.db -c contigs.db --list-collections
```

The quality of each bin was assessed based on:
- **High-Quality MAGs**: >90% completeness and <5% contamination.
- **Medium-Quality MAGs**: >50% completeness and <10% contamination.

This step ensures that the MAGs are of sufficient quality to proceed with downstream analyses, such as functional annotation, metabolic reconstruction, and phylogenetic classification.

## Questions
1- Which binning strategy gives you the best quality for the Archeal bins?
- MetaBAT2

2- How many Archeal bins do you get that are of High Quality? How many Bacterial bins do you get that are of High Quality?
- 0 Archeal bins in high-quality, 13 Bacterial bins with high-quality