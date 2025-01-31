# Genome Assembly
## Quality Control
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

## Quality Control

In order to create a metagenome assembly, we will first check and clean our paired-end read sequences of three metagenomic samples. For quality control, I used `fastqc`, and this provides us fallowing insights into our reads:
1. Per-Base Sequence Quality
2. GC Content
3. Sequence Duplication Level
4. Adapter Contamination
5. Over-Represented Sequences

This step was essential for identifying any potential quality issues in the raw data before proceeding to cleaning. High-quality reads are crucial for accurate metagenome assembly, as metagenomic datasets are inherently more complex due to the presence of multiple species.

```bash
# Quality Control of Raw Read
for file in /work_beegfs/sunam227/metagenomics/0_raw_reads/*.gz; do
    fastqc "${file}" -o /work_beegfs/sunam227/metagenomics/0_raw_reads/quality_reads
done
```


Then I used `fastp` to process the raw reads. This cleaning process involves:
1. Removal of Low-Quality Reads
2. Trimming of Poor Quality Bases from Sequence Ends
3. Removal of Adapter Sequences
4. Filtering of Short Reads for Noise Reduction

This step ensured that only high-quality, clean reads were used for metagenome assembly, which helps reduce errors and biases in the assembly process.

```bash
# Cleaning and Trimming of Raw Reads
# Processing of Raw Sample 1 Reads
fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130305_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130305_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp305_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz -t 6 -q 20

# Processing of Raw Sample 2 Reads
fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130527_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130527_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp527_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz -t 6 -q 20

# Processing of Raw Sample 3 Reads
fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130708_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130708_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp708_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz -t 6 -q 20
```

## Contig Assembling

The processed reads need to be assambled into contigs. Here I used `MEGAHIT`, a tool specifically designed for metagenome assembly. `MEGAHIT` is efficient for handling complex, multi-species datasets and produces high-quality assemblies.

```bash
# Assembly Construction
cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean

# Sample 1
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12

# Sample 2
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12

# Sample 3
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
```

To visualize the assembly graph, I converted the `.fa` file containing our final contigs into a `.fastg` file using the MEGAHIT toolkit's `contig2fastg` tool.

```bash
megahit_toolkit contig2fastg 99 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fa > /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/megahit_out/final.contigs.fastg

# ##----------------- End -------------
module purge
jobinfo
```

then using Bandage, contig graph can be visualised

## Questions
1- Please submit your generated figure and explain in your own words what you can see (keep it short).
![image](./resources/metagenome_graph.png)