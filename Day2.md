# Genome Assembly
## Quality Control
```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
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

# loop quality control for all of the fastqc files in the path
for file in /work_beegfs/sunam227/metagenomics/0_raw_reads/*.gz; do
    fastqc "${file}" -o /work_beegfs/sunam227/metagenomics/0_raw_reads/quality_reads
done
```


## Read Trimming and Cleaning

```bash
#!/bin/bash
#SBATCH --job-name=log
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --partition=base
#SBATCH --time=5:00:00
#SBATCH --reservation=biol217

fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130305_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130305_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp305_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz -t 6 -q 20
fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130527_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130527_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp527_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz -t 6 -q 20
fastp -i /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130708_mapped_R1.fastq.gz -I /work_beegfs/sunam227/metagenomics/0_raw_reads/BGR_130708_mapped_R2.fastq.gz -R /work_beegfs/sunam227/metagenomics/0_raw_reads/reports/fastp708_report -o /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -O /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz -t 6 -q 20

# this also can be done with loops but this one was more easier so I didn't even try...
```
## Contig Assembling
```bash
#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --output=assembly.out
#SBATCH --error=assembly.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --partition=base
#SBATCH --time=5:00:00
#SBATCH --reservation=biol217

cd /work_beegfs/sunam227/metagenomics/0_raw_reads/clean
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130305_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130527_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12
megahit -1 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R1.fastq.gz -2 /work_beegfs/sunam227/metagenomics/0_raw_reads/clean/BGR_130708_clean_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 /work_beegfs/sunam227/metagenomics/0_raw_reads/assembly/ -t12

# This one will probably take hours to complete, stop worrying about it.

# ##----------------- End -------------
module purge
jobinfo
```