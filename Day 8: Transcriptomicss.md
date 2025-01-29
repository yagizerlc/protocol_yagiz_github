# Transcriptomics
The transcriptome can be defined as the complete set of (RNA) transcripts in a cell, or a population of cells, for a specific developmental stage or physiological condition.

The transcriptome can be summarised as the RNA content of a cell.
- Can include messenger RNA (mRNA), non-coding RNA (ncRNA), small RNA (miRNA, piRNA etc.)
- Can also include transfer RNA (tRNA) and ribosomal RNA (rRNA)
- Nearly all RNA is single-stranded (ssRNA)

We are always dealing with a snapshot of a dynamic process at the time we collected data

## RNA-Seq
Transcriptome sequencing or RNA-seq is a next-generation sequencing based approach to profiling and analysing RNA. This technique delivers unbiased information without the need for prior knowledge of the genome or transcriptome.

Transcriptome sequencing is often the method of choice for analysis of differentially expressed genes, as well as for RNA editing and profiling of allele-specific gene expression. RNA-seq can also be used to investigate splicing patterns, splicing variants, gene isoforms, single nucleotide polymorphisms and post transcriptional modifications.

Transcriptome sequencing begins with isolating RNA and converting this RNA into so-called complementary DNA (cDNA). Total RNA or RNA types like mRNA or small RNA can be investigated. Strand-specific or random-primed libraries can be created based on how the cDNA is synthesised. Once adapters are ligated to the cDNA fragments, the cDNA undergoes single-read or paired-end sequencing. The sequencing reads are either assembled de novo or they are aligned to a reference genome or transcriptome.

There are different types of RNA-Seq techniques:
1. Total RNA Sequencing
2. mRNA Sequencing
3. RNA Exome Sequencing
4. Targeted RNA Sequencing
5. Single Cell Sequencing

### Steps in RNA-Seq Library Preperation
1. RNA Isolation and Purification
2. mRNA Enrichment
3. Fragmentation
4. cDNA Synthesis by Reverse Transcription
5. Adoptor Ligation
6. Amplification of Library
7. Sequencing Through a Sequencing Platform

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

#micromamba activate .micromamba/envs/00_anvio/
```

upload the required modules

```bash
module load gcc12-env/12.1.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"
micromamba activate $WORK/.micromamba/envs/reademption


#set proxy environment to download the data and use the internet in the backend
export http_proxy=http://relay:3128
export https_proxy=http://relay:3128
export ftp_proxy=http://relay:3128

# create folders
reademption create --project_path READemption_analysis --species salmonella="Salmonella Typhimurium"

# Download the files
FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna

#rename the files similar to the genome naming
sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa
wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz

# unzip the file
gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2

#read alignment
reademption align -p 4 --poly_a_clipping --project_path READemption_analysis

# read coverage
reademption coverage -p 4 --project_path READemption_analysis

# gene quantification
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 -c InSPI2,InSPI2,LSP,LSP -r 1,2,1,2 --libs_by_species salmonella=InSPI2_R1,InSPI2_R2,LSP_R1,LSP_R2 --project_path READemption_analysis

# visualzation
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
```








```bash
# ##----------------- End -------------
module purge
jobinfo
```

10 fold change in leuD NC.016810.1
12 fold change in leuB NC_016810.1
15 fold change in psiE



```bash
#load necessary modules
module load gcc12-env/12.1.0
module load micromamba
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba
micromamba activate $WORK/.micromamba/envs/10_grabseqs


reademption align --project_path READemption_analysis \
	--processes 32 --segemehl_accuracy 95 \
	--poly_a_clipping \
	--fastq --min_phred_score 25 \
	--progress

#4- Coverage
reademption coverage --project_path READemption_analysis \
	--processes 32

#5- Performing gene wise quantification
reademption gene_quanti --project_path READemption_analysis \
	--processes 32 --features CDS,tRNA,rRNA 

#6- Performing differential gene expression analysis 

####NOTE:: Change the names according to your file names in the READemption_analysis/input/reads/ directory
reademption deseq --project_path READemption_analysis \
	--libs sRNA154_mutant_R1.fastq.gz,sRNA154_mutant_R2.fastq.gz,wildtype_R1.fastq.gz,wildtype_R2.fastq.gz \
	--conditions sRNA154_mutant_R1,sRNA154_mutant_R2,wildtype_R1,wildtype_R2 --replicates 1,2,1,2 \
	--libs_by_species metanosarcina=sRNA154_mutant_R1,sRNA154_mutant_R2,wildtype_R1,wildtype_R2

#7- Create plots 
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis

#micromamba activate .micromamba/envs/00_anvio/

# ##----------------- End -------------
module purge
jobinfo
```