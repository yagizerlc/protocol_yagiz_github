# Pangenomics

A pangenome is the `complete set of genomic information for a group of organisms`. It encompasses all the natural variation in the DNA of a species or a group of organisms.

There are two main contexts in which the term pangenome is used:
1. **Biological Pangenome:** Refers to the `natural variation in the DNA` of a species or group of organisms, representing the total genetic repertoire of that group.
2. **Computational Pangenome:** Refers to a pangenome that has been `catalogued and constructed` from genomic sequencing data using computational tools. It organizes and represents the genetic variation observed in the biological pangenome.

Pangenomics is the study of the pangenome. It involves:
- Analyzing the `variation within` and `between` the genomes of a group of organisms.
- Constructing `computational pangenomes` using sequencing data.
- Interpreting the biological significance of the genetic variation to understand evolution, function, and adaptation.

In essence, pangenomics focuses on `cataloguing, analyzing, and utilizing` the genomic diversity within a group of organisms.

## Types of Biological Pangenomes


![image](./resources/Screenshot%202025-01-28%20at%2007.30.07.png)

## Types of Computational Pangenomes
### Presence–Absence Variation (PAV) Pangenome
PAV pangenomes quantify the `presence` and `absence` of genes within a population. They identify the `core genome`, which includes all of the genes present in all members of the population, and the `accessory/dispensable genome`, which includes all of the genes present in a subset of the population.

Core gene functions are generally under high selective pressure and are highly conserved within the population. They tend to be older and essential for survival, while accessory genes tend to be less conserved and responsible for variations in lifestyle and evolutionary trajectories.

There are two main strategies for constructing a PAV pangenome:
1. **Homology-based Strategy:** This strategy involves annotating genomes individually, extracting gene sequences, and clustering them based on sequence similarity
2. **map-to-pan Strategy:** Whole-genome sequencing reads are aligned to an annotated representative sequence pangenome to determine gene presence or absence.


### Representative Sequence Pangenome
Representative sequence pangenomes are a collection of genomic sequences that minimize the inclusion of homologous loci while still representing as much genomic diversity from the population as possible.

They are usually composed of a reference genome and a number of other sequences called `nonredundant reference (NRR) sequences`. NRR sequences are sequences that are found in at least one member of the population but are not represented in the reference.

A representative sequence pangenome is constructed by identifying genomic sequences that aren’t already present in the reference genome. These sequences are appended to the reference genome as additional contigs to form a pangenome reference. The pangenome reference may then be optionally annotated.

Methods for Identifying NNR Sequences
1. **Metagenome-like Assembly of Unaligned:** Unaligned reads are treated as metagenomic data and assembled together as a single dataset.
2. **Independent Assembly of Unaligned Reads:** Unaligned reads from each sample are assembled independently and allows for the identification of unique sequences specific to each sample.
3. **Iterative Assembly of Unaligned Reads:** Reads from a single sample that don’t align with the reference are de novo assembled into contigs. These contigs (NRR sequences) are then appended to the reference genome, and this updated reference is used for processing the next sample.
4. **Independent Whole-Genome Assembly:** Reads from each sample are de novo assembled separately into contigs, which are then aligned to the reference genome. Unaligned contigs are pooled and clustered by sequence similarity. The longest sequence (an NRR sequence) is taken from each cluster and appended to the reference genome to form the pangenome.


### Pangenome Graph (Graphical Pangenome)
A sequence-oriented pangenome graph models the `location of genomic variation` within a species with respect to either a reference sequence or to the other sequences comprising the pangenome. They are composed of a set of `nodes` and `edges`. Nodes are segments of genomic sequence, and edges join these segments together.
1. **Sequence-Oriented Graph:** Models sequence variation at the nucleotide level, including its location relative to other sequences.
2. **Gene-Oriented Graph:** Models gene presence/absence and their relative order within the population.

There are three main methods for constructing sequence-oriented pangenome graphs
1. Predetermined variants: Requires a reference genome/sequence and a set of predetermined variants. Using a reference as the base of the graph, each variant is added to the graph as a `bubble`, resulting in a directed acyclic graph ordered along the reference genome.
2. Multiple sequence alignment: Constructed through the alignment of genomic sequences directly against each other.
3. De Bruijn graphs: reads are split into k-mers to form the nodes of the graph, and nodes are joined together based on their k-1 sequence length overlap with one another.

![image](./resources/Pangenometypes.png)

It’s impossible to capture the full extent of variation for larger organisms due to practical constraints in sequencing and analysis.

Pangenomes typically represent only a subset of variation to remain functional.

### Applications
- Estimating Evolutionary Trajectories
- Estimating gene functionality
- Study of functional specialisation and redundancy
- Community Dynamics and Population Genetics
- Genotyping and variant calling
- Haplotype inference
- Functional pangenomics
- Graph pangenome challenges
