# West Nile virus

This GitHub repository aims to sum up an evolutionary analysis performed on the [West Nile virus](https://en.wikipedia.org/wiki/West_Nile_virus).

This work is part of the Biozentrum Research Summer project 2023, carried out at Richard Neher lab.

## Index

The analysis consists in the following steps:

- Nextstrain workflow on west nile virus
- Analysis of the mutations of the tree
    - GTR model
    - Dn/Ds
    - Research of secondary structures

## Nextstrain Workflow

[Nextstrain](https://nextstrain.org/) is an open-source project to harness the scientific and public health potential of pathogen genome data.

We created a workflow using the West nile virus sequences that can be found on [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?VirusLineage_ss=Viruses,%20taxid:10239&SeqType_s=Nucleotide), selecting for the sequences longer than 8000 basepairs. (Before starting the workflow a preliminary analysis on the sequences was performed, i.e. we roughly aligned and created a tree of the sequences to have a first impression about the data.)

We worked on these [sequences](wnv/data/sequences.fasta) and this [metadata](wnv/data/metadata.csv), using this [reference sequence](wnv/config/reference.gb).

The metadata was not totally suitable for the nextstrain workflow, so we processed it a bit (mainly to make the Division information explicit) through the [metadata_processing.py](wnv/metadata_processing.py) script.

The workflow is in the form of a [Snakefile](wnv/Snakefile), the nextstrain functions that we used are:

- index
- filter: sample a certain maximum number of sequences for each country
- align
- tree
- refine: refines the tree using time data. Unfortunately in our case the time data was not very informative, the diversity of wnv is pretty large and we have only data from the last decades, probably wnv evolution started some centuries ago so it is impossible to infer the root of the tree thanks to time. We set the clock-rate to 0.0004 and we rooted by midpoint.
- ancestral
- translate
- traits: we set Country, Division and Host traits
- export: to allow Auspice to correctly place the Divisions on the map we used the [lat_longs.tsv](wnv/config/lat_longs.tsv) file of the SARS-Cov2 Nextstrain workflow.

A detailed description of each funciton can be found [here](https://docs.nextstrain.org/projects/augur/en/stable/).
You can find the output of each tool in the [result](wnv/results/) folder, these are all the files that we will use in our successive analyses and that were used by auspice to build the graphical representation of the workflow. To look at the finished Nextstrain workflow you should install nextstrain on your device and run the "nextstrain view" function on the [auspice](wnv/auspice/) folder.

## Analysis of the mutations

After having created the nextstrain workflow we shifted our focus on the mutations that have occurred on the West Nile virus sequences. Thanks to Augur, we have all of the estimated ancestral sequences, together with their mutations, of the whole tree, they are stored in the [ancestral_node_data.json](wnv/results/ancestral_node_data.json) file.

Just by loading the tree and going though each node, we can see each mutation ever happened in the history of wnv sequences. In order to understand the evolution that happened in west nile virus we want to look at the regions that are free to change, with no selective pression, a rough way to do so, is to consider only the mutations happening in the third position of each codon.

Each mutation is indicated as starting_nucleotide-position-mutated_nucleotide. We computed the mutation rates of each directed mutation (starting_nucleotide-->mutated_nucleotide), taking special care of the first set of mutations from the root: the root of a nextstrain tree corresponds to the sequence of one of the two branches (the longest one) and all of the mutations are considered in one direction from a child to the other, to have a more coherent count, we should count these mutations as occurring in both directions (1/2 for each direction).

### GTR model

A lot of formulas to write

[eigen decomposition explanation](eigen_decomposition.ipynb)

![GTR](images/Figure_3.jpeg)

### Dn/Ds

Deepened analyses were carried out on the mutations, a rigorous way of defining synonymous and non-synonymous mutations was implemented.

With the [KaKs_ratio.py](KaKs_ratio.py) script, we compute the ratio in a sliding window over the genome.

![Dn/Ds](images/Figure_2.jpeg)

### Secondary structures

With the [secondary_structures.py](secondary_structures.py) script, we compute the number of synonymous mutations in a sliding window over the genome.

The regions that have less synonymous mutations are the ones that for some reason, not related to the aminoacid they encode, are strongly conserved. The most common interpretation to these regions is that they are part of important RNA secondary structures, in such case, nucleotides needs to pair and every type of mutation is bad for the hairpin.

![secondary_structures](images/Figure_1.jpeg)