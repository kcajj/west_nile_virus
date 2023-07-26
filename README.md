# West Nile virus

This GitHub repository aims to sum up an evolutionary analysis performed on the West Nile virus.

This work is part of the Biozentrum Research Summer project 2023, at Richard Neher lab.

## Index

The analysis consists in the following steps:

- Nextstrain workflow on west nile virus
- Analysis of the mutations on the tree
- GTR model
- Research of secondary structures

## Nextstrain

[Nextstrain](https://nextstrain.org/) is an open-source project to harness the scientific and public health potential of pathogen genome data.

We created a workflow using the West nile virus sequences that can be found on [NCBI_virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?VirusLineage_ss=Viruses,%20taxid:10239&SeqType_s=Nucleotide), selecting for the sequences longer than 8000 basepairs.

(before starting the workflow a preliminary analysis on the sequences was performed, i.e. we roughly aligned and created a tree of the sequences to have a first impression about the data)

We worked on these [sequences](wnv/data/sequences.fasta) and this [metadata](wnv/data/metadata.csv); using [this](wnv/config/reference.gb) as reference sequence.

The metadata was not totally suitable for the nextstrain workflow, so we processed it a bit (mainly to make the Division information explicit) through a [python script](wnv/metadata_processing.py).

The workflow is in the form of a [Snakefile](wnv/Snakefile), the nextstrain functions that we used are:

- index
- filter: sample a certain maximum number of sequences for each country
- align
- tree
- refine: refines the tree using time data, we set the clock-rate to 0.0004
- ancestral
- translate
- traits: we set Country, Division and Host traits
- export: to allow Auspice to correctly place the Divisions on the map we used the [lat_longs.tsv](wnv/config/lat_longs.tsv) file of the SARS-Cov2 Nextstrain workflow.

A detailed description of each funciton can be found [here](https://docs.nextstrain.org/projects/augur/en/stable/).
You can find the output of each tool in the [result](wnv/results/) folder, these are all the files that we will use in our successive analyses and that were used by auspice to build the graphical representation of the workflow. To look at the finished Nextstrain workflow you should install nextstrain on your device and run the "nextstrain view" function on the [auspice](wnv/auspice/) folder.

