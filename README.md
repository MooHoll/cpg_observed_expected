# Scripts for calculating the CpG observed/expected ratio 

This is a pipeline to take any genome fasta and annotation (in gff3 format) and produce a graph of the CpG o/e for genes, exons and introns. 

## Background
DNA methylation (5mC) causes spontaneous deanimation of cytosine nucleotides to thymine nucleotides over time. It is therefore possible to determine if DNA methylation is present in an organism via the number of of CpG sites in a given region relative to the GC content.

A CpG observed/expected frequency of 1.0 is predicted when there is no DNA methylation present in a ragion, whereas a CpG o/e < 1.0 may suggest the presence of DNA methylation. However, it is not quite this simple as other evolutionary forces may influence the CpG content in a given region. Generally, the presence of two peaks in CpG o/e, one around 0.5 and another around 1.0 suggest DNA methylation presence.

For more clear information see this publication [Provataris et al. (2018)](https://doi.org/10.1093/gbe/evy066).

## Scripts and usage

**01_make_cpgOE_input_files.sh**
This script takes a genome in fasta format and the corresonding gff3 annotation file. It will produce an output fasta file for genes, introns and exons.

Requires: bedtools.

<code>bash 01_make_cpgOE_input_files.sh your_genome.gff3 your_genome.fa</code>

It produces three output files:
- your_genome_genes.fa
- your_genomes_intron.fa
- your_genome_exons.fa

**TO-DO:** I also want to create an extra annotation of 2000bp upstream of genes and also intergenic regions.

----

**02_XXX.R**
This script takes the above created output files and calculates the CpG o/e ratio for each annotation. There are currently some manual cut-offs which must be defined within the script by eye-balling the CpG o/e distribution in your species. 

The output of this script is a nice plot which shows the CpG o/e distribution for each annotation given.

**TO-DO:** Speed up the length calculation and find a way to remove any manual input for cut-offs.

----

**03_XXX.R**
Finally, this script will take your feature of choice (currently just one at a time), e.g. just genes, and force two mixture normal distributions using the R package *mixtools* to check for the presence of two distributions, i.e. methylated features and non-methylated features. 

**TO-DO:** Find way of adding 95% confidence intervals to the mean of each distribution.