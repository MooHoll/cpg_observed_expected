# Scripts for calculating the CpG observed/expected ratio 

This is a pipeline to take any genome fasta and annotation (in gff3 format) and produce a graph of the CpG o/e for genes, exons and introns. 

## Background
DNA methylation (5mC) causes spontaneous deanimation of cytosine nucleotides to thymine nucleotides over time. It is therefore possible to determine if DNA methylation is present in an organism via the number of of CpG sites in a given region relative to the GC content.

A CpG observed/expected frequency of 1.0 is predicted when there is no DNA methylation present in a region, whereas a CpG o/e around 0.5 may suggest the presence of DNA methylation. However, it is not quite this simple as other evolutionary forces may influence the CpG content in a given region. Generally, the presence of two peaks in CpG o/e, one around 0.5 and another around 1.0 suggest DNA methylation presence. However, a left-skewed normal distribution could also signal the presense of DNA methylation.

For more clear information see this publication [Provataris et al. (2018)](https://doi.org/10.1093/gbe/evy066).

## Scripts and usage

**01_make_cpgOE_input_files.sh**
This script takes a genome in fasta format and the corresonding gff3 annotation file. It will produce an output fasta file for genes, introns and exons.

Requires: bedtools.

<code>bash 01_make_cpgOE_input_files.sh species.gff3 species.fa</code>

It produces three output files:
- species_genes.fa
- species_exons.fa
- species_introns.fa

**TO-DO:** I also want to create an extra annotation of 2000bp upstream of genes and also intergenic regions. Consider, creating the CpG o/e dataframe here and then the next script can just be to make the graph. Do all the heavy lifting at once?

Also I need to remove small features of less than 200bp and remove features with a high N content, this will skew the CpG o/e distributions.

----

**02_graph_cpgOE_per_feature.R**
This script takes the above created output files and calculates the CpG o/e ratio for each annotation. 

<code>Rscript 02_graph_cpgOE_per_feature.R species_genes.fa species_exons.fa species_introns.fa </code>

It produces two output files:
- species.pdf (a nice graph with the CpG o/e distribution and peaks for each feature)
- species_data.txt (a file with the CpG o/e data for further analysis, if you want)

**TO-DO:** Tidy up, not very effecient, takes a while to run and will only take the three specificed files. I also now need to fix the peak calling code, it doesn't seem to be picking up the second peak in all species.

----

**03_mixed_distribution_testing.R**
Finally, this script will take the above species_data.txt and force two mixture normal distributions using the R package *mixtools* to check for the presence of two distributions, i.e. methylated features and non-methylated features. 

<code>Rscript 03_mixed_distribution_testing.R species_data.txt </code>

It produces three output files:
- species_genes.pdf
- species_exons.pdf
- species_introns.pdf

**TO-DO:** Find way of adding 95% confidence intervals to the mean of each distribution.
