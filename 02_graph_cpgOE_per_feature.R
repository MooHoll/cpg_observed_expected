#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(echo=TRUE)

# Takes the genes, exons, introns fasta files as command line arguments
#---------------------------------------------
# Identify trends in CpG o/e for different genomic Features
#---------------------------------------------

library(readr)
library(seqinr)
library(stringr)
library(sqldf)
library(ggplot2)
library(ggpmisc)

species <- tools::file_path_sans_ext(basename(args[1]))
species <- gsub("_genes","",species)

#---------------------------------------------
# Count number of C's, G's and CpGs per Feature and get the seq length
#---------------------------------------------
gene_sequence <- read.fasta(args[1])

seq_names <- vector("character",length = length(gene_sequence))
c_count <- vector("numeric",length = length(gene_sequence))
g_count <- vector("numeric",length = length(gene_sequence))
cg_count <- vector("numeric",length = length(gene_sequence))

for(i in seq_along(gene_sequence)){
  seq_names[i] <- names(gene_sequence)[i]
  count <- as.data.frame(count(gene_sequence[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
  count2 <- as.data.frame(count(gene_sequence[[i]],2))
  cg_count[i] <- count2[7,"Freq"]
}
lengths <- getLength(gene_sequence)

gene_data <- data.frame(seq_names,c_count,g_count, cg_count, lengths)
gene_data$Feature <- "Genes"

#---------------------------------------------
gene_sequence <- read.fasta(args[2])

seq_names <- vector("character",length = length(gene_sequence))
c_count <- vector("numeric",length = length(gene_sequence))
g_count <- vector("numeric",length = length(gene_sequence))
cg_count <- vector("numeric",length = length(gene_sequence))

for(i in seq_along(gene_sequence)){
  seq_names[i] <- names(gene_sequence)[i]
  count <- as.data.frame(count(gene_sequence[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
  count2 <- as.data.frame(count(gene_sequence[[i]],2))
  cg_count[i] <- count2[7,"Freq"]
}
lengths <- getLength(gene_sequence)

exon_data <- data.frame(seq_names,c_count,g_count, cg_count, lengths)
exon_data$Feature <- "Exons"

#---------------------------------------------
gene_sequence <- read.fasta(args[3])

seq_names <- vector("character",length = length(gene_sequence))
c_count <- vector("numeric",length = length(gene_sequence))
g_count <- vector("numeric",length = length(gene_sequence))
cg_count <- vector("numeric",length = length(gene_sequence))

for(i in seq_along(gene_sequence)){
  seq_names[i] <- names(gene_sequence)[i]
  count <- as.data.frame(count(gene_sequence[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
  count2 <- as.data.frame(count(gene_sequence[[i]],2))
  cg_count[i] <- count2[7,"Freq"]
}
lengths <- getLength(gene_sequence)

intron_data <- data.frame(seq_names,c_count,g_count, cg_count, lengths)
intron_data$Feature <- "Introns"

#---------------------------------------------
# Make one dataframe with all Features in
#---------------------------------------------

all_data <- rbind(gene_data, exon_data, intron_data)

#---------------------------------------------
# Calculate the o/e for each gene
#---------------------------------------------
# CpG o/e = (length2/length)*(CpG count/(C count * G count))
# From Simola et al. (2013) doi:10.1101/gr.155408.113
# Althought this makes no sense as (length*length) / length is just length ???
# From Liu et al (2019) doi:10.3390/genes10020137  (L*#CpG)/(#C*#G)
# Other sources to normalise by length you do: length^2 / length -1
# Just gives the same as Liu et al. so will stick with that

all_data$cpg_ob_ex <- all_data$length * (all_data$cg_count/(all_data$c_count*all_data$g_count))
all_data <- all_data[!is.na(all_data$cpg_ob_ex),]

range(all_data$cpg_ob_ex)
plot(density(all_data$cpg_ob_ex))

# Few genes Cpg o/e > 1, we're only interested in the humps so we can
# cut the data down to remove this tail
all_data_subset <- subset(all_data, cpg_ob_ex < 2)
plot(density(all_data_subset$cpg_ob_ex))

all_data_subset <- subset(all_data_subset, cpg_ob_ex > 0.01)
plot(density(all_data_subset$cpg_ob_ex))

#---------------------------------------------
# Write out the important dataframe in case you want it later
#---------------------------------------------

write.table(all_data_subset, file=paste(species,"data.txt", sep="_"),
            sep="\t", quote = F, col.names = T, row.names = F)

#---------------------------------------------
# Calculate where the peaks are
#---------------------------------------------

# Find the value of the biggest peak for each annotation
max_gene <- which.max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$y) 
density_max_gene <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$x[max_gene] 

max_exon <- which.max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$y) 
density_max_exon <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$x[max_exon] 

max_intron <- which.max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$y) 
density_max_intron <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$x[max_intron] 

# Find the trough so we can find the second peak
DensityY_gene <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$y
DensityX_gene <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$x

DensityY_exon <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$y
DensityX_exon <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$x

DensityY_intron <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$y
DensityX_intron <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$x

# Only looking between the two peaks to avoid the tails
MinYDensity_gene<- min(DensityY_gene[DensityX_gene < 1.1 & DensityX_gene > 0.5]) 
trough_gene <- which(DensityY_gene == MinYDensity_gene) 
density_trough_gene <- DensityX_gene[trough_gene]

MinYDensity_exon<- min(DensityY_exon[DensityX_exon < 1.1 & DensityX_exon > 0.5]) 
trough_exon <- which(DensityY_exon == MinYDensity_exon) 
density_trough_exon <- DensityX_exon[trough_exon]

MinYDensity_intron<- min(DensityY_intron[DensityX_intron < 1.1 & DensityX_intron > 0.5]) 
trough_intron <- which(DensityY_intron == MinYDensity_intron) 
density_trough_intron <- DensityX_intron[trough_intron]

# Find the value of the second biggest peak
MaxY_gene <- max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])
                 $y[density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$x < density_trough_gene])
second_max_gene <- which(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$y == MaxY_gene) 
density_second_max_gene <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"])$x[second_max_gene]

MaxY_exon <- max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])
                 $y[density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$x < density_trough_exon])
second_max_exon <- which(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$y == MaxY_exon) 
density_second_max_exon <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"])$x[second_max_exon]

MaxY_intron <- max(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])
                 $y[density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$x < density_trough_intron])
second_max_intron <- which(density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$y == MaxY_intron) 
density_second_max_intron <- density(all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"])$x[second_max_intron]

#---------------------------------------------
# Make a fancy CpG o/e plot
#---------------------------------------------
pdf(paste(species,".pdf", sep=""))
ggplot(all_data_subset, aes(x=cpg_ob_ex, colour=Feature))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2, bins = 40 )+
  geom_line(stat="density",size=1.5)+
  #geom_vline(aes(xintercept=density_second_max_intron),linetype="dashed", size=1.5,colour="blue")+
  #geom_vline(aes(xintercept=density_max_intron),linetype="dashed",size=1.5,colour="blue")+
  #geom_vline(aes(xintercept=density_second_max_exon),linetype="dashed", size=1.5,colour="red")+
  #geom_vline(aes(xintercept=density_max_exon),linetype="dashed",size=1.5,colour="red")+
  #geom_vline(aes(xintercept=density_second_max_gene),linetype="dashed", size=1.5,colour="green")+
  #geom_vline(aes(xintercept=density_max_gene),linetype="dashed",size=1.5,colour="green")+
  theme_bw()+
  ggtitle(paste(species))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        plot.title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")
dev.off()