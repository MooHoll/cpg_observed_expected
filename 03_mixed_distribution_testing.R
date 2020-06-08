#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(echo=TRUE)

# Takes one argument, the cpg o/e data file generated in 02_graph_cpgOE_per_feature.R

#---------------------------------------------
# Use mixtools to identify statistical mixed distributions
#---------------------------------------------

library(readr)
library(seqinr)
library(stringr)
library(sqldf)
library(ggplot2)
library(ggpmisc)
library(mixtools)

species <- tools::file_path_sans_ext(basename(args[1]))
species <- gsub("_data","",species)

all_data_subset <- read_delim(args[1], "\t", escape_double = FALSE, trim_ws = TRUE)

#---------------------------------------------
# Guassian modelling to see if statistically we might have two peaks
#---------------------------------------------
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

cpg_oe_genes <- all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Genes"]
cpg_oe_exons <- all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Exons"]
cpg_oe_introns <- all_data_subset$cpg_ob_ex[all_data_subset$Feature=="Introns"]

#---------------------------------------------

mixmdl <- normalmixEM(cpg_oe_genes, k = 2)

pdf(paste(species,"genes.pdf", sep="_"))
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", alpha=0.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  geom_vline(aes(xintercept=mixmdl$mu[1]),linetype="dashed", size=1.5,colour="red")+
  geom_vline(aes(xintercept=mixmdl$mu[2]),linetype="dashed",size=1.5,colour="blue")+
  theme_bw()+
  ggtitle(paste(species, ":genes", sep=""))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        plot.title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")
dev.off()

#---------------------------------------------

mixmdl <- normalmixEM(cpg_oe_exons, k = 2)

pdf(paste(species,"exons.pdf", sep="_"))
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", alpha=0.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  geom_vline(aes(xintercept=mixmdl$mu[1]),linetype="dashed", size=1.5,colour="red")+
  geom_vline(aes(xintercept=mixmdl$mu[2]),linetype="dashed",size=1.5,colour="blue")+
  theme_bw()+
  ggtitle(paste(species, ":exons", sep=""))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        plot.title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")
dev.off()

#---------------------------------------------

mixmdl <- normalmixEM(cpg_oe_introns, k = 2)

pdf(paste(species,"introns.pdf", sep="_"))
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.05, colour = "black", alpha=0.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  geom_vline(aes(xintercept=mixmdl$mu[1]),linetype="dashed", size=1.5,colour="red")+
  geom_vline(aes(xintercept=mixmdl$mu[2]),linetype="dashed",size=1.5,colour="blue")+
  theme_bw()+
  ggtitle(paste(species, ":introns", sep=""))+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        plot.title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")
dev.off()

#---------------------------------------------

# This gives confidence intervals for the mu
# Tried 100 and 1000 bootstraps, makes basically no difference
#bootstap <- boot.se(mixmdl, B=100)
#bootstap$mu.se[1] # 0.001018013
#quantile(bootstap$mu[1,],c(0.25,0.975)) #0.7627601 0.7654038
#bootstap$mu.se[2] # 0.003746412
#quantile(bootstap$mu[2,],c(0.25,0.975)) #0.8176026 0.8279798 

# I would prefer confidence intervals for the mean of each distribution 
# as in Bewick et al. (2016), but can't find out how anywhere ...