#!/usr/bin/env Rscript

########################################################
##### Introduction to R for Phylogenomics Tutorial #####
##### Jose Barba #######################################
########################################################

##### part 01 #####
##### data manipulation, analysis and visualization #####

# clear workspace
rm(list=ls())
library(ggplot2)

# set working directory
### NOTE: substitute path with your own
setwd("/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05/ncbi_refseq_genome_statistics")
# print your current working 
getwd()

# download a summary of current NCBI RefSeq genome assemblies
system("wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt")
system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt")

# read the file to edit the table
lines <- readLines("assembly_summary_refseq.txt")

# remove the '#' symbol
cleaned_lines <- gsub("#", "", lines)

# convert cleaned lines to a data frame, skip first line, and use quote = "" to avoid issues with unmatched quotes
g_data <- read.table(text = cleaned_lines, sep = "\t", skip = 1, header = TRUE, fill = TRUE, quote = "")

# export the edited data frame to a new file
write.table(g_data, "assembly_summary_refseq_cleaned.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# count number of columns (variables and rows (geenomes)
nrow(g_data)
ncol(g_data)

# count number of genemoes for each major lineage
genome_counts <- table(g_data$group)
# display the result
print(genome_counts)

##### summary statistics of gc percent #####
summary(g_data$gc_percent)

# statistical messures of gc percent
mean(g_data$gc_percent)
median(g_data$gc_percent)
min(g_data$gc_percent)
max(g_data$gc_percent)

##### create a histogram of gc content #####
hist(g_data$gc_percent, main="GC content", col="red3")

##### create a barplot of genome size #####
# sort genome sizes in decreasing order and preserve the order of organisms
g_data <- g_data[order(g_data$genome_size, decreasing = TRUE), ]
# create the bar plot
barplot(g_data$genome_size[100:37000], main="genome size", xlab="organism", ylab="base pairs")

##### create a stacked barplot of coding and noncoding genes #####
# convert columns to numeric
g_data$protein_coding_gene_count <- as.numeric(as.character(g_data$protein_coding_gene_count))
g_data$non_coding_gene_count <- as.numeric(as.character(g_data$non_coding_gene_count))

# remove rows with na values
g_data_clean <- na.omit(g_data)

# create a data matrix subset
data_matrix <- rbind(g_data_clean$non_coding_gene_count[1:300], g_data_clean$protein_coding_gene_count[1:300]) # NOTE: the maximum number of base pairs is 40054324612 (adjust y-axis as necessary) 

# create the stacked barplot for the subset
barplot(data_matrix, col=c("red3", "blue3"), beside=FALSE, legend.text=c("Non-coding genes", "Protein-coding genes"), main="Number of protein-coding and mon-coding genes per organism", cex.axis=1, cex.names=1, xlab="organism", ylab="genes")

##### create a scatterplot of genome size and number of genes #####
# ensure the necessary packages are loaded
library(ggplot2)

# clean the data to remove NA values
cleaned_df <- na.omit(data.frame(genome_size = as.numeric(g_data$genome_size), total_gene_count = as.numeric(g_data$total_gene_count)))

# create the scatterplot using ggplot2 for better aesthetics
ggplot(cleaned_df, aes(x = genome_size, y = total_gene_count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Genome Size vs Number of Genes",
       x = "Genome Size (base pairs)", 
       y = "Total Gene Count") +
  theme_minimal()

# compute correlation
cor(cleaned_df$total_gene_count, cleaned_df$genome_size, method = c("pearson"))

# liniear regresion through the origin
fit <- lm(cleaned_df$total_gene_count ~ cleaned_df$genome_size + 0)
fit
summary(fit)

##### create a boxplot of genome size for plant and fungi #####
# sort table by group
sorted_table <- g_data[order(g_data$group), ]

# collect plant geneome sizw
rf1 <- (genome_size=c(sorted_table$genome_size[381746:381931]))

# collect fungi genemoe size
rf2 <- (genome_size=c(sorted_table$genome_size[380685:381316]))

# create a boxplot
require(vioplot)
boxplot(rf1, rf2, names = c("plant", "fungi"), ylab = "genome size (bp)", main = "Plant and fungi genome size")

##### create a density plot of percentage difference between non-coding and protein-coding gene counts #####
cleaned_df2 <- na.omit(data.frame(protein_coding_gene_count = g_data$protein_coding_gene_count, non_coding_gene_count = g_data$non_coding_gene_count))

# extract non-coding and protein-coding gene counts
nc <- cleaned_df2$non_coding_gene_count
pc <- cleaned_df2$protein_coding_gene_count

# ensure both vectors are numeric
nc <- as.numeric(nc)
pc <- as.numeric(pc)

# calculate delta_genes as a percentage difference
delta_genes <- ((pc-nc)/nc)*100

# check for NAs in delta_genes after calculation
cat("Number of NAs in delta_genes before cleaning:", sum(is.na(delta_genes)), "\n")

# remove NA values from delta_genes
delta_genes <- na.omit(delta_genes)

# check for NAs again after removal
cat("Number of NAs in delta_genes after cleaning:", sum(is.na(delta_genes)), "\n")

# calculate density of delta_genes
d1 <- density(delta_genes)

# plot the density
plot(d1, main = "Protein-coding and non-coding gene percent difference", xlab = "percent difference", ylab = "density", xlim=c(0, 100000), lwd=2 , col="red3")

######################################################################################
######################################################################################
######################################################################################

##### part 02 #####
##### phylogenomic analysis #####

rm(list=ls())
require(ape)
require(phytools)
require(dendextend)
require(phangorn)
par(mai=c(.3,.3,.3,.3))

# set working directory
### NOTE: substitute path with your own
setwd("/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05")
# print your current working 
getwd()

# download a fasta msa of 44 canid mt genomes
system("wget https://raw.githubusercontent.com/josebarbamontoya/rggs_comparative_genomics_2/main/session_05/44canid_mt_genomes.fas")

# read fasta msa
aln <- read.FASTA("/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05/44canid_mt_genomes.fas", type ="DNA")

##### construct a tree using the nj method #####
# calculate distance matrix
dist_matrix <- dist.dna(aln, model = "JC69")

# construct a tree using the nj method
nj_tree <- nj(dist_matrix)

# ladderize tree
nj_tree <- ladderize(nj_tree, right = FALSE)

# plot tree
plot(nj_tree, main = "neighbor-joining tree")
nodelabels()
edgelabels()

# save newick tree
write.tree(nj_tree, file = "nj_tree.nwk")

##### construct a maximum likelihood tree #####
# create phydat object
phy_data <- phyDat(aln, type = "DNA")
fit <- pml(nj_tree, data = phy_data)

# optimize tree
fit <- optim.pml(fit, model = "GTR", optGamma = TRUE)

# plot the ML tree
ml_tree <- fit$tree
plot(ml_tree, main = "maximum likelihood tree")
nodelabels()
edgelabels()

# save newick tree
write.tree(ml_tree, file = "ml_tree.nwk")

##### sort and compare trees #####
#####  comparephylo and robinson-foulds distance #####
a <- read.tree(file = "nj_tree.nwk")
b <- read.tree(file = "ml_tree.nwk")
a <- root(a, outgroup=c("GF_GrayFox_Vermont"), resolve.root=TRUE)
b <- root(b, outgroup=c("GF_GrayFox_Vermont"), resolve.root=TRUE)
a <- ladderize(a, right = FALSE)
b <- ladderize(b, right = FALSE)
#a <- force.ultrametric(a, method=c("extend"))
#b <- force.ultrametric(b, method=c("extend"))
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE, use.edge.length = TRUE, location = NA)
trees <- c(a,b)
multiRF(trees)
##### norm rf distance = plain_RF/(2*(n-3))
nrf <- 0/(2*(44-3))
nrf

##### tanglegram #####
a <- read.tree(file = "nj_tree.nwk")
b <- read.tree(file = "ml_tree.nwk")
a <- root(a, outgroup=c("GF_GrayFox_Vermont"), resolve.root=TRUE)
b <- root(b, outgroup=c("GF_GrayFox_Vermont"), resolve.root=TRUE)
#a$edge.length <- rep(1, 86)
#b$edge.length <- rep(1, 86)
a <- force.ultrametric(a, method=c("extend"))
b <- force.ultrametric(b, method=c("extend"))
a <- as.dendrogram(a)
b <- as.dendrogram(b)
dend1 <- ladderize(a, right = FALSE)
dend2 <- ladderize(b, right = FALSE)
dl <- dendlist(dend1,dend2)
dl %>% untangle %>% tanglegram(common_subtrees_color_lines=FALSE, lwd=1, highlight_distinct_edges=TRUE, highlight_branches_lwd=FALSE, margin_inner=15, axes=FALSE, lab.cex=1)

##### scatterplot of nj and ml branch lengths #####
# make branch length dataframes
nj_tree <- ladderize(nj_tree, right = FALSE)
ml_tree <- ladderize(ml_tree, right = FALSE)
nj_bl <- nj_tree$edge.length
ml_bl <- ml_tree$edge.length

# make scatterplot
par(mai=c(1,1,1,1))
plot(nj_bl, ml_bl, main = "NJ and ML tree branch length comparison")
abline(lm(ml_bl ~ nj_bl + 0), col="blue3")

# compute correlation
cor(ml_bl, nj_bl, method = c("pearson"))

# liniear regresion through the origin
fit <- lm(ml_bl ~ nj_bl + 0)
fit
summary(fit)
