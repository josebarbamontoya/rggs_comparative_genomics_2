############################################################
############################################################
############################################################
############################################################
rm(list=ls())
require(ape)
require(phytools)
require(dendextend)
par(mai=c(.2,.2,.2,.2))
###### sort and compare trees #####
#pdf("/Users/barba/Desktop/tree_plot.pdf", w = 15, h = 20)
a <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/simulated_dataets/ar_51_tips_true_trees/51taxa_timetree.nwk")
b <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/iqtree_concat_tree/concatenation_alignment.fas.bionj")
a <- root(a, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
b <- root(b, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
a <- ladderize(a, right = FALSE)
b <- ladderize(b, right = FALSE)
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE, use.edge.length = TRUE, location = NA)
trees <- c(a,b)
multiRF(trees)
##### norm rf distance = plain_RF/(2*(n-3)) #####
nrf <- 2/(2*(51-3))
nrf
#dev.off()
#a$edge.length <- rep(1, 100)
#b$edge.length <- rep(1, 100)
a <- force.ultrametric(a, method=c("extend"))
b <- force.ultrametric(b, method=c("extend"))
a <- as.dendrogram(a)
b <- as.dendrogram(b)
dend1 <- ladderize(a, right = FALSE)
dend2 <- ladderize(b, right = FALSE)
#dend1 <- sort(dend1)
#dend2 <- sort(dend2)
dl <- dendlist(dend1,dend2)
dl %>% untangle %>% tanglegram(common_subtrees_color_lines=FALSE, lwd=1, highlight_distinct_edges=TRUE, highlight_branches_lwd=FALSE, margin_inner=8, axes=FALSE, lab.cex=.7)
############################################################
############################################################
############################################################
############################################################
rm(list=ls())
require(ape)
require(phytools)
require(dendextend)
par(mai=c(.2,.2,.2,.2))
###### sort and compare trees #####
#pdf("/Users/barba/Desktop/tree_plot.pdf", w = 15, h = 20)
a <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/simulated_dataets/ar_51_tips_true_trees/51taxa_timetree.nwk")
b <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/iqtree_concat_tree/concatenation_alignment.fas.treefile")
a <- root(a, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
b <- root(b, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
a <- ladderize(a, right = FALSE)
b <- ladderize(b, right = FALSE)
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE, use.edge.length = TRUE, location = NA)
trees <- c(a,b)
multiRF(trees)
##### norm rf distance = plain_RF/(2*(n-3)) #####
nrf <- 0/(2*(51-3))
nrf
#dev.off()
#a$edge.length <- rep(1, 100)
#b$edge.length <- rep(1, 100)
a <- force.ultrametric(a, method=c("extend"))
b <- force.ultrametric(b, method=c("extend"))
a <- as.dendrogram(a)
b <- as.dendrogram(b)
dend1 <- ladderize(a, right = FALSE)
dend2 <- ladderize(b, right = FALSE)
#dend1 <- sort(dend1)
#dend2 <- sort(dend2)
dl <- dendlist(dend1,dend2)
dl %>% untangle %>% tanglegram(common_subtrees_color_lines=FALSE, lwd=1, highlight_distinct_edges=TRUE, highlight_branches_lwd=FALSE, margin_inner=8, axes=FALSE, lab.cex=.7)
############################################################
############################################################
############################################################
############################################################
rm(list=ls())
require(ape)
require(phytools)
require(dendextend)
par(mai=c(.2,.2,.2,.2))
###### sort and compare trees #####
#pdf("/Users/barba/Desktop/tree_plot.pdf", w = 15, h = 20)
a <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/simulated_dataets/ar_51_tips_true_trees/51taxa_timetree.nwk")
b <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/astral_species_tree/astral_tree.nwk")
a <- root(a, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
b <- root(b, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
a <- ladderize(a, right = FALSE)
b <- ladderize(b, right = FALSE)
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE, use.edge.length = TRUE, location = NA)
trees <- c(a,b)
multiRF(trees)
##### norm rf distance = plain_RF/(2*(n-3)) #####
nrf <- 0/(2*(51-3))
nrf
#dev.off()
#a$edge.length <- rep(1, 100)
b$edge.length <- rep(1, 100)
a <- force.ultrametric(a, method=c("extend"))
b <- force.ultrametric(b, method=c("extend"))
a <- as.dendrogram(a)
b <- as.dendrogram(b)
dend1 <- ladderize(a, right = FALSE)
dend2 <- ladderize(b, right = FALSE)
#dend1 <- sort(dend1)
#dend2 <- sort(dend2)
dl <- dendlist(dend1,dend2)
dl %>% untangle %>% tanglegram(common_subtrees_color_lines=FALSE, lwd=1, highlight_distinct_edges=TRUE, highlight_branches_lwd=FALSE, margin_inner=8, axes=FALSE, lab.cex=.7)
############################################################
############################################################
############################################################
############################################################
rm(list=ls())
require(ape)
require(phytools)
require(dendextend)
par(mai=c(.2,.2,.2,.2))
###### sort and compare trees #####
#pdf("/Users/barba/Desktop/tree_plot.pdf", w = 15, h = 20)
a <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/simulated_dataets/ar_51_tips_true_trees/51taxa_timetree.nwk")
b <- read.tree(file = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/mrbayes_trees/phylotree_run1/mrbayes_phylotree.nex.con.tre")
a <- root(a, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
b <- root(b, outgroup=c("Carcharhinidae"), resolve.root=TRUE)
a <- ladderize(a, right = FALSE)
b <- ladderize(b, right = FALSE)
comparePhylo(a, b, plot = TRUE, force.rooted = TRUE, use.edge.length = TRUE, location = NA)
trees <- c(a,b)
multiRF(trees)
##### norm rf distance = plain_RF/(2*(n-3)) #####
nrf <- XYZ/(2*(51-3))
nrf
#dev.off()
#a$edge.length <- rep(1, 100)
b$edge.length <- rep(1, 100)
a <- force.ultrametric(a, method=c("extend"))
b <- force.ultrametric(b, method=c("extend"))
a <- as.dendrogram(a)
b <- as.dendrogram(b)
dend1 <- ladderize(a, right = FALSE)
dend2 <- ladderize(b, right = FALSE)
#dend1 <- sort(dend1)
#dend2 <- sort(dend2)
dl <- dendlist(dend1,dend2)
dl %>% untangle %>% tanglegram(common_subtrees_color_lines=FALSE, lwd=1, highlight_distinct_edges=TRUE, highlight_branches_lwd=FALSE, margin_inner=8, axes=FALSE, lab.cex=.7)
############################################################
############################################################
############################################################
############################################################


