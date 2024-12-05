#!/bin/bash

#change directory
cd /Users/barba/Desktop/phylogenetic_inference_tutorial/iqtree_concat_tree

# find best fit model and infer ml tree
iqtree -s concatenation_alignment.fas -o Carcharhinidae -m MFP -b 10 -nt 10

# other iqtree analyses 
#iqtree -s 101taxa.fasta -o outgrop_name -b 100 -m GTR+G4 -nt 10
#iqtree -s 20taxa.fasta -m HKY+G4 -nt 10
