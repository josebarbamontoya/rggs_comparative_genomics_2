#!/bin/bash

for file in /Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/iqtree_gene_trees/*.fas
do
	# iqtree operation using $file
	/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/iqtree_gene_trees/iqtree -s $file -o Carcharhinidae -m MFP -nt 5
done
