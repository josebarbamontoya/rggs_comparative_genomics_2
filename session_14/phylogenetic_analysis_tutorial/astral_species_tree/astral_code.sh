#!/bin/bash

#change dir
cd /Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_14/phylogenetic_analysis_tutorial/astral_species_tree

#run astral with local posterior prob (lpp)
java -jar astral.5.7.8.jar -i 100_ml_iqtree_tres.nwk -o astral_tree.nwk 2>out.log
