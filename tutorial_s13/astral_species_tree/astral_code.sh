#!/bin/bash

#change dir
cd /Users/barba/Desktop/phylogenetic_inference_tutorial/astral_species_tree

#run astral with local posterior prob (lpp)
java -jar astral.5.7.8.jar -i 100_ml_tres.nwk -o out.nwk 2>out.log
