#!/bin/bash

#PBS -V
#PBS -N mcmctree_parallel
#PBS -q batch
#PBS -S /bin/bash
#PBS -m abe
#PBS -M jbarba@amnh.org
#PBS -e /home/jbarba/HPL_RUNS/error_file.txt
#PBS -o /home/jbarba/HPL_RUNS/output_file.txt
#PBS -l ncpus=5
#PBS -l mem=10G
#PBS -l walltime=48:00:00
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

### my code lines
cd /home/jbarba/nas4/canid_evolution/canid_mtdna_mcmctree_huxley/mcmc01
~/bin/mcmctree mcmctree.ctl > screen.txt &

cd /home/jbarba/nas4/canid_evolution/canid_mtdna_mcmctree_huxley/mcmc02
~/bin/mcmctree mcmctree.ctl > screen.txt &

cd /home/jbarba/nas4/canid_evolution/canid_mtdna_mcmctree_huxley/mcmc03
~/bin/mcmctree mcmctree.ctl > screen.txt &

cd /home/jbarba/nas4/canid_evolution/canid_mtdna_mcmctree_huxley/mcmc04
~/bin/mcmctree mcmctree.ctl > screen.txt &

wait
