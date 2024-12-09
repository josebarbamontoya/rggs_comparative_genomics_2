#!/bin/bash

#PBS -V
#PBS -N mrbayes_parallel
#PBS -q batch
#PBS -S /bin/bash
#PBS -m abe
#PBS -M jbarba@amnh.org
#PBS -e /home/jbarba/HPL_RUNS/error_file.txt
#PBS -o /home/jbarba/HPL_RUNS/output_file.txt
#PBS -l ncpus=10
#PBS -l mem=100G
#PBS -l walltime=175:00:00
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

module load beagle-2.1.2

module load mrbayes-3.2.6 

### my code lines
cd /home/jbarba/nas4/cg2_course/mrbayes_trees/phylotree_run1
mb mrbayes_phylotree.nex &

cd /home/jbarba/nas4/cg2_course/mrbayes_trees/phylotree_run2
mb mrbayes_phylotree.nex &

cd /home/jbarba/nas4/cg2_course/mrbayes_trees/timetree_run1
mb mrbayes_timetree.nex &

cd /home/jbarba/nas4/cg2_course/mrbayes_trees/timetree_run2
mb mrbayes_timetree.nex &

cd /home/jbarba/nas4/cg2_course/mrbayes_trees/timetree_prior_run1
mb mrbayes_timetree_prior.nex &

cd /home/jbarba/nas4/cg2_course/mrbayes_trees/timetree_prior_run2
mb mrbayes_timetree_prior.nex &

wait
