#!/bin/bash

######################################
##### HPC and Scripting Tutorial #####
##### Jose Barba #####################
######################################

##### part 01 #####
##### create a shel script with nano #####

# 1. open nano
nano

# 2. type in the following code (leave the last line empty):
#!/bin/bash
for file in *
do
    echo $file
done

# 3. press ^O to save file
# 4. save file as loop_script.sh and press return 
# 5. to exit nano, press ^X
# 6. modify the file permissions to make the script executable
chmod 777 loop_script.sh 

##### part 02 #####
##### set up the hpc tutorial directory #####

# login to huxley server
### NOTE: substitute the amnh_username with your own
ssh amnh_username@huxley-master.pcc.amnh.org

# copy the hpc_tutorial directory to your nas4 directory
### NOTE: substitute the amnh_username with your own
cp -r /home/jbarba/nas4/cg2_course/hpc_tutorial /home/amnh_username/nas4/

# change into the hpc_tutorial directory
### NOTE: substitute the amnh_username with your own
cd /home/amnh_username/nas4/hpc_tutorial

##### part 03 #####
##### software on huxley #####

# login to huxley server
### NOTE: substitute the amnh_username with your own
ssh amnh_username@huxley-master.pcc.amnh.org

# check software installed on huxley
module avail

# load software
module load RAxML-8.2.11

# confirm the software was loaded
which raxml

# unload software
module unload RAxML-8.2.11

######################################################################################################################################################
##### part 04 #####
##### pbs serial blast #####
### NOTE: substitute the amnh_username and path with your own information

# 1. blast setup:

cd /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search
module load ncbi-blast-2.6.0+
makeblastdb -in GRCh38_latest_genomic.fna -parse_seqids -dbtype nucl

# 2. edit the script below and create a shell script named `pbs_blast_serial.sh` and sav it within `blast_search` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -S /bin/bash
#PBS -N serial-blast
#PBS -l ncpus=1
#PBS -l walltime=00:10:00
cd $PBS_O_WORKDIR
module load ncbi-blast-2.6.0+
date
blastn -query /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/dicty_rna.fa -db /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/GRCh38_latest_genomic.fna -out /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/serial.blastout
date

##### script ends here #####

# 3. submit job
qsub pbs_blast_serial.sh

# 4 check jobs status
qstat

# 5 kill a job
qdel <job_id>

######################################################################################################################################################
##### part 05 #####
##### pbs array blast #####
### NOTE: substitute the amnh_username and path with your own information

# 1. blast setup:

cd /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search
module load ncbi-blast-2.6.0+
makeblastdb -in GRCh38_latest_genomic.fna -parse_seqids -dbtype nucl

# 2. edit the script below and create a shell script named `pbs_blast_array.sh` and sav it within `blast_search` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -S /bin/bash
#PBS -N array-blast
#PBS -l ncpus=1
#PBS -l walltime=00:10:00
#PBS -J 1-3
cd $PBS_O_WORKDIR
module load ncbi-blast-2.6.0+
date
blastn -query /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/dicty_rna.fa.$PBS_ARRAY_INDEX -db /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/GRCh38_latest_genomic.fna -out /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search/serial.blastout.$PBS_ARRAY_INDEX
date

##### script ends here #####

# 3. submit job
qsub pbs_blast_array.sh

# 4 check jobs status
qstat

# 5 kill a job
qdel <job_id>

######################################################################################################################################################
##### part 06 #####
##### pbs threaded raxml #####
### NOTE: substitute the amnh_username and path with your own information

# 1. change into the hpc_tutorial directory
cd /home/jbarba/nas4/cg2_course/hpc_tutorial/raxml_analysis

# 2. edit the script below and create a shell script named `pbs_raxml_threads.sh` and sav it within `raxml_analysis` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -S /bin/bash
#PBS -N threads-raxml
#PBS -l ncpus=10
#PBS -l walltime=00:10:00
cd $PBS_O_WORKDIR
module load RAxML-8.2.11
date
raxml -T 10 -m PROTGAMMAJTTF -f d -N 1 -n best -p 1977 -s /home/jbarba/nas4/cg2_course/hpc_tutorial/raxml_analysis/data/A0JC77.phy
date

##### script ends here #####

# 3. submit job
qsub pbs_blast_array.sh

# 4 check jobs status
qstat

# 5 kill a job
qdel <job_id>

######################################################################################################################################################
##### part 07 #####
##### pbs mpi abyss #####
### NOTE: substitute the amnh_username and path with your own information

# 1. change into the hpc_tutorial directory
cd /home/jbarba/nas4/cg2_course/hpc_tutorial/abyss_genome_assembling

# 2. edit the script below and create a shell script named `pbs_abyss_mpi.sh` and sav it within `abyss_genome_assembling` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash
#PBS -V
#PBS -q batch
#PBS -S /bin/bash
#PBS -N abyss-mpi
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=1:00:00
cd $PBS_O_WORKDIR
module load openmpi-2.1.1 abyss-2.0.2
date
abyss-pe np=16 k=75 v=-v se=/home/jbarba/nas4/cg2_course/hpc_tutorial/abyss_genome_assembling/27924_S15_L001_R2_001.fastq name=test "unitigs"
date

##### script ends here #####

# 3. submit job
qsub pbs_blast_array.sh

# 4 check jobs status
qstat

# 5 kill a job
qdel <job_id>

######################################################################################################################################################
##### part 07 #####
##### gnu parallel #####
### NOTE: substitute the amnh_username and path with your own information

# 1. change into the hpc_tutorial directory
cd /home/jbarba/nas4/cg2_course/hpc_tutorial/gnu_parallel_blast

# 2. edit the script below and create a shell script named `gnu_parallel_blast.sh` and sav it within `gnu_parallel_blast` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash

cd /home/jbarba/nas4/cg2_course/hpc_tutorial/blast_search
module load ncbi-blast-2.6.0+
makeblastdb -in GRCh38_latest_genomic.fna -parse_seqids -dbtype nucl

module load parallel-20171122
cd /home/jbarba/nas4/cg2_course/hpc_tutorial/gnu_parallel_blast

# compress files in the background
for i in {1..10}; do gzip dicty_genome.$i.fa & done
wait  # Wait for all background processes to finish
# decompress all gzipped files
gzip -d *.gz
# parallel gzip using GNU Parallel with 2 jobs at a time
for i in {1..10}; do
  echo "Compressing dicty_genome.$i.fa"
  sem -j2 gzip dicty_genome.$i.fa
done

sem --wait  # wait for all parallel jobs to finish

parallel -j3 'blastn -query ../blast_search/dicty_rna.fa.{} -db ../blast_search/GRCh38_latest_genomic.fna | grep -c "No hits found" > {}.out' ::: {1..3}

##### script ends here #####

# 3. execute the shell script
./gnu_parallel_blast.sh

# 4 kill a job
kill <job_id>

######################################################################################################################################################
##### part 07 #####
##### autogenerate scripts #####
### NOTE: substitute the amnh_username and path with your own information

# 1. change into the hpc_tutorial directory
cd /home/jbarba/nas4/cg2_course/hpc_tutorial/auto_gen_scripts

# 2. edit the script below and create a shell script named `auto_gen_scripts.sh` and sav it within `auto_gen_scripts` directory
# use nano or any other text editor

##### script begins here #####

#!/bin/bash

for i in `cat list.data`
do 
    mkdir $i
    cd $i 
    /home/jbarba/nas4/cg2_course/hpc_tutorial/auto_gen_scripts/generic_pbs.pl -l ncpus=10 -w walltime=1:00:00 -m RAxML-8.2.11 -n $i -b 'raxml -T 10 -m PROTGAMMAJTTF -f d -N 1 -n '$i' -p 1977 -s /home/jbarba/nas4/cg2_course/hpc_tutorial/auto_gen_scripts/data/'$i''
    cd ../ 
done

##### script ends here #####

# 3. execute the shell script
./auto_gen_scripts.sh

# 4 kill a job
kill <job_id>
