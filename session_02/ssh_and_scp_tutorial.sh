#!/bin/bash

################################
##### SSH and SCP Tutorial #####
##### Jose Barba ###############
################################

##### part 01 #####
##### login to remote servers #####

### NOTE: replace the amnh username `jbarba` with your own

# login to amnh huxley server
ssh jbarba@huxley-master.pcc.amnh.org

# explore the structure of the server using the `ls` and `cd` commands
cd /
ls -la
cd ~
ls -la

# logout the server
exit

# login to amnh mendel server
ssh jbarba@mendel.sdmz.amnh.org

# explore the structure of the server using the `ls` and `cd` commands
cd /
ls -la
cd ~
ls -la

# logout the server
exit

##### part 02 #####
##### transport files remotely #####

# create a directory in your home directory
cd ~
mkdir test_directory

# create two text files within `test_directory`
echo "hello :)" > test_directory/f01.txt
echo "bye :(" > test_directory/f02.txt

# display content of all .txt files within `test_directory`
cat test_directory/*.txt

### NOTE: Replace the usernames `barba` and `jbarba` with your own

#1) transport a file or directory from your laptop to huxley server:
scp -r /Users/barba/test_directory jbarba@huxley-master.pcc.amnh.org:/home/jbarba

# login to amnh huxley server
ssh jbarba@huxley-master.pcc.amnh.org

# list of files and directories
ls -la

# create a `README.txt` file within `test_directory`
touch /home/jbarba/test_directory/README.txt

# edit 'README.txt' in nano
nano /home/jbarba/test_directory/README.txt

# within nano, type in:
This directory was modified on Huxley.
# to write out file, press ^O
# to exit nano, press ^X

# display content of 'README.txt' file
cat /home/jbarba/test_directory/README.txt

# logout the server
exit

#2) transport directory from huxley to my laptop (logged in my laptop):
scp -r jbarba@huxley-master.pcc.amnh.org:/home/jbarba/test_directory /Users/barba/Desktop/
