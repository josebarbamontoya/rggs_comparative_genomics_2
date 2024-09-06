#!/bin/bash

##########################################
##### Basic Unix Navigation Tutorial #####
##### Jose Barba #########################
##########################################

# 1. show the current directory
pwd

# 2. list files and directories (including hidden files)
ls -la

# 3. change to home directory
cd ~

# 4. create a new directory called 'test_directory'
mkdir test_directory

# 5. change to 'test_directory'
cd test_directory

# 6. create three empty .txt files
touch f01.txt
touch f02.txt
touch f03.txt

# 7. list files
ls
ls -la

# 8. get help on the 'ls' command and quit
man ls
q

# 9. copy file 'f01.txt'
cp f01.txt f01_copy.txt

# 10. rename 'f01_copy.txt'
mv f01_copy.txt f01_to_delete.txt

# 11. remove 'f01_to_delete.txt'
rm f01_to_delete.txt

# 12. edit 'f01.txt' in nano
nano f01.txt

# 13. within nano, type in:
mammals birds
bat eagle
bear falcon
dolphin hawk
elephant owl
giraffe parrot
human penguin
kangaroo raven
lion sparrow
tiger swallow
whale woodpecker
# to write out file, press ^O
# to exit nano, press ^X

# 14. display the content of 'f01.txt'
cat f01.txt

# 15. rename 'f01.txt'
mv f01.txt f01_edited.txt

# 16. display the top two lines in 'f01_edited.txt' 
head -n2 f01_edited.txt

# 17. display the bottom two lines in 'f01_edited.txt' 
tail -n2 f01_edited.txt

# 18. search for the word 'lion' in 'f01_edited.txt' using grep
grep "lion" f01_edited.txt

# 19. extract the first column of 'f01_edited.txt' using awk (assuming data is space-separated)
awk '{print $1}' f01_edited.txt

# 20. find file 'f01_edited.txt' in the current directory and subdirectories
find . -name "f01_edited.txt"

# 21. change back to the home directory
cd ..

# 22. remove the 'test_directory'
rm -r test_directory
