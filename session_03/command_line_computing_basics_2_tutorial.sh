#!/bin/bash

###########################################
##### Command Line Computing Basics 2 #####
##### Jose Barba ##########################
###########################################

##### part 01 #####
##### shell pipelines #####

# change to home directory
cd ~

# download mammal_data_analysis.zip file from github
wget https://raw.githubusercontent.com/josebarbamontoya/rggs_comparative_genomics_2/main/session_03/mammal_data_analysis.zip
# use curl if wget is not available
#curl -0 https://raw.githubusercontent.com/josebarbamontoya/rggs_comparative_genomics_2/main/session_03/mammal_data_analysis.zip

# uncompress `mammal_data_analysis.zip` using `unzip`
unzip mammal_data_analysis.zip

# or uncompress `mammal_data_analysis.zip` using `tar`
tar -xf mammal_data_analysis.zip

# change to `mammal_data_analysis` directory
cd mammal_data_analysis

# list alignemtns and count the number of alignments
ls *.fas | wc -l

# count the total number of sequences for each FASTA file 
grep -c "^>" *.fas | sed 's/:/ - /'

# check if all the alignments contain the outgroup ornithorhynchus_anatinus
grep "ornithorhynchus_anatinus" *.fas | nl

# search for repetitive segments in all the alignments
grep "AAAAAAAA" *.fas  | nl

# compute the GC content of each .fas file and print the result
ls *.fas | xargs -I {} sh -c 'echo "{}: "; grep -v "^>" "{}" | awk "{gc+=gsub(/[GCgc]/,\"\"); total+=length(\$0)} END {if (total > 0) print (gc/total)*100 \"%\"; else print \"0%\"}"'

##### part 02 #####
##### loops in shell #####

##### for loop examples #####

# for loop: iterate over a list of values
for var in item1 item2 item3
do
    command $var
done

# example1: print each filename in the current directory
for file in *
do
    echo $file
done

# example2: perform an action on multiple files
for file in *.txt
do
    mv "$file" "${file%.txt}.bak"
done

##### while loop examples #####

# while loop: execute as long as a condition is true
while [ condition ]
do
    command
done

# example1: count down from 5 #####
count=5
while [ $count -gt 0 ]
do
    echo $count
    count=$((count - 1))
done

# example2: prepend line numbers to each line in a text file
count=1
while IFS= read -r line
do
    echo "$count: $line" >> m2778_subsampled_gb_lines_numbered.txt
    count=$((count + 1))
done < m2778_subsampled_gb.fas

###### for loop exercises #####

# 1) count the total number of sequences for 5 selected fasta files within mammal_data_analysis directory
for file in m2778_subsampled_gb.fas m3085_subsampled_gb.fas m4281_subsampled_gb.fas m11337_subsampled_gb.fas m14339_subsampled_gb.fas
do
  count=$(grep -c "^>" "$file")
  echo "$file: $count sequences"
done

# 2) count the total number of sequences for each fasta file within mammal_data_analysis directory
for file in *.fas; do
  count=$(grep -c "^>" "$file")
  echo "$file: $count sequences"
done

# 3) iqtree loop mammal_data_analysis 
cp mammal_data_analysis mammal_data_analysis_gtr
cp mammal_data_analysis mammal_data_analysis_hky

cd mammal_data_analysis_gtr
for file in /Users/barba/mammal_data_analysis_gtr/*.fas
do
  # iqtree operation using $file
  /Users/barba/mammal_data_analysis_gtr/iqtree -s $file -m GTR+G5 -nt 2
done

cd mammal_data_analysis_hky
for file in /Users/barba/mammal_data_analysis_hky/*.fas
do
  # iqtree operation using $file
  /Users/barba/mammal_data_analysis_hky/iqtree -s $file -m HKY+G5 -nt 2
done

cd mammal_data_analysis_gtr
# extract likelihood from iqtree info files
grep "BEST SCORE FOUND :" *.log > iqtree_tree_likelihoods_gtr.txt
# display the content of 'iqtree_tree_likelihoods_gtr.txt' 
cd ..

cd mammal_data_analysis_hky
# extract likelihood from iqtree info files
grep "BEST SCORE FOUND :" *.log > iqtree_tree_likelihoods_hky.txt
# display the content of 'iqtree_tree_likelihoods_hky.txt' 
cd ..

# compare tree likelihoods from diffrent models
cat mammal_data_analysis_gtr/iqtree_tree_likelihoods_gtr.txt
cat mammal_data_analysis_hky/iqtree_tree_likelihoods_hky.txt

# 4) find a list of strings, loop through files in a directory, and create a prescence list
# create a script with the loop code below like `find_string_loop_files_make_presence_list.sh`
# define the three words you want to search for
word1="notamacropus_eugenii"
word2="sarcophilus_harrisii"
word3="monodelphis_domestica"

# directory where your files are located
directory="/Users/barba/mammal_data_analysis"

# create an empty file to store the list of matching files
output_file="matching_files.txt"
> "$output_file"

# loop through each file in the directory
for file in "$directory"/*; do
    # check if the file contains any of the three words
    if grep -qE "$word1|$word2|$word3" "$file"; then
        # if it contains any of the words, append the filename to the output file
        echo "$file" >> "$output_file"
    fi
done

# print the list of matching files
cat "$output_file"

###### while loop exercises #####

# 1) print sequence headers
while IFS= read -r line
do
    if [[ $line == ">"* ]]; then
        echo "$line"  # print the header
    fi
done < m2778_subsampled_gb.fas

# 2) count the number of bases in each sequence (ignoring headers)
sequence=""
while IFS= read -r line
do
    if [[ $line == ">"* ]]; then
        if [[ -n $sequence ]]; then
            length=${#sequence}
            echo "Total sequence length: $length" >> output.txt
            sequence=""  # reset sequence for the next one
        fi
        echo "$line" >> output.txt  # print the header to output.txt
    else
        sequence+="$line"  # concatenate sequence lines
    fi
done < m2778_subsampled_gb.fas

# process the last sequence in the file
if [[ -n $sequence ]]; then
    length=${#sequence}
    echo "Total sequence length: $length" >> output.txt
fi

##### part 03 #####
##### combining loops and pipelines #####

# 1) search for occurrences of "homo_sapiens" in all .fas files
for file in *.fas; do
    count=$(grep "homo_sapiens" "$file" | wc -l)
    echo "$file: $count"
done

# 2) extract ornithorhynchus_anatinus sequence from all the alignments:
awk 'BEGIN { FS="\n"; RS=">"; ORS="" } NR > 1 { header = $1; seq = ""; for (i = 2; i <= NF; i++) { seq = seq $i; } print ">" header "\n" seq "\n" FILENAME; }' *.fas | sed 's/ //g' | grep -A 1 'ornithorhynchus_anatinus'

# 3) calculate GC Content for each sequence
sequence=""
while IFS= read -r line
do
    if [[ $line == ">"* ]]; then
        if [[ -n $sequence ]]; then
            gc_count=$(echo "$sequence" | grep -o "[GCgc]" | wc -l)  # count Gs and Cs
            total_count=${#sequence}  # get total number of bases
            gc_percentage=$(echo "scale=2; ($gc_count/$total_count)*100" | bc)
            echo "GC Content: $gc_percentage%" >> gc_output.txt
            sequence=""  # reset sequence for the next one
        fi
        echo "$line" >> gc_output.txt  # Print the header to output.txt
    else
        sequence+="$line"  # append sequence lines
    fi
done < m2778_subsampled_gb.fas

# for the last sequence in the file
if [[ -n $sequence ]]; then
    gc_count=$(echo "$sequence" | grep -o "[GCgc]" | wc -l)
    total_count=${#sequence}
    gc_percentage=$(echo "scale=2; ($gc_count/$total_count)*100" | bc)
    echo "GC Content: $gc_percentage%" >> gc_output.txt
fi
