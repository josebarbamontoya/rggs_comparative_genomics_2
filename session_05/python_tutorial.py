#!/usr/bin/python3

#############################################################
##### Introduction to Python for Phylogenomics Tutorial #####
##### Jose Barba ############################################
#############################################################

##### part 01 #####
##### data manipulation, analysis and visualization #####

import os
import subprocess
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.stats import pearsonr
import seaborn as sns
from scipy.stats import gaussian_kde

# set working directory
# NOTE: Substitute path with your own
working_directory = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05/ncbi_refseq_genome_statistics"
os.chdir(working_directory)

# display current working directory
os.getcwd()

# download files using wget
#subprocess.run("wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt", shell=True)
#subprocess.run("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt", shell=True)

# open and read the file
with open("assembly_summary_refseq.txt", "r") as file:
    g_data = file.read()

# remove the '#' symbol from each line
cleaned_lines = [line.replace("#", "") for line in g_data]

# convert cleaned lines to a dataframe, skipping the first line
# use stringio to read the cleaned lines as a cvs-like structure
data_string = "".join(cleaned_lines)
g_data = pd.read_csv(StringIO(data_string), sep="\t", skiprows=1, quoting=3, engine='python')

# export the edited dataframed to a new file
g_data.to_csv("assembly_summary_refseq_cleaned.txt", sep="\t", index=False, quoting=3)

# read the cleaned dataframe with low_memory set to false
g_data = pd.read_csv("assembly_summary_refseq_cleaned.txt", sep="\t", low_memory=False)

# count number of rows (genomes) and columns (variables) 
# number of rows
num_rows = g_data.shape[0]
# number of columns
num_cols = g_data.shape[1]
# print numbers
print(f"Number of genomes: {num_rows}")
print(f"Number of variables: {num_cols}")

# count number of genomes for each major lineage
genome_counts = g_data['group'].value_counts()
print("\nGenome counts for each major lineage:")
print(genome_counts)

##### summary statistics of gc percent #####
print("\nSummary statistics for GC percent:")
gc_percent_summary = g_data['gc_percent'].describe() 
print(gc_percent_summary)

# statistical measures of gc percent
mean_gc = g_data['gc_percent'].mean()
median_gc = g_data['gc_percent'].median()
min_gc = g_data['gc_percent'].min()
max_gc = g_data['gc_percent'].max()

# print statistics
print(f"\nMean GC percent: {mean_gc}")
print(f"Median GC percent: {median_gc}")
print(f"Min GC percent: {min_gc}")
print(f"Max GC percent: {max_gc}")

##### create a histogram of gc content #####
plt.figure(figsize=(10, 6))
plt.hist(g_data['gc_percent'], bins=30, color='blue', edgecolor='black')
plt.title('GC Content Histogram')
plt.xlabel('GC Percentage')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.show()

# compute maximum number of base pairs
g_data['genome_size'].max()

# compute average number of base pairs
g_data['genome_size'].mean()

##### create a barplot of genome size #####
g_data_sorted = g_data.sort_values(by='genome_size', ascending=False)
plt.figure(figsize=(20, 6))
plt.bar(g_data_sorted['organism_name'], g_data_sorted['genome_size'], color='green')
plt.title('Genome Size Barplot')
plt.xlabel('Organism')
plt.ylabel('Base Pairs')
plt.ylim(0, 100000000) # NOTE: the maximum number of base pairs is 40054324612 (adjust y-axis as necessary) 
plt.xticks([], rotation=45) # hide x axis labels
plt.grid(axis='y', alpha=0.75)
plt.tight_layout()
plt.show()

##### create a stacked barplot of coding and noncoding genes #####
##### create the stacked barplot for a subset #####
# ensure the columns are numeric, replacing non-numeric entries with NaN and dropping rows with missing data
g_data_clean = g_data[['protein_coding_gene_count', 'non_coding_gene_count']].apply(pd.to_numeric, errors='coerce').dropna()
# compute maximum number protein coding genes 
g_data_clean['protein_coding_gene_count'].max()
# select rows 10 to 1000 (python uses zero-based indexing, so this is equivalent to rows 1 to 1000 in r)
subset = g_data_clean.iloc[10:1000] # NOTE: the total number of refseq genomes is 397713 (adjust set as necessary) 
# create the data matrix (stacked values for each organism in the subset)
data_matrix = np.array([subset['protein_coding_gene_count'], subset['non_coding_gene_count']])
# create the stacked bar plot
fig, ax = plt.subplots(figsize=(10, 6))
organisms = np.arange(subset.shape[0])
# plot the stacked bars (non-coding genes first, then protein-coding genes stacked on top)
bar1 = ax.bar(organisms, data_matrix[1], color='red', label='Non-coding genes')
bar2 = ax.bar(organisms, data_matrix[0], bottom=data_matrix[1], color='blue', label='Protein-coding genes')
ax.set_title("Number of Protein-coding and Non-coding Genes (geneomes 1 to 1000)")
ax.set_ylabel("Gene Count")
### set dynamic y-limit based on data
#ax.set_ylim(0, data_matrix.sum(axis=0).max() * 1.1)
#### set hard y-limit
plt.ylim(0, 45000) # NOTE: the total number of protein_coding_gene_count 103787 (adjust set as necessary)
plt.xlabel('Organism')
ax.legend()
plt.xticks([], rotation=45) # hide x axis labels
plt.tight_layout()
plt.show()

##### create a scatterplot of genome size and number of genes #####
# ensure the columns are numeric and remove missing values
g_data_clean = g_data[['genome_size', 'total_gene_count']].apply(pd.to_numeric, errors='coerce').dropna()

# create the scatterplot of genome size vs number of genes 
plt.figure(figsize=(8, 6))
plt.scatter(g_data_clean['genome_size'], g_data_clean['total_gene_count'], alpha=0.5, color='blue')
plt.title("Genome size vs Number of Genes (Subset)")
plt.xlabel("Genome Size")
plt.ylabel("Total Gene Count")

# fit a linear regression line through the origin 
X = g_data_clean['genome_size'].values.reshape(-1, 1)
Y = g_data_clean['total_gene_count'].values 

# perform linear regression through the origin (no constant term)
model = sm.OLS(Y, X)  # No constant term here, which implies regression through the origin
results = model.fit()

# plot the regression line 
plt.plot(g_data_clean['genome_size'], results.predict(), color='red', linewidth=2)
plt.show()

# compute the correlation (pearson correlation coefficient)
correlation, p_value = pearsonr(g_data_clean['genome_size'], g_data_clean['total_gene_count'])
print(f"Pearson Correlation: {correlation:.3f}, p-value: {p_value:.3e}")

# display the linear regression summary
print(results.summary())

##### create a boxplot of genome size for plant and fungi #####
# sort table by group
sorted_table = g_data.sort_values(by='group')

# collect plant genome size (rows 381746 to 381931 in the sorted table)
plant_genome_size = sorted_table.iloc[381748:381934]['genome_size']

# collect fungi genome size (rows 380685 to 381316 in the sorted table)
fungi_genome_size = sorted_table.iloc[380685:381317]['genome_size']

# create a boxplot
plt.figure(figsize=(10, 6))
plt.boxplot([plant_genome_size, fungi_genome_size], labels=['Plant', 'Fungi'])
plt.ylabel("Genome size (bp)")
plt.title("Plant and fungi genome size")
plt.show()

##### create a density plot of percentage difference between non-coding and protein-coding gene counts #####
# remove rows with missing values in protein_coding_gene_count and non_coding_gene_count columns
cleaned_df2 = g_data[['protein_coding_gene_count', 'non_coding_gene_count']].dropna()

# extract non-coding and protein-coding gene counts and convert to numeric, coercing errors to NaN
nc = pd.to_numeric(cleaned_df2['non_coding_gene_count'], errors='coerce')
pc = pd.to_numeric(cleaned_df2['protein_coding_gene_count'], errors='coerce')

# calculate delta_genes as the percentage difference
delta_genes = ((pc - nc) / nc) * 100

# remove NA values from delta_genes
delta_genes = delta_genes.dropna()

# check if there are any NAs left after cleaning
num_na_after_cleaning = delta_genes.isna().sum()
if num_na_after_cleaning > 0:
    print(f"Warning: {num_na_after_cleaning} NAs remain in delta_genes after cleaning.")

# check for infinite values
if np.isinf(delta_genes).any():
    print("Warning: Infinite values detected in delta_genes. Replacing with NaN.")
    delta_genes.replace([np.inf, -np.inf], np.nan, inplace=True)

# remove any remaining NaN values
delta_genes = delta_genes.dropna()

# check if there are any NAs left after all cleaning steps
num_na_after_all_cleaning = delta_genes.isna().sum()
print(f"Number of NAs in delta_genes after all cleaning: {num_na_after_all_cleaning}")

# proceed only if there are valid values left
if len(delta_genes) > 0:
    # calculate the density of delta_genes
    density = gaussian_kde(delta_genes)
    x_vals = np.linspace(0, 100000, 1000)
    density_vals = density(x_vals)

    # plot the density
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, density_vals, lw=1)
    plt.title("Protein-coding and non-coding gene percent difference")
    plt.xlabel("Percent difference")
    plt.ylabel("Density")
    plt.xlim(0, 100000)
    plt.grid()
    plt.show()
else:
    print("No valid values to plot.")

######################################################################################
######################################################################################
######################################################################################

##### part 02 #####
##### phylogenomic analysis #####

import os
import requests
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from scipy import stats

# set working directory
# NOTE: substitute path with your own
working_directory = "/Users/barba/Desktop/rggs-comparative_genomics_2_course/session_05"

# create the directory if it does not exist
os.makedirs(working_directory, exist_ok=True)

# change to the working directory
os.chdir(working_directory)

# download a FASTA MSA of 44 canid mt genomes
url = "https://raw.githubusercontent.com/josebarbamontoya/rggs_comparative_genomics_2/main/session_05/44canid_mt_genomes.fas"
response = requests.get(url)
response.raise_for_status()  # Raise an error for bad responses
with open("44canid_mt_genomes.fas", "wb") as file:
    file.write(response.content)

# tead FASTA MSA
alignment = AlignIO.read("44canid_mt_genomes.fas", "fasta")

# calculate distance matrix using the identity model (equivalent to JC69 for DNA)
calculator = DistanceCalculator('identity')
dist_matrix = calculator.get_distance(alignment)

# construct a tree using the neighbor-joining method
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dist_matrix)

# plot the tree
Phylo.draw(nj_tree, do_show=True)

# dave the Newick tree
Phylo.write(nj_tree, "nj_tree.nwk", "newick")

##### construct a maximum likelihood tree #####
# import necessary r packages
ape = importr('ape')
phangorn = importr('phangorn')

# set the path to the fasta file
fasta_file_path = os.path.join(working_directory, "44canid_mt_genomes.fas")

# read the fasta multiple sequence alignment
aln = robjects.r['read.FASTA'](fasta_file_path, type="DNA")

# create a phyDat object
phy_data = robjects.r['phyDat'](aln, type="DNA")

# read the neighbor-joining tree from a Newick file
nj_tree = robjects.r['read.tree']("nj_tree.nwk")

# fit the maximum likelihood model
fit = robjects.r['pml'](nj_tree, data=phy_data)

# optimize the tree
fit = robjects.r['optim.pml'](fit, model="GTR", optGamma=True)

# extract the optimized tree
ml_tree_r = fit.rx2('tree')

# save the maximum likelihood tree as a Newick file
robjects.r['write.tree'](ml_tree_r, file="ml_tree.nwk")

# read the Newick tree back into Biopython
ml_tree = Phylo.read("ml_tree.nwk", "newick")

# plot the maximum likelihood tree
Phylo.draw(ml_tree, do_show=True)

##### sort and compare trees #####
##### comparephylo and robinson-foulds distance #####

# activate automatic conversion of pandas and numpy objects to R objects
pandas2ri.activate()

# import necessary R packages
ape = importr('ape')

# read trees from Newick files
a = ape.read_tree("nj_tree.nwk")
b = ape.read_tree("ml_tree.nwk")

# convert outgroup list to R character vector
outgroup = r['c']('GF_GrayFox_Vermont')

# root trees with the specified outgroup
a = ape.root(a, outgroup=outgroup, resolve_root=True)
b = ape.root(b, outgroup=outgroup, resolve_root=True)

# ladderize the trees
a = ape.ladderize(a, right=False)
b = ape.ladderize(b, right=False)

# print summaries of the trees to inspect their structures using R's summary function
print("NJ tree summary:")
r['print'](r['summary'](a))

print("ML tre summary:")
r['print'](r['summary'](b))

# import necessary r packages
ape = importr('ape')
phangorn = importr('phangorn')
phytools = importr('phytools') 

# compare trees and plot
compare_results = ape.comparePhylo(a, b, plot=True, force_rooted=False, use_edge_length=True)

# print the entire comparison results to understand its structure
print(compare_results)

# combine the trees for distance calculation
trees = robjects.r['c'](a, b)

# calculate the multiRobinsonFoulds distance
rf_results = phytools.multiRF(trees)

# normalized rf distance calculation
n = 44 
plain_rf = rf_results[0]  # Get the plain RF distance
nrf = plain_rf / (2 * (n - 3))

# print the normalized rf distance
print(f"Normalized RF distance: {nrf[0]}")

##### Scatterplot of NJ and ML branch lengths #####
# Load trees (replace with actual paths)
try:
    nj_tree = Phylo.read("nj_tree.nwk", "newick")
    ml_tree = Phylo.read("ml_tree.nwk", "newick")
except Exception as e:
    print(f"Error loading trees: {e}")
    exit(1)

# Manually root trees with the specified outgroup
outgroup = "GF_GrayFox_Vermont"
nj_tree.root_with_outgroup(outgroup)
ml_tree.root_with_outgroup(outgroup)

# Manually ladderize the trees
def ladderize(tree):
    terminal_nodes = [clade for clade in tree.get_terminals()]
    terminal_nodes.sort(key=lambda x: x.name)
    for clade in tree.get_nonterminals():
        clade.clades.sort(key=lambda x: terminal_nodes.index(x) if x in terminal_nodes else float('inf'))
    return tree

nj_tree = ladderize(nj_tree)
ml_tree = ladderize(ml_tree)

# Extract branch lengths
def extract_branch_lengths(tree):
    lengths = []
    for clade in tree.get_nonterminals() + tree.get_terminals():
        if clade.branch_length is not None:  # Use branch_length instead of length
            lengths.append(clade.branch_length)
    return np.array(lengths)

nj_bl = extract_branch_lengths(nj_tree)
ml_bl = extract_branch_lengths(ml_tree)

# Check lengths of branch length arrays
print(f"NJ branch lengths count: {len(nj_bl)}")
print(f"ML branch lengths count: {len(ml_bl)}")

# if lengths are not equal, find the minimum length and truncate the longer one
if len(nj_bl) != len(ml_bl):
    min_length = min(len(nj_bl), len(ml_bl))
    nj_bl = nj_bl[:min_length]
    ml_bl = ml_bl[:min_length]
    print(f"Truncated lengths to: {min_length}")

# create a DataFrame for easy plotting
branch_lengths_df = pd.DataFrame({
    'NJ': nj_bl,
    'ML': ml_bl
})

# make scatterplot
plt.figure(figsize=(8, 8))
plt.scatter(branch_lengths_df['NJ'], branch_lengths_df['ML'], alpha=0.7)
plt.title("NJ and ML Tree Branch Length Comparison")
plt.xlabel("NJ Branch Length")
plt.ylabel("ML Branch Length")

# add diagonal line for reference
max_length = max(branch_lengths_df['NJ'].max(), branch_lengths_df['ML'].max())
plt.plot([0, max_length], [0, max_length], color='red', linestyle='--')  # y=x line

# add grid and show the plot
plt.grid(True)
plt.axis('equal')  # Equal aspect ratio for better comparison
plt.xlim(0, max_length)
plt.ylim(0, max_length)
plt.show()

