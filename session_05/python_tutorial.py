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

# see table content
#g_data 

# remove the '#' symbol from each line
cleaned_lines = [line.replace("#", "") for line in g_data]

# convert cleaned lines to a DataFrame, skipping the first line
# use stringio to read the cleaned lines as a cvs-like structure
data_string = "".join(cleaned_lines)
g_data = pd.read_csv(StringIO(data_string), sep="\t", skiprows=1, quoting=3, engine='python')

# export the edited DataFrame to a new file
g_data.to_csv("assembly_summary_refseq_cleaned.txt", sep="\t", index=False, quoting=3)

# read the cleaned dataframe with low_memory set to false
g_data = pd.read_csv("assembly_summary_refseq_cleaned.txt", sep="\t", low_memory=False)

# count number of rows and columns
# number of rows
num_rows = g_data.shape[0]
# number of columns
num_cols = g_data.shape[1]
# print numbers
print(f"Number of rows: {num_rows}")
print(f"Number of columns: {num_cols}")

# count number of genomes for each major lineage
genome_counts = g_data['group'].value_counts()
print("\nGenome counts for each major lineage:")
print(genome_counts)

##### summary statistics of gc percent #####
print("\nSummary statistics for GC percent:")
gc_percent_summary = g_data['gc_percent'].describe()  # Summary statistics
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

# compute maximum nuber of base pairs
g_data['genome_size'].max()

# compute average number of base pairs
g_data['genome_size'].mean()

##### create a bar plot of genome size #####
#plt.figure(figsize=(25, 15))
#plt.bar(g_data['organism_name'], g_data['genome_size'], color='green')
#plt.title('Genome Size Barplot')
#plt.xlabel('Organism')
#plt.ylabel('Base Pairs')
#plt.ylim(0, 100000000)  # Set y-axis limits from 0 to 1 billion (adjust as necessary)
#plt.xticks(rotation=45)
#plt.grid(axis='y', alpha=0.75)
#plt.tight_layout()
#plt.show()

# create a bar plot of genome size
plt.figure(figsize=(25, 15))
plt.bar(range(len(g_data)), g_data['genome_size'], color='green')  # Use indices for the x-axis
plt.title('Genome Size Barplot')
plt.ylabel('Base Pairs')
plt.ylim(0, 100000000)
plt.box(on=None)
plt.xticks([])
plt.grid(axis='y', alpha=0.75)
plt.tight_layout()
plt.show()

##### create a stacked barplot of coding and noncoding genes #####
# ensure the columns are numeric, replacing non-numeric entries with NaN and dropping rows with missing data
g_data_clean = g_data[['protein_coding_gene_count', 'non_coding_gene_count']].apply(pd.to_numeric, errors='coerce').dropna()

# create the data matrix (stacked values for each organism)
data_matrix = np.array([g_data_clean['protein_coding_gene_count'], g_data_clean['non_coding_gene_count']])

# create the stacked bar plot
# define the figure size and the positions of the bars
fig, ax = plt.subplots(figsize=(10, 6))
organisms = np.arange(g_data_clean.shape[0])

# plot the stacked bars (first the non-coding genes, then add protein-coding genes on top)
bar1 = ax.bar(organisms, data_matrix[1], color='red', label='Non-coding genes')  # Non-coding genes
bar2 = ax.bar(organisms, data_matrix[0], bottom=data_matrix[1], color='blue', label='Protein-coding genes')
# add titles and labels
ax.set_title("Number of Protein-coding and Non-coding Genes per Organism")
#ax.set_xlabel("Organism")
ax.set_ylabel("Gene Count")
plt.ylim(0, 132906)
# add legend
ax.legend()
# edit the x-axis to show the organism labels if available (or else show indices)
#ax.set_xticks(organisms)
#ax.set_xticklabels(g_data_clean.index, rotation=90)
# Show the plot
plt.tight_layout()
plt.show()

##### create the stacked barplot for the subset #####
# ensure the columns are numeric, replacing non-numeric entries with NaN and dropping rows with missing data
g_data_clean = g_data[['protein_coding_gene_count', 'non_coding_gene_count']].apply(pd.to_numeric, errors='coerce').dropna()

# select rows 10 to 100 (python uses zero-based indexing, so this is equivalent to rows 11 to 100 in r)
subset = g_data_clean.iloc[10:100]
# create the data matrix (stacked values for each organism in the subset)
data_matrix = np.array([subset['protein_coding_gene_count'], subset['non_coding_gene_count']])
# create the stacked bar plot
fig, ax = plt.subplots(figsize=(10, 6))
organisms = np.arange(subset.shape[0])
# plot the stacked bars (non-coding genes first, then protein-coding genes stacked on top)
bar1 = ax.bar(organisms, data_matrix[1], color='red', label='Non-coding genes')
bar2 = ax.bar(organisms, data_matrix[0], bottom=data_matrix[1], color='blue', label='Protein-coding genes')
# add titles and labels
ax.set_title("Number of Protein-coding and Non-coding Genes (Rows 10 to 100)")
ax.set_ylabel("Gene Count")
# set dynamic y-limit based on data
#ax.set_ylim(0, data_matrix.sum(axis=0).max() * 1.1)
# set hard y-limit
plt.ylim(0, 132906)
# add legend
ax.legend()
# uncomment the following to add organism labels if necessary (index can be used here)
# ax.set_xticks(organisms)
# ax.set_xticklabels(subset.index, rotation=90)
# Sshow the plot
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

# compute the correlation (Pearson correlation coefficient)
correlation, p_value = pearsonr(g_data_clean['genome_size'], g_data_clean['total_gene_count'])
print(f"Pearson Correlation: {correlation:.3f}, p-value: {p_value:.3e}")

# display the linear regression summary
print(results.summary())

##### create a boxplot or a vioplot of genome size for plant and fungi #####
# sort table by group
sorted_table = g_data.sort_values(by='group')

# collect plant genome size (rows 381746 to 381931 in the sorted table)
plant_genome_size = sorted_table.iloc[381746:381931]['genome_size']

# collect fungi genome size (rows 380685 to 381316 in the sorted table)
fungi_genome_size = sorted_table.iloc[380685:381316]['genome_size']

# create a boxplot
plt.figure(figsize=(10, 6))
plt.boxplot([plant_genome_size, fungi_genome_size], labels=['Plant', 'Fungi'])
plt.ylabel("Genome Size (bp)")
plt.title("Plant and Fungi Genome Size - Boxplot")
plt.show()

# alternatively, create a violin plot
#plt.figure(figsize=(10, 6))
#sns.violinplot(data=[plant_genome_size, fungi_genome_size])
#plt.xticks([0, 1], ['Plant', 'Fungi'])
#plt.ylabel("Genome Size (bp)")
#plt.title("Plant and Fungi Genome Size - Violin Plot")
#plt.show()

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
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri
import rpy2.robjects as robjects
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio import Phylo

# set working directory
# NOTE: substitute path with your own
working_directory = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05"
os.chdir(working_directory)

# download a fasta msa of 44 canid mt genomes
url = "https://raw.githubusercontent.com/josebarbamontoya/rggs_comparative_genomics_2/main/session_05/44canid_mt_genomes.fas"
response = requests.get(url)
with open("44canid_mt_genomes.fas", "wb") as file:
    file.write(response.content)

# read FASTA MSA
alignment = AlignIO.read("44canid_mt_genomes.fas", "fasta")

# calculate distance matrix using the JC69 model
calculator = DistanceCalculator('identity')  # 'identity' corresponds to JC69 for DNA sequences
dist_matrix = calculator.get_distance(alignment)

# construct a tree using the neighbor-joining method
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dist_matrix)

# plot the tree
Phylo.draw(nj_tree, do_show=True)

# save the Newick tree
Phylo.write(nj_tree, "nj_tree.nwk", "newick")

##### construct a maximum likelihood tree #####
# import necessary r packages
ape = importr('ape')
phangorn = importr('phangorn')

# set the path to the fasta file
fasta_file_path = "/Users/barba/Desktop/rggs-copmarative_genomics_2_course/session_05/44canid_mt_genomes.fas"

# read the fasta multiple sequence alignment
aln = robjects.r['read.FASTA'](fasta_file_path, type="DNA")

# create a phydat object
phy_data = robjects.r['phyDat'](aln, type="DNA")

# create a neighbor-joining tree (replace 'nj_tree' with your actual tree if needed)
nj_tree = robjects.r['nj'](robjects.r['dist.dna'](aln, model="raw"))  # Assuming nj_tree needs to be created from the alignment

# fit the maximum likelihood model
fit = robjects.r['pml'](nj_tree, data=phy_data)

# optimize the tree
fit = robjects.r['optim.pml'](fit, model="GTR", optGamma=True)

# plot the ml tree
ml_tree = fit.rx2('tree')
robjects.r['plot'](ml_tree, main="Maximum Likelihood Tree")
robjects.r['nodelabels']()
robjects.r['edgelabels']()

# save the newick tree
robjects.r['write.tree'](ml_tree, file="ml_tree.nwk")

##### sort and compare trees #####
##### comparephylo and robinson-foulds distance #####

# activate automatic conversion of pandas and numpy objects to r objects
pandas2ri.activate()

# import necessary R packages
ape = importr('ape')

# read trees from Newick files
a = ape.read_tree("nj_tree.nwk")
b = ape.read_tree("ml_tree.nwk")

# convert outgroup list to r character vector
outgroup = robjects.StrVector(['GF_GrayFox_Vermont'])

# root trees with the specified outgroup
a = ape.root(a, outgroup=outgroup, resolve_root=True)
b = ape.root(b, outgroup=outgroup, resolve_root=True)

# ladderize the trees
a = ape.ladderize(a, right=False)
b = ape.ladderize(b, right=False)

# compare trees and plot
compare_results = ape.comparePhylo(a, b, plot=True, force_rooted=True, use_edge_length=True)

# print the entire comparison results to understand its structure
print(compare_results)

# combine the trees for distance calculation
trees = robjects.r('c(a, b)')

# calculate the multiRobinsonFoulds distance
rf_results = robjects.r('multiRF(trees)')

# norm RF distance calculation
n = 44  # Number of taxa
plain_rf = rf_results[0]  # Get the plain RF distance
nrf = plain_rf / (2 * (n - 3))

# print the normalized RF distance
print(f"Normalized RF distance: {nrf[0]}")

##### scatterplot of nj and ml branch lengths #####
# load trees (replace with actual paths)
try:
    nj_tree = Phylo.read("nj_tree.nwk", "newick")
    ml_tree = Phylo.read("ml_tree.nwk", "newick")
except Exception as e:
    print(f"Error loading trees: {e}")
    exit(1)

# manually root trees with the specified outgroup
outgroup = "GF_GrayFox_Vermont"
nj_tree.root_with_outgroup(outgroup)
ml_tree.root_with_outgroup(outgroup)

# manually ladderize the trees
def ladderize(tree):
    terminal_nodes = [clade for clade in tree.get_terminals()]
    terminal_nodes.sort(key=lambda x: x.name)
    for clade in tree.get_nonterminals():
        clade.clades.sort(key=lambda x: terminal_nodes.index(x) if x in terminal_nodes else float('inf'))
    return tree

nj_tree = ladderize(nj_tree)
ml_tree = ladderize(ml_tree)


# extract branch lengths
def extract_branch_lengths(tree):
    lengths = []
    for clade in tree.get_nonterminals() + tree.get_terminals():
        if clade.branch_length is not None:  # Use branch_length instead of length
            lengths.append(clade.branch_length)
    return np.array(lengths)

nj_bl = extract_branch_lengths(nj_tree)
ml_bl = extract_branch_lengths(ml_tree)

# create a dataframe for easy plotting
branch_lengths_df = pd.DataFrame({
    'NJ': nj_bl,
    'ML': ml_bl
})

# make scatterplot
plt.figure(figsize=(8, 8))
plt.scatter(branch_lengths_df['NJ'], branch_lengths_df['ML'], alpha=0.7)
plt.title("NJ and ML tree branch length comparison")
plt.xlabel("NJ Branch Length")
plt.ylabel("ML Branch Length")

# fit linear regression through the origin
slope, intercept, r_value, p_value, std_err = stats.linregress(branch_lengths_df['NJ'], branch_lengths_df['ML'])
plt.plot(branch_lengths_df['NJ'], slope * branch_lengths_df['NJ'], color='red', label='Fit: y = {:.2f}x'.format(slope))
plt.legend()

# show plot
plt.grid()
plt.show()

# compute correlation
correlation = branch_lengths_df.corr(method='pearson').iloc[0, 1]
print(f"Pearson correlation coefficient: {correlation:.4f}")

# linear regression summary
print(f"Linear regression through the origin: y = {slope:.4f}x")
print(f"R-squared: {r_value**2:.4f}, p-value: {p_value:.4f}, standard error: {std_err:.4f}")

