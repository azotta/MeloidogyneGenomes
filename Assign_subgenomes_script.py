#!/usr/bin/env python

## Authors: Ana Zotta and Georgios D Koutsovoulos

import pandas as pd
import sys
import re
import sys
import os
import glob
from __future__ import print_function

## The input files are from the McScanX analysis
## 

# Open gff file
gff_file_path = 'species.gff'

# Define the column names for the GFF file
column_gff = ['contig', 'gene', 'start', 'end']

# Open the ks file
col_file = 'species_values_kaks.txt'

# Define the column names for the collinearity file
column_col = ['gene1', 'gene2', 'ks']

# Read the GFF file into a pandas DataFrame
df_gff = pd.read_csv(gff_file_path, sep='\t', comment='#', names=column_gff)

# Read the collinearity file
df_col = pd.read_csv(col_file, sep='\t', comment='#', names=column_col) 


file_groups = [] # open a list file to be incremented on the for loop
file_list = []
filtered_groups = []
#For each gene in the gff which is present more than 3 time in the collinearity file, print all the match in the same line 
for gene_name in df_gff['gene'].unique():
    # Count occurrences of the gene in the other DataFrame
    
    count = df_col[(df_col['gene1'] == gene_name) | (df_col['gene2'] == gene_name)].shape[0]
    if count >= 3: 
        matching_genes = df_col[df_col['gene1'] == gene_name]['gene2'].tolist()
        result_string = " ".join([str(item)for item in matching_genes])
        tmp_genes = [gene_name]
        tmp_list = [gene_name]+matching_genes # Create a list of a list
        file_list = file_list+[tmp_genes] # Increment the list of genes
        file_groups = file_groups+[tmp_list] # Increment the list file each loop        
        filtered_groups = [row for row in file_groups if len(row) == 4]  ## Filter accordling to the ploidy of the species being analyzed (4 for tetraploid species/3 for triploid)

df = pd.DataFrame(file_groups)  # Put the list file into a DataFrame      
df1 = pd.DataFrame(file_list)
df3 = pd.DataFrame(filtered_groups)

len(df3. index)

len(df. index)


##Create file with all genes from each contig belonging to each block 
## File to be used later on the loops to retrieve gene1 gene2 ks_values
df3.to_csv('result_genes_in_groups.csv', sep='\t', index=False, header=False) 


# Save a list of genes in blocks

df1.to_csv('list_uniq_genes_within_blocks.csv', sep='\t', index=False, header=True) 


## This part will use the file with all the ks values and the second file with the name of all genes to produce all the pairs
## the resulting file "ks_pairs.txt" will be used in the next analyses, to create one file per groups with the ks values, separated by directory

D={}

#
colinearity_size = 4 ##change according to ploidy

with open('species_values_kaks.txt', 'r') as lines_handle:
    for line in lines_handle:
        L = line.rstrip().split("\t")
        gene1 = L[0]
        gene2 = L[1]
        ks = L[2]
        if gene1 not in D:
            D[gene1]={}
        if gene2 not in D:
            D[gene2]={}
        D[gene1][gene2] = ks
        D[gene2][gene1] = ks

ks_pairs = open("ks_pairs.txt","w")
with open('result_genes_in_groups.csv', 'r') as lines_handle:
    for line in lines_handle:
        L = line.rstrip().split("\t")
        if len(L)==colinearity_size:
            i = 0
            while i < colinearity_size:
                k = i + 1
                while k < colinearity_size:
                    if L[k] not in D[L[i]]:
                        ks_pairs.write(L[i] + "\t" + L[k] + "\t" + "NA" + "\n")
                    elif D[L[i]][L[k]] == "-2":  ## To eliminate values -2 in the ks columns
                        ks_pairs.write(L[i] + "\t" + L[k] + "\t" + "NA" + "\n")
                    elif float(D[L[i]][L[k]]) > 1:  ## To eliminate values of ks bigger than 1
                        ks_pairs.write(L[i] + "\t" + L[k] + "\t" + "NA" + "\n")    
                    else:
                        ks_pairs.write(L[i] + "\t" + L[k] + "\t" + D[L[i]][L[k]] + "\n")
                    k += 1
                i += 1

ks_pairs.close()

##This part will produce a list of groups and their correspondence in the file "ks". 



filename = 'result_genes_in_groups.csv' ## The file is the one produced by the first step (df3)

blocks = {}
with open(filename) as file:
    i = 0
    for line in file:
        i += 1
        genes = line.strip().split()
        contigs = []
        for gene in genes:
            parts = gene.split('g', 2)  # the 2 means that we cut twice to search for the second 'g'
            if len(parts) > 1:
                contigs.append(parts[0] + 'g' + parts[1])
        block = "\t".join(contigs)
        if block not in blocks:
            blocks[block] = []
        blocks[block].append(str(i))

for block, fasta_ids in blocks.items():
    fasta_ids_str = "\t".join(fasta_ids)
    #print(f"{block}\t{fasta_ids_str}")
    
sorted_blocks = sorted(blocks.keys(), key=lambda block: len(block.split("\t")), reverse=True)

blocks_file = open("blocks_from_code.txt","w")
for block in sorted_blocks:
    fasta_ids_str = "\t".join(blocks[block])
    #print(f"{block}")  # If I just want to check the number of the groups
    blocks_file.write(f"{block}\t{fasta_ids_str}\n")   ## If I want the line number from the file result_teste.csv, however, this is not the same order as the file ks. So the association is not good. 
blocks_file.close()    

### so to sort uniq this file, I used the script in awk "verify_repeated_columns.txt"

#This code will generate one directory for each block from the file produced above "blocks_from_code.txt"



res_filename = 'ks_pairs.txt'
filename = 'blocks_from_code.txt'

with open(res_filename, 'r') as handle:
    divs = [line.strip() for line in handle]

with open(filename, 'r') as infile:
    count = 0
    for line in infile:
        count += 1
        p = line.strip().split()
        for i in range(4):
            p[i] = p[i].replace('g\d+', '').replace('Mjav_v4_contig_', '')

        num = len(p) - 4
        dir_name = f"BLOCK_{count}_contigs_{p[0]}_{p[1]}_{p[2]}_{p[3]}_num_{num}"
        os.mkdir(dir_name)

        with open(f"{dir_name}/divergence.txt", 'w') as outfile:
            for i in range(4, len(p)):
                loc = (int(p[i]) - 1) * 6  #need to be 6 because we have 6 connections for each gene pairs (need to change to three if its a triploid species)
                for k in range(6):
                    outfile.write(f"{divs[loc + k]}\n")


# For each directory (corresponding to each synteny block), disconsider the name of genes, and print only the number of contigs, and the ks values




list_of_directories = glob.glob('/Users/azotta/DocumentsLocal/JupyterNotebooks/MeloidogyneGenomes/SeparateGenomes/Mjav/BLOCK_*/')

for directory in list_of_directories:
    output_file_path = directory + 'result.txt'

    with open(output_file_path, 'w') as output_file:
        list_of_files = sorted(glob.glob(directory + 'divergence.txt'))

        for f in list_of_files:
            with open(f, 'r') as file:
                lines = file.readlines()

            for line in lines:
                # Split the line based on whitespace
                columns = line.split()
                contig_1 = columns[0]
                contig_2 = columns[1]
                divergence = columns[2]

                # Extract the characters before and after the first "g" from the contig names and remove underscores
                contig_1_name = contig_1.split('g', 1)[-1].split('g', 1)[0].replace('_', '')
                contig_2_name = contig_2.split('g', 1)[-1].split('g', 1)[0].replace('_', '')

                result_string = f"{contig_1_name}-{contig_2_name} {divergence}\n"
                output_file.write(result_string)


# Create a file with all the directories and results.txt files


get_ipython().system('ls -d BLOCK_*/result.txt > list_files.txt')


# For each directory, create a file with the ks values for each pair of contigs, and calculate the median for each contig pair
# For triploid species, the final final will be composed of three lines, for tetraploid species, 6 lines




# Set working directory
os.chdir("/Users/azotta/DocumentsLocal/JupyterNotebooks/MeloidogyneGenomes/SeparateGenomes/Mjav/")

# Read list of files
filelist = pd.read_table("list_files.txt", header=None)

for filename in filelist[0]:
    # Read data from file with custom delimiter "-"
    x = pd.read_table(filename, header=None, names=["V1","V2"], sep=" ")
    
    # Check that data is formatted correctly
    if not all(x.columns == ["V1", "V2"]):
        raise ValueError(f"File '{filename}' is not formatted correctly.")
    
    # Create plot
    plotname = filename.replace(".txt", "_cov.pdf")
    plt.figure()
    plt.title(filename)
    for group, data in x.dropna().groupby("V1"):
        plt.hist(data["V2"], alpha=0.7, label=group)
    plt.legend()
    plt.savefig(plotname)
    plt.close()
    
    
    # Calculate median coverage and save results to CSV file
    medians = x.dropna().groupby("V1")["V2"].median().reset_index()
    means_file = filename.replace(".txt", "_medians.csv")
    medians.to_csv(means_file, index=False)





# Select blocks where only the right number of lines is present (three or six, for example), this threshold is important to select only "perfect" synteny relations, to better assign the subgenomes




list_of_files = sorted(glob.glob('/Users/azotta/DocumentsLocal/JupyterNotebooks/MeloidogyneGenomes/SeparateGenomes/Mjav/BLOCK_*/result_medians.csv'))

threshold = 6

result = []
for f in list_of_files:
    df_tmp = pd.read_csv(f)
    #print(df_tmp)
    num_rows = len(df_tmp)
    if num_rows == threshold: 
        block_name = re.search(".*\/(BLOCK_\d*)_.*",f).group(1)
        df_tmp = df_tmp.sort_values(by='V2', ascending=False).reset_index(drop=True).assign(block_name=block_name)
        #print(df_tmp)
        second_lowest = df_tmp.loc[4, 'V1'].split('-')
        #print (second_lowest)
        lowest = df_tmp.loc[5, 'V1'].split('-')
        #print (lowest)
        result.append(df_tmp)
        
        
df_result = pd.concat(result).rename(columns={'V1':'contigs', 'V2':'ks'})



# Save the file


df_result.to_csv('Mjav_block_def.tsv', sep='\t', index=False)


# Use this file as input for genome structure 




