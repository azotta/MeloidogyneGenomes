#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
pd.set_option('display.max_columns', None)
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[8]:


import scipy.cluster.hierarchy as shc
from sklearn.metrics.pairwise import nan_euclidean_distances
import re


# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Data-loading-&amp;-pre-processing" data-toc-modified-id="Data-loading-&amp;-pre-processing-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Data loading &amp; pre-processing</a></span></li><li><span><a href="#Trees-computation" data-toc-modified-id="Trees-computation-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Trees computation</a></span></li><li><span><a href="#Trees-saving" data-toc-modified-id="Trees-saving-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Trees saving</a></span></li><li><span><a href="#Topology-extraction" data-toc-modified-id="Topology-extraction-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Topology extraction</a></span></li></ul></div>

# In[7]:


def create_square_matrix(block):
    """
    Switch from a long format to a square format. 
    
    Parameters :
    ------------
    block  (pandas dataframe)  -- a block dataframe (3 cols : ks, contig_X, contig_Y)
    
    Returns :
    ---------
    (pandas dataframe)  -- square version of the block
    """
    
    ##  Make a list of unique contig names from the concatenation of
    #   contig_X and contig_Y
    contigs = list(set(list(block.contig_X)+list(block.contig_Y)))
    lst_rows = []
    ##  Pairwize combination of contigs names
    for c1 in contigs:
        for c2 in contigs:
            tmp_dict = {'contig_X':c1, 'contig_Y':c2}
            #   case when we are considering the same contig twice (assigned 
            #   value : 0)
            if c1 == c2:
                tmp_dict['ks'] = 0
            else:
                #   check if the paiwize combination exist
                pw = block.loc[(block.contig_X==c1)&(block.contig_Y==c2)]
                #   if it exist, save the value
                if len(pw) > 0:
                    tmp_dict['ks'] = pw.ks.values[0]
                #   else, assign nan.
                else:
                    tmp_dict['ks'] = np.nan
            lst_rows.append(tmp_dict)
    
    return pd.DataFrame(lst_rows).pivot(index='contig_X', columns='contig_Y', values='ks')

def dendrogram_to_newick(dendrogram, labels):
    """
    Transform a clustering dendogram to a newick format with its leaf labels 
    and branch size.
    
    Parameters :
    ------------
    dendogram   (object)    --  the dendogram
    labels  (list)  --  list of dendogram labels (in the same order as the 
                        dendogram)
    
    Returns :
    ---------
    (str)   --  newick formated string
    """
    root_node = sch.to_tree(dendrogram)
    return get_newick(root_node, labels)

def get_newick(node, labels):
    """
    Extract infos from the tree
    """
    if node.is_leaf():
        return f"{labels[node.id]}:{node.dist}"
    else:
        left_newick = get_newick(node.left, labels)
        right_newick = get_newick(node.right, labels)
        return f"({left_newick},{right_newick}):{node.dist}"


# # Data loading & pre-processing

# In[5]:


##  Load the data
df_data = pd.read_csv('Mjav_block_def.tsv', sep='\t')

##  Transform the contig columns in two separated columns names 'contig_X' 
#   and 'contig_Y'
df_data[['contig_X', 'contig_Y']]=df_data.contigs.str.split('-', expand=True)

##  Delete the contig columns which is useless now
df_data.drop(columns=['contigs'], inplace=True)

df_data.head(10)


# # Trees computation

# We can now parse the file block by block and create a tree for each block.

# In[12]:


import scipy.cluster.hierarchy as shc
from sklearn.metrics.pairwise import nan_euclidean_distances
import re
##  Init a dict where the newick trees will be saved for each block
dct_newicks = {}
##  For each block in the dataframe
for blk_name in df_data.block_name.unique():
    ##  Extract the block related info 
    block = df_data[df_data.block_name==blk_name].drop(columns=['block_name'])
    ##  Create a square version of the block dataframe.
    block_sq = create_square_matrix(block)
    ##  Compute an euclidian distance matrix from the block info
    block_dist_matrix = nan_euclidean_distances(block_sq.values, block_sq.values)
    ##  Perform the hierarchical/agglomerative clustering using the ward method.
    block_link_matrix = shc.linkage(block_dist_matrix, method='ward')
    ##  Extract the tree structure as newick tree and save it
    block_newick = dendrogram_to_newick(block_link_matrix, block_sq.index)
    dct_newicks[blk_name] = block_newick


# # Trees saving

# We'll now save each one of the trees as an individual newick file. **Careful ! If file exists, they might be overwritten !** 

# In[13]:


##  Create an output directory 
os.makedirs('newick_trees', mode = 0o755, exist_ok = True)

##  Save the newick as individual files 
for block, tree in dct_newicks.items():
    with open(os.path.join('newick_trees', f'{block}.nwk'), 'w') as f:
        f.write(tree)


# # Topology extraction

# In[14]:


##  Extract the topology for each tree and store it in a dataframe
lst_dct_topo = []
for block, tree in dct_newicks.items():
    topo = re.sub('[0-9:.]','', tree)
    lst_dct_topo.append({'block_name':block, 'topology':topo})

##  Create a dataframe from it
df_topo = pd.DataFrame(lst_dct_topo)

##  Save the topologies
df_topo.to_csv(os.path.join('newick_trees', 'topologie.tsv'), 
               sep='\t', 
               index=False)


# In[ ]:




