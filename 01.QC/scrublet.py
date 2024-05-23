##############################################################################
# Script information                                                      
# Title: Doublet removing
# Author: Erping Long
# Date: 2021-03-15
# Description: None
##############################################################################

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

input_dir = '/filtered_feature_bc_matrix/'
counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))


scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

import pandas as pd

df = pd.DataFrame(predicted_doublets)
df.to_csv('/scrublet.csv')







