
## load libraries
library(MIND)
library(data.table)
library(pbapply)
library(dplyr)
library(tidyverse)


## load in ground truth genes

all_genes_bulk <- read.csv('references/bulk_all_ground_truth_121023.csv', 
                           row.names = 1)

UNN_ground_truth_mtx <- read.table('references/ubituiqtous_and_nonNeuronal_gt_genes_matrix_042222.tsv', sep = '\t')

all_gt_genes <- c(rownames(UNN_ground_truth_mtx), rownames(all_genes_bulk)) |> sort()


## load in bulk data

bulk <- read.table('Data/Bulk_data_bsn12_231211.tsv', sep = '\t')


## load in proportions estimations

proportions <- read.table('Data/NNLS_average_across_100_bootstraps.log1p.30Cells.231211.tsv', sep = '\t')

## load pseudobulk data

seurat_pseudobulk_tissue <- read.table('Data/seurat_pseudobulk_tissue_072922.tsv', sep = '\t')


sc_meta = data.frame(row.names = colnames(seurat_pseudobulk_tissue),
                     sample = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,2],
                     cell_type = str_split_fixed(colnames(seurat_pseudobulk_tissue), '__', 2)[,1])


counts_prior <- get_prior(sc = seurat_pseudobulk_tissue, meta_sc = sc_meta, filter_pd = F)

bMIND.common.genes <- intersect(all_gt_genes, rownames(counts_prior$profile))



bMIND_genes <- list()
counter = 0

for(gene in bMIND.common.genes){
  
  bMIND_genes[[gene]] <- 
    tryCatch(
      bMIND(bulk = log2(1+bulk[gene,]),
            frac = t(proportions[colnames(counts_prior$profile),]),
            profile = counts_prior$profile[gene,]),
      error=function(e) NULL)
  counter <- counter + 1
  print(counter)
  
}

bMIND_genes[[1]]$A['Neuron',]
saveRDS(bMIND_genes, 'bMIND_bsn12_deconvolution_list_using_prior_only_version1_031024.RDS')

bMIND_genes |> names()

bMIND_neuron <- sapply(bMIND_genes |> names(), function(gene){
  
  bMIND_genes[[gene]]$A['Neuron',]
  
}) |> t()


