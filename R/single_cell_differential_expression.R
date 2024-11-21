### single cell differential expression testing


library(edgeR)
library(stringr)
library(dplyr)
library(pbapply)

#### load in bulk data


bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv')

bulk_cells <- str_split_fixed(colnames(bulk_data), 'r', 2)[,1] |> unique()



### load in single cell data


sc_bioreps <- readRDS('Data/CeNGEN_biorep_aggregate_counts_data_240709.RDS')

prop2count <- read.table('Data/sc_biorep_prop2count_050124.tsv.gz')
prop2count <- prop2count[order(rownames(prop2count)),]

sc_counts <- sc_bioreps$counts[,sc_bioreps$nCells >= 10]
sc_counts <- sc_counts[order(rownames(sc_counts)), ]


colnames(sc_counts) <- sapply(colnames(sc_counts), function(x){gsub('-', '\\.', x)})

common_reps <- intersect(colnames(sc_counts), colnames(prop2count))

sc_counts <- sc_counts[,common_reps]
prop2count <- prop2count[,common_reps]

waldo::compare(colnames(sc_counts), colnames(prop2count))


replicate_groups <- str_split_fixed(common_reps, '__', 2)[,1]
cell_types <- unique(replicate_groups)

common.cells <- intersect(bulk_cells, cell_types)
common.cells <- c(common.cells, 'VD_DD')

## verify > 1 replicate per cell type
table(replicate_groups)

keep_genes <- rowSums(sc_counts > 2) > 1
keep_replicates <- replicate_groups %in% common.cells

edgeR_counts <- DGEList(counts = sc_counts[keep_genes, keep_replicates],
                        group = replicate_groups[keep_replicates])
edgeR_counts <- calcNormFactors(edgeR_counts, method = 'TMM')

group <- factor(cell_types)
design_no_int <- model.matrix(~0+group, data=edgeR_counts$samples)
colnames(design_no_int) <- levels(edgeR_counts$samples$group)


edgeR_counts <- estimateDisp(edgeR_counts, design_no_int)
fit <- glmQLFit(edgeR_counts, design_no_int, robust=TRUE)

saveRDS(fit, 'data_out/single_cell_counts_edgeR_fit.rds')

fit <- readRDS('single_cell_counts_edgeR_fit.rds')

### make contrasts ----

contrasts_names <- factor()
contrast_cells <- fit$samples$group |> unique() |> sort()
cells_shrink <- contrast_cells
for(f in contrast_cells){
  cells_shrink <- cells_shrink[(cells_shrink != f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names
contrasts_names
length(contrasts_names)

contrasts_edgeR <- lapply(contrasts_names, function(v){ cont <-
  makeContrasts(contrasts = contrasts_names[[v]], levels = fit$design)
})



qlf_tables <- pblapply(contrasts_edgeR, function(v){
  #print(v)
  qlf <- glmQLFTest(glmfit = fit, contrast = v)
  
  qlf <- topTags(qlf, n = nrow(qlf), adjust.method =  'BH')[[1]]
  
  qlf <- qlf[order(rownames(qlf)),]
  
  return(qlf)
})


saveRDS(qlf_tables, 'data_out/sc_counts_qlf_tables.rds')



### run for prop2count


edgeR_prop2count <- DGEList(counts = prop2count[keep_genes, keep_replicates],
                            group = replicate_groups[keep_replicates])
edgeR_prop2count <- calcNormFactors(edgeR_prop2count, method = 'TMM')

group <- factor(cell_types)
design_no_int <- model.matrix(~0+group, data=edgeR_prop2count$samples)
colnames(design_no_int) <- levels(edgeR_prop2count$samples$group)


edgeR_prop2count <- estimateDisp(edgeR_prop2count, design_no_int)
fit <- glmQLFit(edgeR_prop2count, design_no_int, robust=TRUE)

saveRDS(fit, 'data_out/sc_prop_edgeR_fit.rds')

fit <- readRDS('data_out/sc_prop_edgeR_fit.rds')

### make contrasts ----

contrasts_names <- factor()
contrast_cells <- fit$samples$group |> unique() |> sort()
cells_shrink <- contrast_cells
for(f in contrast_cells){
  cells_shrink <- cells_shrink[(cells_shrink != f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names
contrasts_names
length(contrasts_names)

contrasts_edgeR <- lapply(contrasts_names, function(v){ cont <-
  makeContrasts(contrasts = contrasts_names[[v]], levels = fit$design)
})


qlf_tables <- pblapply(contrasts_edgeR, function(v){
  #print(v)
  qlf <- glmQLFTest(glmfit = fit, contrast = v)
  
  qlf <- topTags(qlf, n = nrow(qlf), adjust.method =  'BH')[[1]]
  
  qlf <- qlf[order(rownames(qlf)),]
  
  return(qlf)
})


saveRDS(qlf_tables, 'data_out/sc_prop_qlf_tables.rds')



