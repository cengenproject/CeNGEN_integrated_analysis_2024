### expression of protein coding genes not detected in any single cell clusters


### datasets needed, and assumed already loaded: ----

# CeNGEN_TPM: bulk annotation aggregates for single cell data
# aggr_raw_TMM: aggregated raw bulk GeTMM values, arithmetic mean per cell type
# GeTMM_cor_s_contaminant: genes x contaminants spearman correlations
# ws289: wormbase gene metadata, from wbData package
# nn_genes: a list of non-neuronal genes, curated from wormbase using simpleMine



get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}


### load data, normalize ----

bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv')


bulk_meta <- read.table('Data/bulk_bsn12_metadata.tsv', sep = '\t')

bulk_meta <- bulk_meta[rownames(bulk_data),]

dim(bulk_data)

dim(bulk_meta)



bulk_raw_TMM <- DGEList(bulk_data, group = str_split_fixed(colnames(bulk_data), 'r', 2)[,1])
bulk_raw_TMM <- calcNormFactors(bulk_raw_TMM)
bulk_raw_TMM <- cpm(bulk_raw_TMM, normalized.lib.sizes = T)

dim(bulk_raw_TMM)


### get average profile per cell type ----
dim(bulk_raw_TMM)
aggr_raw_TMM <- bulk_raw_TMM
colnames(aggr_raw_TMM) <-str_split_fixed(colnames(aggr_raw_TMM),"r",2)[,1]
aggr_raw_TMM <- data.frame(vapply(unique(colnames(aggr_raw_TMM)), function(x)
  rowMeans(aggr_raw_TMM[,colnames(aggr_raw_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_TMM)) ))


### load Alexis single cell thresholds
ws289 <- data.frame(wbData::wb_load_gene_ids(289))
rownames(ws289) <- ws289$gene_id
ws289_prc <- ws289[ws289$biotype=='protein_coding_gene',]

sc_gene_thresholds <- readRDS('references/211028_genes_categorized_by_pattern.rds')

nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, rownames(ws289_prc))


CeNGEN_TPM <- cengenDataSC::cengen_TPM_bulk


not_detected_Seurat <- setdiff(rownames(aggr_raw_TMM), rownames(CeNGEN_TPM))
not_detected_Seurat <- intersect(not_detected_Seurat, rownames(ws289_prc))
sum(nondetected_threshold_genes %in% not_detected_Seurat)
length(not_detected_Seurat)

nondetected_combined <- c(not_detected_Seurat, nondetected_threshold_genes) |> unique()

### subset to the genes of interest given either criteria and combine the matrices ----
#aggr_not_detected_Seurat <- aggr_raw_TMM[not_detected_Seurat,]
aggr_nondetected_threshold_genes <- aggr_raw_TMM[nondetected_combined,]
dim(aggr_not_detected_Seurat)
dim(aggr_nondetected_threshold_genes)

#aggr_nondetected_all <- rbind(aggr_nondetected_threshold_genes, aggr_not_detected_Seurat)
aggr_nondetected_all <- na.omit(aggr_nondetected_threshold_genes)
dim(aggr_nondetected_all)


#NNLS_reduce_store <- NNLS_reduce
NNLS_reduce <- read.table('Data/NNLS_average_across_100_bootstraps.log1p.30Cells.231211.tsv')


dim(bulk_data)
dim(NNLS_reduce)



bulk_raw_TMM_prc <- bulk_raw_TMM[intersect(rownames(bulk_raw_TMM), rownames(ws289_prc)),]

waldo::compare(colnames(bulk_raw_TMM_prc), colnames(NNLS_reduce))


TMM_cor_s_contaminant_genes <- pbsapply(rownames(bulk_raw_TMM_prc), function(gene){
  apply(NNLS_reduce[2:nrow(NNLS_reduce),colnames(bulk_raw_TMM_prc)], 1, function(tissue){
    stats::cor.test(bulk_raw_TMM_prc[gene,colnames(NNLS_reduce)], as.numeric(tissue), method = 'pearson')$estimate
  })
}) |> t()

TMM_cor_s_contaminant_genes <- na.omit(TMM_cor_s_contaminant_genes)



max_TMM_contaminant_corr_prc <- data.frame(row.names = rownames(TMM_cor_s_contaminant_genes),
                                             corr = MatrixGenerics::rowMaxs(as.matrix(TMM_cor_s_contaminant_genes)))




sum(rownames(max_TMM_contaminant_corr_prc)%in% ws289_prc$gene_id)
rownames(max_TMM_contaminant_corr_prc)
intersect(rownames(max_TMM_contaminant_corr_prc), ws289_prc$gene_id)


### subset the correlation to contaminants ----
TMM_cor_s_contaminant_prc <- max_TMM_contaminant_corr_prc[intersect(rownames(max_TMM_contaminant_corr_prc), ws289_prc$gene_id),, drop = F]

TMM_cor_s_contaminant_prc_nondetect_all <- TMM_cor_s_contaminant_prc[intersect(rownames(TMM_cor_s_contaminant_prc), 
                                                                                   rownames(aggr_nondetected_all)),,drop=F]


### plot the maximum correlation per gene
data.frame( max_corr = apply(TMM_cor_s_contaminant_prc_nondetect_all, 1, max)) |> ggplot() + geom_density(aes(x = max_corr))


### fit two gaussians

#devtools::install_github("yuliadm/mixComp")
library(mixComp)



#~~~~ details mixtools ----

val <- TMM_cor_s_contaminant_prc_nondetect_all$corr
lambda = .2
mu = c(0.1, 1.1)
sig = c(0.8, 0.1)
mixt <- mixtools::normalmixEM2comp(x=val,
                                   lambda = lambda,
                                   mu = mu,
                                   sigsqrd = sig, verb = T)


plot(mixt, density = TRUE, whichplots = 2)


d1 <- function(x) mixt$lambda[1]*dnorm(x, mean = mixt$mu[1], sd = mixt$sigma[1])
d2 <- function(x) mixt$lambda[2]*dnorm(x, mean = mixt$mu[2], sd = mixt$sigma[2])


bayestestR::auc(seq(0,1,0.01),d1(seq(0,1,0.01)))
bayestestR::auc(seq(0,1,0.01),d2(seq(0,1,0.01)))
bayestestR::auc(seq(0,.1757959,0.00000001),d2(seq(0,.1757959,0.00000001)))/bayestestR::auc(seq(0,1,0.01),d2(seq(0,1,0.01)))

plot(seq(0,1,0.01),d2(seq(0,1,0.01)))

temp_nc_corr.df <- data.frame(corr = TMM_cor_s_contaminant_prc_nondetect_all$corr)

values_ <- seq(min(temp_nc_corr.df$corr), max(temp_nc_corr.df$corr),length = nrow(temp_nc_corr.df)/10)

df <- data.frame(values = values_,
                 d1 = d1(values_),
                 d2 = d2(values_)
)


contaminant_correlation_threshold = .1757959

genes_low_max_contaminant_corr <- rownames(max_TMM_contaminant_corr_prc)[max_TMM_contaminant_corr_prc[,1] < contaminant_correlation_threshold]




### panel A

ggplot() + geom_density(data = temp_nc_corr.df, aes(x = corr, y = , fill = ''), alpha = 0.4) +
  scale_fill_manual(values = ('purple3')) +
  geom_line(data = df, aes(x = values, y = d1), color = 'blue', linetype = 'dashed') +
  geom_line(data = df, aes(x = values, y = d2), color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) + theme(legend.position = '') +
  geom_vline(xintercept = .1757959, color = 'darkred', size = 2) +
  ylab('density') +
  xlab('Highest Observed correlation to contaminant tissue, per gene\nCandidate missing genes only')
ggsave('figures/ncRNA/supp-low expressed genes/Max_correlation_to_contaminants_ProteinCoding_not_Detected_singleCell_240616.pdf')

### keep just the genes with a low correlation
low_corr_genes_keep <- rownames(TMM_cor_s_contaminant_prc_nondetect_all)[apply(TMM_cor_s_contaminant_prc_nondetect_all, 
                                                                                 1, max) < .1757959]



### Check that most non-neuronal genes are removed by thresholding on the correlation to contaminants

likely_non_neuronal <- read.table('references/simpleMine_non_neuronal_genes_101521.tsv', sep = '\t',header = T)
likely_non_neuronal <- likely_non_neuronal[likely_non_neuronal$Putative.match.=='Yes',]

likely_non_neuronal_gt <- data.frame(row.names = likely_non_neuronal$WormBase.Gene.ID, 
                                     matrix(0, ncol=ncol(aggr_raw_TMM), nrow=nrow(likely_non_neuronal)))
colnames(likely_non_neuronal_gt) <- colnames(aggr_raw_TMM)

nn_genes <- rownames(likely_non_neuronal_gt)

length(intersect(nn_genes, rownames(aggr_nondetected_all)))
low_expression_nn_genes2 <- intersect(nn_genes, low_corr_genes_keep)
length(low_expression_nn_genes2)

neurons <- colnames(aggr_raw_TMM)

### set up ground truth
low_expression_likely_non_neuronal_gt_cut2 <- likely_non_neuronal_gt[low_expression_nn_genes2, neurons]
aggr_raw_TMM_low2_nnplot <- aggr_raw_TMM[low_expression_nn_genes2, neurons]


### calculate FPR rate for the remaining non-neuronal genes
diags_aggr_raw_TMM_low2_nnplot <- tibble(threshold = seq(5,7, 0.0001),
                                         FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_low2_nnplot, 
                                                                           low_expression_likely_non_neuronal_gt_cut2, .x)),
                                         counts = "aggr_raw_TMM_low2_nnplot")

ggplot(diags_aggr_raw_TMM_low2_nnplot[diags_aggr_raw_TMM_low2_nnplot$threshold < 100,], aes(x = threshold, y = FPR)) + 
  geom_point() +xlim(0,100) + geom_hline(yintercept = 0) + theme_classic(base_size = 20) +
  xlab('TMM threshold') + ylab('Non-Neuronal FPR')


## find a threshold that excludes all known non-neuronal genes
diags_aggr_raw_TMM_low2_nnplot[diags_aggr_raw_TMM_low2_nnplot$FPR == 0,]


### calculate the new genes per cell type ----
new_genes <- sapply(colnames(aggr_nondetected_all), function(cell){
  sum(aggr_nondetected_all[low_corr_genes_keep,cell] >= 6.01)
})
new_genes
quantile(new_genes)

mean(new_genes)
conf.model <- lm(new_genes ~ 1)
mean(new_genes) - confint(conf.model)[[1]]


new_genes_list
new_genes_list <- sapply(colnames(aggr_nondetected_all), function(cell){
  low_corr_genes_keep[aggr_nondetected_all[low_corr_genes_keep,cell] > 4.59]
})
new_genes_list$ADL %>% clipr::write_clip()



## bar plot ----
data.frame(new_genes,
           cell_type = names(new_genes)) |>
  ggplot() +
  geom_col(aes(x = cell_type, y = new_genes)) +
  theme_classic(base_size = 15)+ theme(axis.text.x = element_text(angle = 60, hjust = 1, face = 'bold')) +
  xlab('Cell Type') + ylab('"New" Genes Detected in Bulk\n(Never Detected in Single Cell)')
ggsave('genes_never_seen_in_singleCell_barplot.pdf', width = 9)



###### ----




## calling gene expression ----

GT_thresholded_FDR_M <- sapply(colnames(aggr_bulk_TMM), function(cell){
  x1 <- rownames(aggr_bulk_TMM[sc_unexpr_1_list[[cell]],][aggr_bulk_TMM[sc_unexpr_1_list[[cell]], cell] > 45.3,])
  x2 <- rownames(TMM_cor_s_contaminant_prc)
  x3 <- intersect(x1, x2)
  length(x3)
})

GT_thresholded_FDR_M.gene.list <- sapply(colnames(aggr_bulk_TMM), function(cell){
  x1 <- rownames(aggr_bulk_TMM[sc_unexpr_1_list[[cell]],][aggr_bulk_TMM[sc_unexpr_1_list[[cell]], cell] > 45.3,])
  x2 <- rownames(TMM_cor_s_contaminant_prc)
  x3 <- intersect(x1, x2)
  return(x3)
})

names(GT_thresholded_FDR_M.gene.list) <- colnames(aggr_bulk_TMM)


sapply(GT_thresholded_FDR_M.gene.list, length)

GT_thresholded_FDR_M.gene.list$PHA |> clipr::write_clip()

quantile(GT_thresholded_FDR_M)


GT_thresholded_FDR_M.df <- data.frame(protein_genes = as.numeric(GT_thresholded_FDR_M),
                                      cells = names(GT_thresholded_FDR_M))

GT_thresholded_FDR_M.gene.list$OLQ %>% clipr::write_clip(.)

sapply(GT_thresholded_FDR_M.gene.list, length)

sc_size <- read.csv('references/sc_size_032322.csv')
cell_names <- rownames(sc_size)
sc_size <- sc_size$x
names(sc_size) <- cell_names

sc_size[names(GT_thresholded_FDR_M.gene.list)][order(sc_size[names(GT_thresholded_FDR_M.gene.list)])]

IL1_show_TMM <- aggr_bulk_TMM[GT_thresholded_FDR_M.gene.list$IL1 ,]


rownames(IL1_show_TMM) <- wbData::i2s(GT_thresholded_FDR_M.gene.list$IL1 , wb_load_gene_ids('277'))
#View(IL1_show_TMM[,c('AFD', 'I5', 'IL1', 'IL2')])


gene_ <- 'WBGene00010982'
rbind(aggr_bulk_TMM[gene_,], CeNGEN_TPM[gene_,colnames(aggr_bulk_TMM)])


conf.model <- lm(GT_thresholded_FDR_M ~ 1)
confint(conf.model)
mean(GT_thresholded_FDR_M) - confint(conf.model)[[1]]

##### bar plot of "new" genes ----
GT_thresholded_FDR_M.df %>% ggplot(data = ., aes(x = cells, y = protein_genes)) + geom_col() +
  #scale_y_continuous(transform = 'log10') +
  theme_classic(base_size = 15) + theme(axis.text.x = element_text(angle = 60, hjust = 1, face = 'bold', color = 'black')) +
  xlab('Cell Type') + ylab('Bulk Genes undectected in specific single cell clusters') 
ggsave('figures/ncRNA/supp-low expressed genes/New_LowExpressed_Protein_coding_genes_240630.pdf', width = 9, height = 6)




new_genes.vs.size.df <- data.frame(row.names = GT_thresholded_FDR_M.df$cells,
                                   raw_protein_genes = GT_thresholded_FDR_M.df$protein_genes,
                                   sc_size = sc_size[GT_thresholded_FDR_M.df$cells])
new_genes.vs.size.df['OLL',]

new_genes.vs.size.df %>%
  ggplot() +
  geom_point(aes(x = log10(sc_size), log10(raw_protein_genes)), size = 5) +
  geom_smooth(aes(x = log10(sc_size), log10(raw_protein_genes)), method = 'lm') 


library(nlstools)
nls_new_genes_decay <- nls((raw_protein_genes) ~ m + (M-m)*exp(-sc_size/alpha), data = new_genes.vs.size.df,
                           list(M = 150, alpha = 10, m = 18))


boots_nls_new_genes_decay <- nlsBoot(nls_new_genes_decay, niter=2000)
boots_nls_new_genes_decay_conf <- nlsBootPredict(boots_nls_new_genes_decay, interval = 'confidence')


new_genes.df <- cbind(new_genes.vs.size.df, boots_nls_new_genes_decay_conf)

new_genes.df[order(new_genes.df$sc_size),]
new_genes.df$cell <- rownames(new_genes.df)
ggplot(new_genes.df) +
  theme_classic(base_size = 25) +
  geom_point(aes(x=sc_size, y = raw_protein_genes), size = 3) +
  geom_line(aes(x=sc_size, y= Median), color = "red3", linetype = "dashed") +
  geom_line(aes(x=sc_size, y= `2.5%`), color = "blue3", linetype = "dashed") +
  geom_line(aes(x=sc_size, y= `97.5%`), color = "blue3", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(x=sc_size, y = raw_protein_genes, label = cell), max.overlaps = 2) +
  xlab("Single Cell Cluster size") +
  ylab("'new' protein coding genes found in bulk")
ggsave('figures/ncRNA/supp-low expressed genes/new_genes_versus_sc_cluster_size_240803.pdf', width = 8, height = 8)

new_genes_decay












