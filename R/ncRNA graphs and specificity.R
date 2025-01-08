## ncRNA
library(wbData)
library(edgeR)
library(MatrixGenerics)
library(tidyverse)
library(reshape)
library(ComplexHeatmap)
library(dendsort)
library(ggbeeswarm)
library(pbapply)
library(bayestestR)

### PEM function ----

#adapted from Kryuchkova-Mostacci, et al., 2017

# input is a genes x cell types matrix. I used the arithmetic mean across samples within cell types
fPem <- function(x){
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score

    x[x<1] <- 1
    x <- log10(x)

    x<- abs(x)
    res <- x[-nrow(x),-ncol(x)]
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}


#setwd('~/Bioinformatics/bsn5/')



### load data, normalize ----

bulk_data <- read.table('Data/Bulk_data_bsn12_231211.tsv')

NNLS_reduce <- read.table('Data/NNLS_average_across_100_bootstraps.log1p.30Cells.231211.tsv')


NNLS_reduce |> Heatmap(cluster_columns = F, row_order = order(rowMeans(NNLS_reduce), decreasing = T))

dim(bulk_data)
dim(NNLS_reduce)


bulk_meta <- read.table('Data/bulk_bsn12_metadata.tsv', sep = '\t')

bulk_meta <- bulk_meta[rownames(bulk_data),]


bulk_data_pk <- (bulk_data/bulk_meta$Length) * 1000


bulk_raw_GeTMM <- DGEList(bulk_data_pk, group = str_split_fixed(colnames(bulk_data_pk), 'r', 2)[,1])
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM)
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)

dim(bulk_raw_GeTMM)


### get average profile per cell type ----
dim(bulk_raw_GeTMM)
aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x)
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))


## get reference, subset to ncRNAs of interest ----


ws289 <- wb_load_gene_ids('289') %>% data.frame(row.names = .$gene_id, .)

unique(ws289$biotype)

ws289_allnc <- ws289[ws289$biotype != 'gene' & ws289$biotype != 'protein_coding_gene' &
                       ws289$biotype != 'miRNA_gene' & ws289$biotype != 'piRNA_gene' &
                       ws289$biotype != 'antisense_lncRNA_gene' & ws289$biotype != 'transposable_element_gene' &
                       ws289$biotype != 'rRNA_gene' & ws289$biotype != 'scRNA_gene', ]


ws289_possible_artifacts <- read.table('../references/ws289_CeNGEN_marker_possible_artifact_genes.tsv', sep = '\t')


aggr_nc_raw_GeTMM <- aggr_raw_GeTMM[rownames(aggr_raw_GeTMM) %in% rownames(ws289_allnc) ,]
dim(aggr_nc_raw_GeTMM)

aggr_nc_raw_GeTMM <- aggr_nc_raw_GeTMM[!rownames(aggr_nc_raw_GeTMM) %in% ws289_possible_artifacts$gene_id,]


### calculate PEM scores, genes x cell types

GeTMM_PEM <- fPem(aggr_nc_raw_GeTMM)
GeTMM_PEM <- data.frame(GeTMM_PEM)

GeTMM_PEM_ncRNA <- GeTMM_PEM[intersect(rownames(GeTMM_PEM), rownames(ws289_allnc)),]

### calculate correlation to contaminants

NNLS_reduce <- NNLS_reduce[,intersect(colnames(NNLS_reduce), colnames(bulk_raw_GeTMM))]

bulk_raw_GeTMM_nc <- bulk_raw_GeTMM[intersect(rownames(bulk_raw_GeTMM), rownames(ws289_allnc)),]

bulk_raw_GeTMM_nc <- bulk_raw_GeTMM_nc[!rownames(bulk_raw_GeTMM_nc) %in% ws289_possible_artifacts$gene_id,]

waldo::compare(colnames(bulk_raw_GeTMM_nc), colnames(NNLS_reduce))


## for each gene calculate the correlation of the contaminant gene expression to the contamination score from the subsampled NNLS
GeTMM_cor_s_contaminant_noncoding_genes <- t(pbsapply(rownames(bulk_raw_GeTMM_nc), function(gene){
  apply(NNLS_reduce[2:nrow(NNLS_reduce),colnames(bulk_raw_GeTMM_nc)], 1, function(tissue){
    stats::cor.test(bulk_raw_GeTMM_nc[gene,colnames(NNLS_reduce)], as.numeric(tissue), method = 'pearson')$estimate
  })
}))

#remove NA genes (445, 3.9%)
GeTMM_cor_s_contaminant_noncoding_genes <- na.omit(GeTMM_cor_s_contaminant_noncoding_genes)


### this is to generate an example figure of a few correlation estimates for all tissues for a few genes
# GeTMM_cor_s_contaminant_noncoding_genes |> head(3) |>
#   reshape2::melt() |>
#   mutate(`contaminant correlation` = value) |>
#   mutate(`contaminant tissue` = Var2) |>
#   mutate(genes = Var1) |>
#   ggplot() +
#   geom_point(aes(x = genes, y = `contaminant correlation`, color = `contaminant tissue`), size = 8) +
#   theme_classic(base_size = 20) +
#   theme(axis.text.x = element_blank(),
#         legend.position = '')


## take the highest correlation to a contaminant tissue type for each gene
max_GeTMM_contaminant_corr_noncoding <- data.frame(row.names = rownames(GeTMM_cor_s_contaminant_noncoding_genes),
                                         MatrixGenerics::rowMaxs(as.matrix(GeTMM_cor_s_contaminant_noncoding_genes)))



max_GeTMM_contaminant_corr_noncoding_genes <- intersect(rownames(GeTMM_cor_s_contaminant_noncoding_genes),
                                                        rownames(ws289_allnc))
max_GeTMM_contaminant_corr_noncoding <- data.frame(row.names = max_GeTMM_contaminant_corr_noncoding_genes,
                                                   ncRNA_corr = max_GeTMM_contaminant_corr_noncoding[max_GeTMM_contaminant_corr_noncoding_genes,1])

#### calculate binary expression using static threshold ----

nc_GeTMM_bin <- sapply(unique(str_split_fixed(colnames(bulk_raw_GeTMM_nc), 'r', 2)[,1]), function(cell){
  print(cell)
  bulk <- bulk_raw_GeTMM_nc[,str_split_fixed(colnames(bulk_raw_GeTMM_nc), 'r', 2)[,1]==cell]
  bulk <- bulk[intersect(rownames(bulk),
                         rownames(ws289_allnc)),]
  bin <- bulk > 5 ## greater than 5 GeTMM
  bin <- Matrix::rowMeans(bin)
  return((bin > 0.65) * 1) ## greater than threshold in >65% of samples
})

colSums(nc_GeTMM_bin)

ncRNA_expressed_somewhere <- rownames(nc_GeTMM_bin[Matrix::rowSums(nc_GeTMM_bin) > 0,])
length(ncRNA_expressed_somewhere)



### fit two gaussians to the maximum correlation to contaminants

#devtools::install_github("yuliadm/mixComp")
library(mixComp)
max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,]



#~~~~ details mixtools ----

val <- max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,]
lambda = .2
mu = c(0.1, 1.1)
sig = c(0.8, 0.2)
mixt <- mixtools::normalmixEM2comp(x=val,
                                   lambda = lambda,
                                   mu = mu,
                                   sigsqrd = sig, verb = T)


plot(mixt, density = TRUE, whichplots = 2)


d1 <- function(x) mixt$lambda[1]*dnorm(x, mean = mixt$mu[1], sd = mixt$sigma[1])
d2 <- function(x) mixt$lambda[2]*dnorm(x, mean = mixt$mu[2], sd = mixt$sigma[2])


bayestestR::auc(seq(0,1,0.01),d1(seq(0,1,0.01)))
bayestestR::auc(seq(0,1,0.01),d2(seq(0,1,0.01)))
bayestestR::auc(seq(0,0.221,0.00001),d2(seq(0,0.221,0.00001)))/bayestestR::auc(seq(0,1,0.01),d2(seq(0,1,0.01)))


temp_nc_corr.df <- data.frame(ncRNA_corr = max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,])


## generate values to plot for the distributions, lazy conversion for ggplot
values_ <- seq(min(temp_nc_corr.df$ncRNA_corr), 
               max(temp_nc_corr.df$ncRNA_corr),
               length = nrow(temp_nc_corr.df)/10)

df <- data.frame(values = values_,
                 d1 = d1(values_),
                 d2 = d2(values_)
)

## set threshold for maximum correlation to contaminants
contaminant_correlation_threshold = 0.221


## select genes that are below that threshold
genes_low_max_contaminant_corr <- rownames(max_GeTMM_contaminant_corr_noncoding)[max_GeTMM_contaminant_corr_noncoding[,1] < contaminant_correlation_threshold]




### panel A

ggplot() + geom_density(data = temp_nc_corr.df, aes(x = ncRNA_corr, y = , fill = ''), alpha = 0.4) +
  scale_fill_manual(values = ('purple3')) +
  geom_line(data = df, aes(x = values, y = d1), color = 'blue', linetype = 'dashed') +
  geom_line(data = df, aes(x = values, y = d2), color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) + theme(legend.position = '') +
  geom_vline(xintercept = 0.221, color = 'darkred', size = 2) +
  ylab('density') +
  xlab('Highest Observed correlation to contaminant tissue, per gene\nCandidate missing genes only')

ggsave('figures/noncoding RNA correlation fit 051524.pdf', width = 8, height = 6)






#### genes passing thresholds ----
ncRNA_corr_cut.df <- pbsapply(colnames(aggr_nc_raw_GeTMM), function(cell){

  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]

  bulk <- bulk_raw_GeTMM_nc[, samples_cell_type %in% c(cell)]


  sapply(unlist(unique(ws289_allnc$biotype)), function(ncRNA){
    #print(genes)
    genes <- rownames(ws289_allnc[ws289_allnc$biotype==ncRNA,])
    genes <- intersect(genes, rownames(bulk_raw_GeTMM_nc))
    genes <- intersect(genes, genes_low_max_contaminant_corr)
    #print('genes')
    bulk_bin <- bulk[genes,] > 5

    bulk_bin_ave <- rowMeans(bulk_bin)

    sum(bulk_bin_ave > 0.65)


  })
})

ncRNA_corr_cut.df



apply(ncRNA_corr_cut.df, 1, mean)
apply(ncRNA_corr_cut.df, 1, function(type){

  mean(type) - confint(lm(type ~ 1))[[1]]
})

apply(ncRNA_corr_cut.df, 2, sum)

mean(apply(ncRNA_corr_cut.df, 2, sum))
mean(apply(ncRNA_corr_cut.df, 2, sum)) - confint(lm(apply(ncRNA_corr_cut.df, 2, sum)~1))


#### panel B

reshape::melt(ncRNA_corr_cut.df) %>%
  data.frame(ncRNA_type = .[,1],
             cell_type = .[,2],
             total = .[,3]) %>%
  ggplot() + geom_col(aes(x = cell_type, y = total, fill = ncRNA_type)) +
  scale_y_continuous(breaks = c(0,400,800,1200), labels = c('0','400','800','1200')) +
  theme_classic(base_size = 20) + xlab('Cell Type') +
  theme(axis.text.x = element_text(face = 'bold', angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
        axis.text.y = element_text(face = 'bold', color = 'black'),
        axis.title = element_text(face = 'bold', color = 'black'))
ggsave('figures/noncoding RNAs by type barplot 051624.pdf', width = 20, height = 15)

cell = 'AFD'
samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM_nc), 'r', 2)[,1]
bulk <- bulk_raw_GeTMM_nc[, samples_cell_type %in% c(cell)]

genes <- intersect(rownames(bulk_raw_GeTMM_nc), genes_low_max_contaminant_corr)
bulk <- bulk[genes,]

bulk_bin <- bulk[genes,] > 5

bulk_bin_ave <- rowMeans(bulk_bin)

ncRNA_corr_cut.list.df <- lapply(colnames(aggr_nc_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM_nc), 'r', 2)[,1]
  bulk <- bulk_raw_GeTMM[, samples_cell_type %in% c(cell)]

  genes <- intersect(rownames(bulk_raw_GeTMM_nc), genes_low_max_contaminant_corr)
  bulk <- bulk[genes,]

  bulk_bin <- bulk[genes,] > 5

  bulk_bin_ave <- rowMeans(bulk_bin)


  return(names(bulk_bin_ave[bulk_bin_ave > 0.65]))
})



names(ncRNA_corr_cut.list.df) <- colnames(aggr_nc_raw_GeTMM)
ncRNA_corr_cut.list.df

ncRNA_corr_cut.list.df_table <- table(unlist(ncRNA_corr_cut.list.df))
ncRNA_corr_cut.list.df_table[ncRNA_corr_cut.list.df_table > 40]

ncRNA_corr_cut.list.df_old <- ncRNA_corr_cut.list.df[!names(ncRNA_corr_cut.list.df) %in% c('CEP', 'SIA', 'PVQ', 'PVP', 'DVB', 'RME', 'HSN')]

unique(ncRNA_corr_cut.list.df_old |> unlist()) |> length()

unique(ncRNA_corr_cut.list.df |> unlist()) |> length()


## how many ubiquitous ncRNAs:
sum(ncRNA_corr_cut.list.df_table >= 37)


# panel C
data.frame(hist = seq(1,41,1),
           totals = sapply(seq(1,41,1), function(threshold){
  sum(ncRNA_corr_cut.list.df_table == threshold)
})) %>% ggplot() + geom_col(aes(x = hist, y = log10(totals)), color = 'black', fill = 'grey') +
  scale_y_continuous(breaks = seq(0,3,1), labels = c('0', '10', '100', '1000')) +
  geom_vline(xintercept = 36.5, col = 'darkred', size = 1) +
  xlab('Total Cell types expressing ncRNA') + ylab('total ncRNAs') +
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text = element_text(face = 'bold', size = 20), axis.title = element_text(face = 'bold', size = 25)
        )
ggsave('figures/pan-neuronal noncoding RNA histogram 051624.pdf')


ncRNA_corr_cut.list.df_table
highly_specific_ncRNAS <- aggr_nc_raw_GeTMM[names(ncRNA_corr_cut.list.df_table),][apply(GeTMM_PEM_ncRNA[names(ncRNA_corr_cut.list.df_table),],
                                                                                        1, max) > 0.65,]
highly_specific_ncRNAS <- highly_specific_ncRNAS[rowSums(highly_specific_ncRNAS > 3) <= 10,]
pan_ncRNAs <- aggr_nc_raw_GeTMM[names(ncRNA_corr_cut.list.df_table[ncRNA_corr_cut.list.df_table >=37]),]

dim(highly_specific_ncRNAS)


sum(highly_specific_ncRNAS$RIM > 6)

RIM_specific_genes <- highly_specific_ncRNAS |> filter(RIM > 6) |> rownames()


ws289_gene_coords <- wbData::wb_load_gene_coords('289')

ws289_gene_coords_RIM <- ws289_gene_coords |> filter(gene_id %in% RIM_specific_genes)

ws289_gene_coords_RIM

ws289_gene_coords_RIM <- ws289_gene_coords_RIM[, 1:3]

ws289_gene_coords_RIM_boot <- rsample::bootstraps(ws289_gene_coords_RIM, times = 5000)

ws289_gene_coords_RIM_boot$splits[1]

ws289_gene_coords_RIM_boot$splits[[1]] |> as.data.frame()

ws289_gene_coords_RIM_boot_per_chr <- t(pbsapply(1:5000, function(x){

  tmp <- ws289_gene_coords_RIM_boot$splits[[x]] |> as.data.frame()

  I = sum(tmp$chr == 'I')
  II = sum(tmp$chr == 'II')
  III = sum(tmp$chr == 'III')
  IV = sum(tmp$chr == 'IV')
  V = sum(tmp$chr == 'V')
  X = sum(tmp$chr == 'X')

  vec <- c(I, II, III, IV, V, X)
  names(vec) <- c('I', 'II', 'III', 'IV', 'V', 'X')

  return(vec)

}))

summary(ws289_gene_coords_RIM_boot_per_chr)


apply(ws289_gene_coords_RIM_boot_per_chr, 2, coxed::bca)

ws289_gene_coords_RIM.df <- cbind(table(ws289_gene_coords_RIM$chr) |> data.frame(), t(apply(ws289_gene_coords_RIM_boot_per_chr, 2, coxed::bca)))
colnames(ws289_gene_coords_RIM.df) <- c('chr', 'Freq', 'lower.ci', 'upper.ci')

ws289_gene_coords_RIM.df$expected_Freq <- sum(ws289_gene_coords_RIM.df$Freq)/6

ws289_gene_coords_RIM.df |>
  ggplot(aes(x = chr, y = Freq)) +
  geom_col() +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.2,
                position=position_dodge(.9)) +
  xlab('Chromosome') + ylab('RIM specific ncRNA genes') +
  theme_classic(base_size = 25)
ggsave('figures/RIM specific ncRNA chromosomes 051624.pdf')



ws289_gene_coords_RIM.df <- table(ws289_gene_coords_RIM$chr) |>
  data.frame()
colnames(ws289_gene_coords_RIM.df) <- c('chr', 'Freq')

data.frame(chr = names(table(ws289_gene_coords_RIM$chr)),
           genes = table(ws289_gene_coords_RIM$chr))


ws289_gene_coords_RIM_X <- ws289_gene_coords_RIM |> filter(chr == 'X')

ws289_gene_coords_RIM_X$start

ws289_gene_coords_RIM_X |>
  ggplot(aes( x = start )) +
  geom_histogram()



#pan_ncRNAs

library(circlize)


## set up colors ----
col_fun = colorRamp2(c(0, 3), c("#DEDCDC", "#Bb361d"))

col_fun2 = colorRamp2(c(0, 5), c("#DEDCDC", "#Bb361d"))

tmp.df <- t(scale(t(log10(1+highly_specific_ncRNAS))))
row_dend = dendsort(hclust(dist(tmp.df)))
col_fun = colorRamp2(c(0, 3), c("#DEDCDC", "#920000"))
#col_fun = colorRamp2(c(min(unlist(tmp.df)),0, max(unlist(tmp.df))), c('blue', "#DEDCDC", "#Bb361d"))

dim(highly_specific_ncRNAS)
ws289[rownames(highly_specific_ncRNAS), 'biotype'] %>% clipr::write_clip()


## subset to highly specific genes in reference

highly_specific_ncRNAS_ncRNA_genes <- highly_specific_ncRNAS[rownames(highly_specific_ncRNAS) %in% ws289_allnc[ws289_allnc$biotype=='ncRNA_gene', 'gene_id'], ]

highly_specific_ncRNAS_ncRNA_genes[order(highly_specific_ncRNAS_ncRNA_genes$ADL, decreasing = T),]


### plot the highly specific genes ----
library(tidyHeatmap)

highly_specific_ncRNAS_tidy <- log10(1+t(highly_specific_ncRNAS)) |>
  as_tibble(rownames="Cell_type") |>
  pivot_longer(cols = -Cell_type, names_to = 'gene', values_to = 'value')

highly_specific_ncRNAS_tidy$biotype <- ws289[highly_specific_ncRNAS_tidy$gene,'biotype']
highly_specific_ncRNAS_tidy$Modality <- ws289[highly_specific_ncRNAS_tidy$gene,'biotype']

summary(highly_specific_ncRNAS_tidy$value)

highly_specific_ncRNAS_tidy$biotype |> unique()

highly_specific_ncRNAS_tidy |>
  group_by(biotype) |>
  tidyHeatmap::heatmap(gene, Cell_type, value,
          col = col_fun, 
          palette_grouping = list(c("#009292", "#ff6db6", "#b66dff", "#006ddb", "#920000")),
          column_names_gp = gpar(fontsize = 15, fontface = 'bold'),
          row_dend_width = unit(0, "cm")) #|> save_pdf('figures/specific_ncRNAs_heatmap_051624.pdf', width = 8, height = 8)




highly_specific_ncRNAS_ncRNA_genes.names <- highly_specific_ncRNAS_ncRNA_genes
rownames(highly_specific_ncRNAS_ncRNA_genes.names) <- i2s(rownames(highly_specific_ncRNAS_ncRNA_genes.names), ws289)



### plot the "pan-neuronal" genes ----
hmap2 <- Heatmap(
  log10(1+pan_ncRNAs),
  col = col_fun2,
  name = " ",
  show_row_names = FALSE,
  show_column_names = T,
  cluster_rows = T,
  cluster_columns = F,
  #column_order = order(neuron_meta_cut.df$Modality_collapsed),
  show_column_dend = F,
  show_row_dend = F,
  column_names_gp = gpar(fontsize = 15, fontface = 'bold'),
  #top_annotation=colAnn
  )

pdf("pan_ncRNAs_heatmap.pdf",width=8,height=5)
draw(hmap2, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()




nrow(highly_specific_ncRNAS)


### plot max PEM specificity scores ----
max_pem.df <- apply(GeTMM_PEM_ncRNA[rowSums(GeTMM_PEM_ncRNA > 3) < 10 & rownames(GeTMM_PEM_ncRNA) %in% names(ncRNA_corr_cut.list.df_table),],
      1, max) %>% data.frame(max_PEM = .)

ggplot(max_pem.df) + geom_density(aes(x = (max_PEM)), color = 'black', fill = '#77bfc2', binwidth = 0.01) +
  geom_vline(xintercept = (0.65), size = 1, color = 'darkred') +
  theme_classic(base_size = 20) +
  xlab('Highest per gene PEM specificity score')
ggsave('noncoding RNA PEM scores.pdf')




#### how many highly specific genes per cell type

#highly_specific_ncRNAS_expressed <- sa
highly_specific_ncRNAS

Specific_ncRNA_per_celltype <- data.frame(ncRNA = sapply(colnames(highly_specific_ncRNAS), function(cell){
  specific_genes <- rownames(highly_specific_ncRNAS)
  sum(GeTMM_PEM_ncRNA[specific_genes,cell] > 0.6)
  }),
           cell = names(highly_specific_ncRNAS)
           #modality = sapply(names(highly_specific_ncRNAS), function(cell){neuron_meta_cut.df[cell, 'Modality_collapsed']})
  )

mean(Specific_ncRNA_per_celltype$ncRNA)
mean(Specific_ncRNA_per_celltype$ncRNA) - confint(lm(Specific_ncRNA_per_celltype$ncRNA ~ 1))[[1]]


Specific_ncRNA_per_celltype |> 
  mutate(cell = factor(cell, levels = cell[order(ncRNA, decreasing = T)])) |>
  ggplot() +
  geom_col(aes(x = cell, y = ncRNA))
  



#### plot the number of pan neuronal and highly specific genes by RNA classs, pie charts ----
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

RNA_Classes <- unique(ws289_allnc$biotype)
specific_genes <- rownames(highly_specific_ncRNAS)

specific_genes_class <- sapply(RNA_Classes, function(class){
  sum(ws289_allnc[specific_genes, 'biotype'] == class)
})
specific_genes_class/sum(specific_genes_class)
sum(specific_genes_class)
specific_genes_class_proportions <- sapply(RNA_Classes, function(class){
  sum(ws289_allnc[specific_genes, 'biotype'] == class)/length(specific_genes)
})

specific_genes_class_proportions <- data.frame(proportions = specific_genes_class_proportions,
           group = names(specific_genes_class_proportions))
table(ws289_allnc$biotype)['tRNA_gene']/sum(table(ws289_allnc$biotype))

pct <- round(specific_genes_class_proportions$proportions*100, digits = 1)
lbls <- specific_genes_class_proportions$group
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(specific_genes_class_proportions$proportions, labels = lbls, col=brewer.pal(length(lbls), 'Paired'),
    main="RNA classes in Cell Type Specific noncoding RNAs")


dim(pan_ncRNAs)
pan_genes <- rownames(pan_ncRNAs)
pan_genes_class_proportions <- sapply(RNA_Classes, function(class){
  sum(ws289_allnc[pan_genes, 'biotype'] == class)/length(pan_genes)
})

pan_genes_class_proportions <- data.frame(proportions = pan_genes_class_proportions,
                                               group = names(pan_genes_class_proportions))

pan_genes_class_proportions

pan_genes_classes <- sapply(unique(ws289_allnc$biotype), function(type){
  sum(ws289_allnc[pan_genes, 'biotype'] == type)
})




pan_genes_classes[order(names(pan_genes_classes))]/table(ws289_allnc$biotype)
pan_genes_classes

pan_genes_classes
sum(table(ws289_allnc$biotype))-table(ws289_allnc$biotype)['lincRNA_gene']

pan_genes_classes/sum(pan_genes_classes)

pan_genes_class_proportions$expected <- sapply(pan_genes_class_proportions$group, function(RNA_type){
  table(ws289_allnc$biotype)[RNA_type]/length(ws289_allnc$biotype)
})
pan_genes_class_proportions$proportions_over_expected <- pan_genes_class_proportions$proportions/pan_genes_class_proportions$expected

library(RColorBrewer)

pct <- round(pan_genes_class_proportions$proportions*100, digits = 1)
lbls <- pan_genes_class_proportions$group
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels

pdf('figures/ncRNA/pan_genes_051624.pdf')
pie(pan_genes_class_proportions$proportions, labels = lbls, col=brewer.pal(length(lbls), 'Paired'),
    main="RNA classes in Pan-Neuronal noncoding RNAs")
dev.off()



pdf('figures/ncRNA/specific_genes_051624.pdf')
pie(specific_genes_class_proportions$proportions, labels = lbls, col=brewer.pal(length(lbls), 'Paired'),
    main="RNA classes in Cell Type Specific noncoding RNAs")
dev.off()







