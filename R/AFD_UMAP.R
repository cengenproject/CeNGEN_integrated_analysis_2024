##### get AFD UMAP
library(Seurat)

sc_object <- readRDS('032224_L4_all_cells_Seurat5.rds')

sc_object <- sc_object[,sc_object$Cell.type=='AFD']

sc_object


sc_object <- ScaleData(sc_object)
sc_object <- RunPCA(sc_object)
sc_object <- RunUMAP(sc_object, dims = 1:20, reduction = 'pca')
PCAPlot(sc_object)
sc_object@meta.data

UMAPPlot(sc_object, cols = 'WBGene00006652')
dimplot(sc_object, )


AFD_umap <- sc_object@reductions$umap@cell.embeddings |> data.frame()
ttx_1<- sc_object@assays$RNA@counts['WBGene00006652',]
AFD_umap$`ttx-1` <- ttx_1
AFD_umap$nonZero <- ((AFD_umap$`ttx-1` > 0) *1) 
AFD_umap$nonZero[AFD_umap$nonZero == 1] <- 'detected'
AFD_umap$nonZero[AFD_umap$nonZero == 0] <- 'not detected'


ggplot(AFD_umap) +
  ggrastr::geom_point_rast(aes(x = umap_1, y = umap_2, color = log1p(`ttx-1`)), size = 3) +
  scale_color_continuous(low = "grey", high = "brown") +
  ylab('') +
  xlab('') +
  theme_classic() +
  theme(axis.text = element_blank())
ggsave('figures/AFD_ttx1_logcounts_240804.pdf', width = 6, height = 5)

ggplot(AFD_umap) +
  ggrastr::geom_point_rast(aes(x = umap_1, y = umap_2, color = nonZero), size = 3) +
  scale_color_manual(values = c('#D81B60', '#1E88E5')) +
  ylab('') +
  xlab('') +
  theme_classic() +
  theme(axis.text = element_blank())
ggsave('figures/AFD_ttx1_detection_240804.pdf', width = 6, height = 5)
sum(AFD_umap$proportion)/245

ttx_1 |> length()
