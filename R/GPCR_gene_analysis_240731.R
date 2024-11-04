##### load libraries ----





#### load gene reference tables ----

GPCR_genes.df <- read.table('references/GPCR_Gene_table.csv', sep = ',', header = 1)

GPCR_pseudogenes <- read.table('references/GPCR_pseudogene_table.csv', sep = ',', header = 1)


#### load modality metadata ----

modality <- read.table('references/Hammarlund_Kaitlyn_Alec_neuron_annotation_070921.csv', sep = ',', header = T)

#### load data ----

ncRNA_expr_list <- readRDS('Data/noncoding_RNA_expression_list_240731.rds')

bulk_integrated_aggregate_threshold_2 <- read.table('Data/Integrated_thresholded/240721_Integrated_bsn12_cpm_threshold_2.csv',
                                                    sep = ',')



#####


gpcr_pseudogenes_table <- sapply(ncRNA_corr_cut.list.df, function(x){sum(x %in% 
                                                                           SR_class_pseudogenes_joint)}) |> 
  data.frame(sum = _)

gpcr_pseudogenes_table$cell <- rownames(gpcr_pseudogenes_table)

gpcr_pseudogenes_table$Modality <- sapply(gpcr_pseudogenes_table$cell, function(x){
  modality[modality$Neuron == x, 'Modality']
})


ggplot(gpcr_pseudogenes_table) + 
  geom_boxplot(aes(x = Modality, y = sum, color = Modality)) +
  ggrepel::geom_text_repel(aes(x = Modality, y = sum, color = Modality, label = cell), box.padding = 0.7,
                           color = 'black', max.overlaps = 1, fontface = 'bold') +
  geom_beeswarm(aes(x = Modality, y = sum, color = Modality), size = 3) +
  ylab('Serpentine Receptor Pseudogenes per cell') +
  ggtitle('GPCR Pseudogenes') +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', color = 'black'),
        plot.title = element_text(hjust = 0.5, face = 'bold'))  
ggsave("figures/Figure5_Serpentine_Receptor_Pseudogenes_240630.pdf", width = 8, height = 8)



bulk_integrated_aggregate_threshold_2_GPCR <- bulk_integrated_aggregate_threshold_2[GPCR_genes.df$wormbase_id,] |> na.omit()


modality['VD_DD',] <- modality['VD',]


GPCR_med_df <- sapply(colnames(bulk_integrated_aggregate_threshold_2_GPCR), function(x){
  gpcr <- sum(bulk_integrated_aggregate_threshold_2_GPCR[,x] > 0)
  m <- modality[modality$Neuron == x, 'Modality'][1]
  if(x == 'VD_DD'){m = 'Motor'}
  return(c('GPCR_total' = as.numeric(gpcr), 'modality' = m, cell = x))
}) #|> t() |> data.frame()

GPCR_med_df$GPCR_total <- as.numeric(GPCR_med_df$GPCR_total)


ggplot(data = GPCR_med_df, aes(x = modality, y = GPCR_total, color = modality)) +
  geom_boxplot() +
  ggrepel::geom_text_repel(aes(label = cell), box.padding = 0.7,
                           color = 'black', max.overlaps = 1, fontface = 'bold') +
  geom_beeswarm(size = 4) + 
  ylab('Protein Coding GPCRs Detected') +
  xlab('') +
  ggtitle('GPCR Protein Coding Genes') +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', color = 'black'),
        plot.title = element_text(hjust = 0.5, face = 'bold'))  

ggsave('figures/GPCR_Protein_per_cell_240730.pdf', height = 8, width = 8)



