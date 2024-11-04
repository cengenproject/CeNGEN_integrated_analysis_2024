#### filtering out likely artifact genes

## libraries ----


library(wbData)
library(dplyr)
library(GenomicRanges)

#### read in marker genes ----


marker_list <- c("bnc-1", "cat-4", "ceh-17", "ceh-24", "ceh-28", "ceh-36", "daf-7", "dat-1", 
                 "eat-4", "F23H12.7", "rgef-1", "F49C12.10", "flp-10", "flp-11", 
                 "flp-12", "flp-13", "flp-17", "flp-19", "flp-20", "flp-3", "flp-33", 
                 "flp-1", "gcy-35", "gcy-5", "gcy-7", "gcy-8", "glr-1", "glr-3", "gpa-4", 
                 "hlh-34", "klp-6", "lad-2", "lgc-35", "lgc-8", "mec-17", "mec-3", "mgl-1",
                 "nlp-17", "nmr-1", "ocr-3", "pept-3", "rab-3", "rig-3", "rol-6", "ser-2", 
                 "srg-13", "srg-8", "srh-142", "srv-3", "srw-119", "str-1", "tbh-1", 
                 "tdc-1", "tph-1", "ttll-9", "ttr-39", "ttx-3", "unc-122", "unc-25", 
                 "unc-4", "unc-47", "unc-53", "W02A2.5")

ws289_coord <- wbData::wb_load_gene_coords('289') |> data.frame()


## verifying that all genes are named correctly in WS289 ----
marker_list[!marker_list %in% ws289_coord$name]
sum(marker_list %in% ws289_coord$name)/length(marker_list)



## generate coordinate tables for 4k regions upstream of marker genes ----

ws289 <- wb_load_gene_ids('289') %>% data.frame(row.names = .$gene_id, .)

marker_list_id <- wbData::s2i(marker_list, ws289)


marker_gene_coord <- sapply(marker_list_id, function(gene){
  
  name <- ws289_coord[ws289_coord$gene_id == gene, 'name']
  chr <- ws289_coord[ws289_coord$gene_id == gene, 'chr']
  strand <- ws289_coord[ws289_coord$gene_id == gene, 'strand']
  
  start <- ws289_coord[ws289_coord$gene_id == gene, 'start']
  end <- ws289_coord[ws289_coord$gene_id == gene, 'end']
  
  if(strand == '-'){

    list <- c('gene_id' = gene, 
              'name' = name,
              'chr' = chr,
              'start' = as.numeric(end-1000),
              'end' = as.numeric(end+4000))
  } else{
    
    list <- c('gene_id' = gene, 
                   'name' = name,
                   'chr' = chr,
                   'start' = as.numeric(start-4000),
                   'end' = as.numeric(start+1000))
  }
  
  
  
}) |> 
  t() |> 
  data.frame() |> 
  tibble() |> 
  mutate(across(4, as.integer)) |> 
  mutate(across(5, as.integer)) |>
  mutate(position = paste0(chr, ':', start, '-', end))

marker_gene_coord

write.table(marker_gene_coord, '../references/marker_gene_upstream_5k_051524.tsv', sep = '\t', quote = F, row.names = F)


marker_gene_coord <- read.table('../references/marker_gene_upstream_5k_051524.tsv', sep = '\t', header = 1)

marker_gene_coord
marker_gene_coord

ws289_coord_GR <- makeGRangesFromDataFrame(ws289_coord)
ws289_coord_GR@metadata <- ws289_coord

ws289_coord_GR@ranges@NAMES <- ws289_coord$gene_id


marker_gene_coord_GR <- makeGRangesFromDataFrame(marker_gene_coord)
marker_gene_coord_GR@ranges@NAMES <- marker_gene_coord$gene_id

tmp <- subsetByOverlaps(ws289_coord_GR, marker_gene_coord_GR)
tmp@ranges@NAMES

possible_artifact_genes <- tmp@ranges@NAMES


marker_gene_coord

ws289_coord |> filter(gene_id %in% possible_artifact_genes) |>
  write.table('../references/ws289_CeNGEN_marker_possible_artifact_genes.tsv', sep = '\t')



