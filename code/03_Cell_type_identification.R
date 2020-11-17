#!/usr/bin/env Rscript

### title: Identification of cell types in lung samples (7 COVID-19 and 1 control donors)
### authors: Yiping Wang, Jana Biermann
### date: 09/29/20

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)


# Read in file from script #2
seu <- readRDS('data/lungs_all/data_lungs_all.rds')


#### Cell type identification draft assignment ####

## Add module scores for signatures: 
# Canonical markers to identify AT1/AT2 cells, ionocytes, tuft cells, myofibroblasts, and lipofibroblasts, along with gene signatures based on differential gene expression to identify alveolar macrophages (AM), were obtained from Travaglini et al. (Travaglini et al., 2020) 
# Signatures to identify subtypes of macrophages were obtained from Mould et al. (Mould et al., 2019)
# Signatures to identify DATP cells, primed and cycling AT2 cells were obtained from Choi et al. (Choi et al., 2020)
# Signatures to identify tuft cell subtypes were obtained from Montoro et al. (Montoro et al., 2018)
# Here only used for draft assignments
sigs <- read.csv('data/covid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = T)
}

# Cell type draft assignment
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(0, 1, 5, 7, 10, 13, 23), 'Myeloid', 'NA')
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(9, 14, 16, 22, 24), 'Epithelial cells', seu@meta.data$cell_type_main)
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(8, 21), 'Endothelial cells', seu@meta.data$cell_type_main)
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(4, 11, 12, 15, 18, 25, 26), 'Fibroblasts', seu@meta.data$cell_type_main)
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(19), 'Mast cells', seu@meta.data$cell_type_main)
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(6, 20), 'B cells', seu@meta.data$cell_type_main)
seu@meta.data$cell_type_main <- ifelse(seu@meta.data$integrated_snn_res.0.8 %in% c(2, 3, 17), 'T cells', seu@meta.data$cell_type_main)


#### Overall classification assignment ####

combined = seu
userevisedtuft = TRUE
combined = combined[,!is.na(combined$celltype_bped_main)]
source("11_Add_gene_signature_module_scores.R")
seu = combined

# Save Seurat object
saveRDS(seu, paste0('data/lungs_all/data_', pat, '_with_covid_sigs.rds'))


#### Cell type identification in epithelial cells ####

# Subset
seu <- subset(seu, cell_type_main == 'Epithelial cells')

# Seurat workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(object = seu, dims = 1:20)

# Rerun module scores for signatures on subset of epithelial cells
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}


print(Sys.time())
