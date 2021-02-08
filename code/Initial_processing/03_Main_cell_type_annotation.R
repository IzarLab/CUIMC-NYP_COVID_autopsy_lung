#!/usr/bin/env Rscript

#### Main cell type annotation of integrated Seurat object
library(dplyr)
library(Seurat)


# Read in object
seu <- readRDS('data/lungs_all/data_lungs_all_initial_obj.rds')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_rpca.csv', row.names = F)

# Main cell type assignment after manual annotation based on DGE
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(3, 4, 6, 7, 9, 10), 'Myeloid cells', 'NA')
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(2, 5, 11, 12, 15, 18, 23, 25), 'Epithelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(14, 28, 29), 'Endothelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(0, 13, 16, 17, 27, 31, 32), 'Fibroblasts', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(8, 24), 'B cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(1, 19, 20, 21), 'T cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(22), 'Neuronal cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(26), 'Mast cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(30), 'APC-like cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.8 %in% c(33), 'remove', seu$cell_type_main)

# remove cluster 33
seu <- subset(seu, integrated_snn_res.0.8 != 33)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)

# Add clinical data to meta.data
clin <- read.csv('data/clinical_data_covid.csv', na.strings = 'NA')
seu$barcode <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, clin, by = 'orig.ident')
rownames(seu@meta.data) <- seu$barcode

# Save object
saveRDS(seu, 'data/lungs_all/data_lungs_all.rds')