#!/usr/bin/env Rscript

#### Fine cell type annotation of B cell subset
library(dplyr)
library(Seurat)

seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Subset to B cells
seu <- subset(seu, cell_type_main == 'B cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:20)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/lungs_all/markers_lungs_all_bcells.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(5), 'B cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(13, 4, 9, 8), 'Activated B cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 %in% c(10, 11, 12, 0, 2, 6, 7, 3, 1), 'Plasma cells', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/lungs_all/data_lungs_all_bcells.rds')
celltype <- seu@meta.data %>% select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/lungs_all/data_lungs_all_bcell_celltype.csv')