#!/usr/bin/env Rscript

#### Final fine cell type annotation of integrated Seurat object
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)

# Read in object
seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Read in cell type annotations from first round of subsets
myeloid <- read.csv('data/lungs_all/data_lungs_all_myeloid_celltype.csv')
epithelial <- read.csv('data/lungs_all/data_lungs_all_epithelial_celltype.csv')
tcells <- read.csv('data/lungs_all/data_lungs_all_tcells_celltype.csv')
bcells <- read.csv('data/lungs_all/data_lungs_all_bcell_celltype.csv')
fibros <- read.csv('data/lungs_all/data_lungs_all_fibro_celltype.csv')
endo <- read.csv('data/lungs_all/data_lungs_all_endothelial_celltype.csv')

# Add cell type and cell cycle annotations from subsets and tidy up annotation
celltypes <- rbind.data.frame(myeloid, epithelial, tcells, bcells, fibros, endo)
seu$barcode <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, celltypes, by = 'barcode')
rownames(seu@meta.data) <- seu$barcode
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine) == T, seu$cell_type_main, seu$cell_type_fine)
seu$group <- ifelse(seu$group == 'cov', 'COVID-19', 'Control')
seu$initial_clustering <- seu$integrated_snn_res.0.8

# Remove macrophage clusters 3 & 16
seu <- subset(seu, cell_type_fine != 'remove')
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)

# Save object and run DC analyses
saveRDS(seu, 'data/lungs_all/data_lungs_all.rds')


### Add more granular fine cell subtypes after DC analyses
# Read in object
seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Read in fine cell type annotations for macrophages and tuft cells
macs <- read.csv('data/lungs_all/data_lungs_all_macs_celltype.csv')
tuft <- read.csv('data/lungs_all/data_lungs_all_tuft_celltype.csv')
celltypes2 <- rbind.data.frame(macs, tuft)
seu@meta.data <- left_join(seu@meta.data, celltypes2, by = 'barcode')
rownames(seu@meta.data) <- seu$barcode
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine.y) == T, seu$cell_type_fine.x, seu$cell_type_fine.y)

# Save object
saveRDS(seu, 'data/lungs_all/data_lungs_all.rds')