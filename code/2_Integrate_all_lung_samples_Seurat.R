#!/usr/bin/env Rscript

### title: Seurat integration for lung samples (7 COVID-19 and 1 control donors; input for diffusion component analysis)
### author: Jana Biermann
### date: 09/29/20

# Libraries
library(dplyr)
library(Seurat)

# Patients to be included
pat_list <- c('lungs_02_cov', 'lungs_03_cov', 'lungs_05_cov', 'lungs_06_cov', 'lungs_08_cov', 'lungs_11_cov', 'lungs_12_cov', 'lungs_14_ctr')

# Load objects
for (pat in pat_list) {
  assign(pat, readRDS(paste0('data/', pat, '/data_', pat, '_cb.rds')))
}

# Create object.list for integration
object.list <- NULL
for (obj in grep('lungs_', ls(), value = T)) {
  object.list <- c(object.list, eval(parse(text = obj)))
}

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:20)

# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:20)

# Seurat workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu)
seu <- RunUMAP(object = seu, dims = 1:20)

# Cell cycle scoring
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'

# Save object
ifelse(!dir.exists(file.path('data/lungs_all/')), dir.create(file.path('data/lungs_all/')), FALSE)
saveRDS(seu, file = 'data/lungs_all/data_lungs_all.rds')

print(Sys.time())