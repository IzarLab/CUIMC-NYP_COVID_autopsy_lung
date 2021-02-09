#!/usr/bin/env Rscript

#### Seurat analysis for cellbender output and doublet rate provided as an argument
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(DropletUtils)
source('misc/ReadCB_h5.R')


# Set patient and doublet rate variables
pat <- commandArgs()[6]
doublet_rate <- as.numeric(commandArgs()[7])

# Load data set
seu.data <- ReadCB_h5(paste0('data/', pat, '/', pat, '_filtered.h5'))

# Initialize the Seurat object with the raw (non-normalized) data
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

# Annotate MT, RPS, and RPL genes
seu_raw[['percent.mt']] <- PercentageFeatureSet(seu_raw, pattern = '^MT-')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPS', col.name = 'percent.rps')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPL', col.name = 'percent.rpl')
seu_raw$percent.rp <- seu_raw$percent.rps + seu_raw$percent.rpl

# Annotate patient info
seu_raw[['patient']] <- pat
seu_raw[['organ']] <- 'Lung'
seu_raw[['group']] <- substr(pat, 4, 6)
seu_raw[['ID']] <- substr(pat, 1, 3)

# Identify doublets using scrublet
writeMM(seu_raw@assays$RNA@counts, paste0('data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste('python3 misc/scrublet_code.py', pat, doublet_rate))
doublets <- read.table(paste0('data/', pat, '/doublets_', pat, '_raw.txt'), header = T)
seu_raw[['predicted_doublets']] <- doublets$predicted_doublets
seu_raw[['doublet_scores']] <- doublets$doublet_scores
system(paste0('rm data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste0('rm data/', pat, '/doublets_', pat, '_raw.txt'))

# Subset
minFeature <- 200
maxFeature <- 7500
minCount <- 400
maxCount <- 40000
maxMT <- 10
seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
                nCount_RNA < maxCount & percent.mt < maxMT & predicted_doublets == F)

# Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, reduction = 'pca', dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:30)

# Save object
saveRDS(seu, file = paste0('data/', pat, '/data_', pat, '_cb.rds'))