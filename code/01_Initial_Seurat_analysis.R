#!/usr/bin/env Rscript

### title: Seurat analysis for CellBender output with sample name and expected doublet rate provided in arguments
### author: Jana Biermann
### date: 09/29/20

# Run for each patient individually

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(purrr)
library(rscrublet)
library(DropletUtils)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)


pat <- commandArgs()[6]  # sample name
doubs <- as.numeric(commandArgs()[7])  # expected doublet rate (0.096 or 0.04)
print(pat)
print(doubs)

#### Loading, QC, dimension reduction ####
# Load dataset
seu.data <- Read10X_h5(paste0('data/', pat, '/', pat, '_filtered.h5'), use.names = TRUE, unique.features = TRUE)

# Initialize the Seurat object
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

# Annotate MT genes
seu_raw[['percent.mt']] <- PercentageFeatureSet(seu_raw, pattern = '^MT-')

# Annotate patient info
seu_raw[['patient']] <- pat
seu_raw[['organ']] <- strsplit(pat, '_')[[1]][1]
seu_raw[['group']] <- strsplit(pat, '_')[[1]][3]
seu_raw[['ID']] <- strsplit(pat, '_')[[1]][2]
seu_raw[['sequencing']] <- 'sn'

# Identify doublets using scrublet
scrub = scrubDoublets(as.matrix(seu_raw@assays$RNA@counts), expected_doublet_rate = doubs)
seu_raw[['ScrubDoublet']] <- scrub$scrubDoublets
seu_raw[['ScrubDoublet_score']] <- scrub$doublet_scores_obs

# Subset using cut-offs below
minFeature <- 200
maxFeature <- 7500
minCount <- 100
maxCount <- 40000
maxMT <- 10
seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & ScrubDoublet == F)

# Seurat workflow
seu <- NormalizeData(seu, normalization.method = 'LogNormalize', scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu)
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, reduction = 'pca', dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:20)

#### Cell type identification using SingleR ####
seu_sce <- as.SingleCellExperiment(seu)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels

iced <- DatabaseImmuneCellExpressionData()
pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
pruneScores(pred_iced_main)
seu[['celltype_iced_main']] <- pred_iced_main$pruned.labels
pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
pruneScores(pred_iced_fine)
seu[['celltype_iced_fine']] <- pred_iced_fine$pruned.labels

hpca <- HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']] <- pred_hpca_main$pruned.labels
pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']] <- pred_hpca_fine$pruned.labels

mid <- MonacoImmuneData()
pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
pruneScores(pred_mid_main)
seu[['celltype_mid_main']] <- pred_mid_main$pruned.labels
pred_mid_fine <- SingleR(test = seu_sce, ref = mid, labels = mid$label.fine)
pruneScores(pred_mid_fine)
seu[['celltype_mid_fine']] <- pred_mid_fine$pruned.labels


#### Save object ####
saveRDS(seu, file = paste0('data/', pat, '/data_', pat, '_cb.rds'))

print(paste('End:', Sys.time()))
