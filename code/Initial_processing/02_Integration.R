#!/usr/bin/env Rscript

#### Seurat intergation of individual objects (RPCA)
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(future)
library(future.apply)


# Set up future for parallelization
plan('multiprocess', workers = 40)
options(future.globals.maxSize = 2e+05 * 1024^2)

# Patients that will be integrated in final object
pat_list <- c('C51ctr', 'C52ctr', 'C53ctr', 'C54ctr', 'C55ctr', 'C56ctr', 'C57ctr', 'L01cov', 'L03cov', 
              'L04cov', 'L04covaddon', 'L05cov', 'L06cov', 'L07cov', 'L08cov', 'L09cov', 'L10cov', 'L11cov',
              'L12cov', 'L13cov', 'L15cov', 'L16cov', 'L17cov', 'L18cov', 'L19cov', 'L21cov', 'L22cov')

# Load objects
for (pat in pat_list) {
  assign(pat, readRDS(paste0('data/', pat, '/data_', pat, '_cb.rds')))
}

# Prepare object.list for integration
object.list <- NULL
for (obj in grep(paste(c('cov', 'ctr'), collapse = '|'), ls(), value = T)) {
  object.list <- c(object.list, eval(parse(text = obj)))
}

object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = object.list)
object.list <- future_lapply(X = object.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors using references (RPCA)
anchors <- FindIntegrationAnchors(object.list = object.list, reference = c(1, 2, 3, 9, 20, 23), 
                                  reduction = 'rpca', dims = 1:50)

# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:50)

# Normal workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu)
seu <- RunUMAP(object = seu, dims = 1:30)

# Save object
ifelse(!dir.exists(file.path('data/lungs_all/')), dir.create(file.path('data/lungs_all/')), FALSE)
saveRDS(seu, file = 'data/lungs_all/data_lungs_all_initial_obj.rds')