#!/usr/bin/env Rscript

### title: Analysis of B cells in lung samples (7 COVID-19 and 1 control samples). Code for Figure2efh, Extended Data Figure 3acdf, Extended Data Table 4.
### authors: Jana Biermann, Yiping Wang
### date: 10/05/20

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)


pat <- 'lungs_all'

#### UMAP analysis (Figure 2e) ####
# Read in file from script #3
patient <- readRDS(paste0('data/', pat, '/data_', pat, '_with_covid_sigs.rds'))
patient@meta.data$group <- ifelse(patient@meta.data$group == 'cov', 'COVID-19', 'Control')

# Subset to B cells as identified in cell_type_main and celltype_bped_fine (SingleR)
patient <- subset(patient, cell_type_main %in% c('B cells'))
patient <- subset(patient, celltype_bped_fine %in% c('Class-switched memory B-cells', 'Memory B-cells', 'naive B-cells', 'Plasma cells'))

# Rerun Seurat workflow
patient <- ScaleData(object = patient)
patient <- RunPCA(object = patient)
patient <- FindNeighbors(patient, dims = 1:20)
patient <- FindClusters(patient, resolution = 1.8)
patient <- RunUMAP(object = patient, dims = 1:20)

# Save UMAP (Figure 2e)
pdf(file = paste0('data/', pat, '/b_cells/Figure_2e.pdf'))
plot <- DimPlot(patient, reduction = 'umap', label = F, group.by = 'celltype_bped_fine', repel = T, cols = viridis(4))
AugmentPlot(plot, dpi = 300)
dev.off()


#### Diffusion component analysis (Figure 2f) ####
es <- as.ExpressionSet(as.data.frame(t(patient@assays$integrated@data)))
es@phenoData@data <- patient@meta.data

# Make a diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Save diffusion map (Figure 2f)
pdf(paste0('data/', pat, '/b_cells/Figure_2f.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(viridis(4))
plot(dm, col = as.factor(es@phenoData@data$celltype_bped_fine), pch = 20)
legend('bottomright', inset = c(-0.07, -0.2), legend = levels(as.factor(es@phenoData@data$celltype_bped_fine)), 
       pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$celltype_bped_fine))), bty = 'n')
dev.off()


#### Extended Data Figure 3acd ####
pdf(file = paste0('data/', pat, '/b_cells/Extended_Data_Figure_3acd.pdf'))
DimPlot(patient, reduction = 'umap', label = F, group.by = 'group', repel = T, cols = viridis(3))
FeaturePlot(patient, features = c('SDC1'), min.cutoff = 'q05', max.cutoff = 'q95')
FeaturePlot(patient, features = c('MS4A1'), min.cutoff = 'q05', max.cutoff = 'q95')
dev.off()


#### Identfication of variable and constant IG genes (Figure 2h, Extended Data Table 4, Extended Data Figure 3f)
# List all IG genes
ig <- grep('^IG', rownames(patient@assays$RNA@counts), value = T)

# Create gene expression dataframe with only IG genes, delete empty rows and columns
tab <- as.data.frame(patient@assays$RNA@data[ig, ])
tab <- tab[, colSums(tab) > 0]
tab <- tab[rowSums(tab) > 0, ]

# List heavy and light variable regions
ighv <- grep('^IGHV', rownames(tab), value = T)
iglv <- grep('^(IGLV|IGKV)', rownames(tab), value = T)

# Set up empty dataframe for results
res <- matrix(data = 0, nrow = length(iglv), ncol = length(ighv))
colnames(res) <- ighv
rownames(res) <- iglv
res <- as.data.frame(res)

# List constant regions, manually remove the ones that don't follow the right pattern, and set up
# result dataframe for constant chains
ig_con <- grep('^(IGHG|IGHA|IGHM|IGHE)', rownames(tab), value = T)
ig_con <- ig_con[-1]
ig_con <- ig_con[-4]
ig_con <- ig_con[-9]
res_con <- as.data.frame(matrix(data = NA, nrow = 500, ncol = 4))
colnames(res_con) <- c('barcode', 'heavy', 'light', 'constant')

# Loop to identify cells with variable and constant chains expressed
i <- 1
for (c in 1:ncol(tab)) {
  # Remove values from previous loop
  rm('igl')
  rm('igh')
  # Go through gene matrix cell by cell
  cell <- c(tab[, c])
  names(cell) <- rownames(tab)
  
  # Find cell's highest expressed heavy chain
  heav <- cell[ighv]
  igh <- names(sort(heav, decreasing = T)[1])
  
  # Find cell's highest expressed light chain
  lig <- cell[iglv]
  igl <- names(sort(lig, decreasing = T)[1])
  
  # Find cell's highest expressed constant region
  con <- cell[ig_con]
  con <- names(sort(con, decreasing = T)[1])
  
  # Only add cell to variable regions result dataframe if it has both heavy and light variable chains expressed
  if (cell[igh] > 0 & cell[igl] > 0) {
    res[igl, igh] <- res[igl, igh] + 1
  }
  
  # Only add cell to constant regions result dataframe if it has heavy and light variable chains and a constant chain expressed
  if (cell[igh] > 0 & cell[igl] > 0 & cell[con] > 0) {
    res_con$barcode[i] <- colnames(tab)[c]
    res_con$heavy[i] <- igh
    res_con$light[i] <- igl
    res_con$constant[i] <- con
    i <- i + 1
  }
}

# Remove columns and rows that are empty
res <- res[, colSums(res) > 0]
res <- res[rowSums(res) > 0, ]

# Add IG subtype info
res_con$subtype <- substr(res_con$constant, 1, 4)
write.csv(res_con, paste0('data/', pat, '/b_cells/Extended_Data_Table_4.csv'), row.names = F)

# Save Figure 2h
pdf(paste0('data/', pat, '/b_cells/Figure_2h.pdf'))
heatmap.2(as.matrix(res), trace = 'none', scale = 'none', col = rev(heat.colors(100)), margins = c(6, 6), cexRow = 0.7, 
          cexCol = 0.7, hclustfun = function(d) {hclust(d, method = 'average')}, xlab = 'Heavy chains', ylab = 'Light chains', 
          lhei = c(1.7, 5.3))
dev.off()


# Run Code to make Extended Data Figure 3f

source("15_IG_combination_heatmaps.R")

##### Save B cell image #####
save.image(paste0('data/', pat, '/b_cells/image_b_cells.RData'))


print(Sys.time())
