# Author: Somnath Tagore, Ph.D. Title: Master Regulator Analysis of
# COVID data Script Name: Master_Regulator_Analysis_COVID.R Last
# Updated: 02/05/2021

# Packages required for this analysis
formatR::tidy_app()
install.packages(c("cluster", "ggplot2", "devtools", "dplyr", "XML", "Seurat", 
                   "pheatmap", "BiocManager", "RColorBrewer", "fpc", "org.Hs.eg.db"))
BiocManager::install("viper")
BiocManager::install("biomaRt")
BiocManager::install("annotate")
devtools::install_github("JEFworks/MUDAN")

library(cluster)
library(ggplot2)
library(viper)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(fpc)
library(ggrepel)
source("Functions_MasterRegulators.R")

# Step 1: Generate sub-matrices from original gene expression data
# based on pre-defined sub-clusters. Load the Seurat Object
seu.data <- readRDS(file = "data_lungs_all_with_cell_type_intermediate.rds")
print(class(seu.data))

# Display the Object
print(dim(seu.data@assays$RNA@counts))
print(table(seu.data$cell_type_main))
print("---------------------------------")
print(table(seu.data$cell_type_fine))
print(table(seu.data$cell_type_macs))

# Extract the sub-clusters
fibroblasts.Adventitial.FB.raw.data <- subset(seu.data, cell_type_fine == 
                                                "Adventitial FB")@assays$RNA@counts
print("fibroblasts.Adventitial.FB.raw.data")
print(dim(fibroblasts.Adventitial.FB.raw.data))
saveRDS(fibroblasts.Adventitial.FB.raw.data, file = "fibroblasts.Adventitial.FB.raw.data.rds")
write.csv(colnames(fibroblasts.Adventitial.FB.raw.data), file = "fibroblasts.Adventitial.FB.raw.data.csv")

fibroblasts.Alveolar.FB.raw.data <- subset(seu.data, cell_type_fine == 
                                             "Alveolar FB")@assays$RNA@counts
print("fibroblasts.Alveolar.FB.raw.data")
print(dim(fibroblasts.Alveolar.FB.raw.data))
saveRDS(fibroblasts.Alveolar.FB.raw.data, file = "fibroblasts.Alveolar.FB.raw.data.rds")
write.csv(colnames(fibroblasts.Alveolar.FB.raw.data), file = "fibroblasts.Alveolar.FB.raw.data.csv")

fibroblasts.Mesothelial.FB.raw.data <- subset(seu.data, cell_type_fine == 
                                                "Mesothelial FB")@assays$RNA@counts
print("fibroblasts.Mesothelial.FB.raw.data")
print(dim(fibroblasts.Mesothelial.FB.raw.data))
saveRDS(fibroblasts.Mesothelial.FB.raw.data, file = "fibroblasts.Mesothelial.FB.raw.data.rds")
write.csv(colnames(fibroblasts.Mesothelial.FB.raw.data), file = "fibroblasts.Mesothelial.FB.raw.data.csv")

fibroblasts.Other.FB.raw.data <- subset(seu.data, cell_type_fine == "Other FB")@assays$RNA@counts
print("fibroblasts.Other.FB.raw.data")
print(dim(fibroblasts.Other.FB.raw.data))
saveRDS(fibroblasts.Other.FB.raw.data, file = "fibroblasts.Other.FB.raw.data.rds")
write.csv(colnames(fibroblasts.Other.FB.raw.data), file = "fibroblasts.Other.FB.raw.data.csv")

fibroblasts.Airway.smooth.muscle.raw.data <- subset(seu.data, cell_type_fine == 
                                                      "Airway smooth muscle")@assays$RNA@counts
print("fibroblasts.Airway.smooth.muscle.raw.data")
print(dim(fibroblasts.Airway.smooth.muscle.raw.data))
saveRDS(fibroblasts.Airway.smooth.muscle.raw.data, file = "fibroblasts.Airway.smooth.muscle.raw.data.rds")
write.csv(colnames(fibroblasts.Airway.smooth.muscle.raw.data), file = "fibroblasts.Airway.smooth.muscle.raw.data.csv")

fibroblasts.Pericyte.raw.data <- subset(seu.data, cell_type_fine == "Pericyte")@assays$RNA@counts
print("fibroblasts.Pericyte.raw.data")
print(dim(fibroblasts.Pericyte.raw.data))
saveRDS(fibroblasts.Pericyte.raw.data, file = "fibroblasts.Pericyte.raw.data.rds")
write.csv(colnames(fibroblasts.Pericyte.raw.data), file = "fibroblasts.Pericyte.raw.data.csv")

fibroblasts.Vascular.smooth.muscle.raw.data <- subset(seu.data, cell_type_fine == 
                                                        "Vascular smooth muscle")@assays$RNA@counts
print("fibroblasts.Vascular.smooth.muscle.raw.data")
print(dim(fibroblasts.Vascular.smooth.muscle.raw.data))
saveRDS(fibroblasts.Vascular.smooth.muscle.raw.data, file = "fibroblasts.Vascular.smooth.muscle.raw.data.rds")
write.csv(colnames(fibroblasts.Vascular.smooth.muscle.raw.data), file = "fibroblasts.Vascular.smooth.muscle.raw.data.csv")

fibroblasts.Pathological.FB.raw.data <- subset(seu.data, cell_type_fine == 
                                                 "Pathological FB")@assays$RNA@counts
print("fibroblasts.Pathological.FB.raw.data")
print(dim(fibroblasts.Pathological.FB.raw.data))
saveRDS(fibroblasts.Pathological.FB.raw.data, file = "fibroblasts.Pathological.FB.raw.data.rds")
write.csv(colnames(fibroblasts.Pathological.FB.raw.data), file = "fibroblasts.Pathological.FB.raw.data.csv")

# Step 2: Run ARACNe using the above sub-matrices
# https://github.com/califano-lab/ARACNe-AP

# Step 2a: Run ARACNe on fibroblasts.Pathological.FB.raw.data.rds
fibroblasts.Pathological.FB.raw.data <- readRDS(file = "fibroblasts.Pathological.FB.raw.data.rds")
dim(fibroblasts.Pathological.FB.raw.data)
fibroblasts.Pathological.FB.raw.data[1:10, 1:5]

rownames(fibroblasts.Pathological.FB.raw.data) <- mapIds(org.Hs.eg.db, 
                                                         rownames(fibroblasts.Pathological.FB.raw.data), "ENTREZID", "SYMBOL")
fibroblasts.Pathological.FB.raw.data[1:20, 1:5]

fibroblasts.Pathological.FB.raw.data <- fibroblasts.Pathological.FB.raw.data[!is.na(rownames(fibroblasts.Pathological.FB.raw.data)), 
]
dim(fibroblasts.Pathological.FB.raw.data)

# log2TPM+1 normalization
fibroblasts.Pathological.FB.raw.data.ges <- t(t(fibroblasts.Pathological.FB.raw.data)/(colSums(fibroblasts.Pathological.FB.raw.data)/1e+06))
fibroblasts.Pathological.FB.raw.data.ges <- log2(fibroblasts.Pathological.FB.raw.data.ges + 
                                                   1)

# sub-sampling for 2000 cells
fibroblasts.Pathological.FB.raw.data.ges <- fibroblasts.Pathological.FB.raw.data.ges[, 
                                                                                     sample(colnames(fibroblasts.Pathological.FB.raw.data.ges), min(ncol(fibroblasts.Pathological.FB.raw.data.ges), 
                                                                                                                                                    2000))]
fibroblasts.Pathological.FB.raw.data.ges <- fibroblasts.Pathological.FB.raw.data.ges[which(rowSums(fibroblasts.Pathological.FB.raw.data.ges) > 
                                                                                             0), ]

# save data as .txt file
write.table(as.matrix(fibroblasts.Pathological.FB.raw.data.ges), file = "fibroblasts.Pathological.FB.raw.data.ges.txt", 
            sep = "\t", col.names = NA, quote = FALSE)

# Step 3 & 4: Generate networks specific to sub-matrices based on Step
# 2 Step 3
fibroblasts.Pathological.FB.raw.data.ges.interactome <- aracne2regulon("network.fibroblasts.Pathological.FB.raw.data.ges.txt", 
                                                                       as.matrix("fibroblasts.Pathological.FB.raw.data.ges.txt"), gene = FALSE, 
                                                                       verbose = TRUE)
saveRDS(fibroblasts.Pathological.FB.raw.data.ges.interactome, file = "fibroblasts.Pathological.FB.raw.data.ges.interactome.rds")

# Step 4a Perform the steps 2a to 3 for all sub-matrices

# Step 4b Integrate the networks based on significant p values.

# Step 5: Protein activity prediction using VIPER. viper without null

# Prune Regulon
pregul.fibroblasts.Pathological.FB.raw.data.ges.interactome <- pruneRegulon(fibroblasts.Pathological.FB.raw.data.ges.interactome, 
                                                                            cutoff = 50)

# Run VIPER
fibroblasts.Pathological.FB.raw.data.ges.vpmat <- viper(eset = as.matrix(fibroblasts.Pathological.FB.raw.data.ges), 
                                                        regulon = pregul.fibroblasts.Pathological.FB.raw.data.ges.interactome, 
                                                        method = "none")

saveRDS(fibroblasts.Pathological.FB.raw.data.ges.vpmat, file = "fibroblasts.Pathological.FB.raw.data.ges.vpmat.rds")

# Step 6: Master regulator analysis.

# load the protein activity matrix
fibroblasts.like.cells.raw.data.ges.vpmat <- readRDS(file = "fibroblasts.Pathological.FB.raw.data.ges.vpmat.rds")

# copy the protein activity into an object
pbmc.vip <- fibroblasts.like.cells.raw.data.ges.vpmat

# Generate a distance matrix and find the optimal clustering using
# silhouette score optimized Louvain clsutering:
pbmc.vip.dist <- as.dist(1 - as.matrix(viperSimilarity(fibroblasts.like.cells.raw.data.ges.vpmat)))
saveRDS(pbmc.vip.dist, file = "pbmc.vip.dist.rds")

# Find the optimal clustering using silhouette score optimized Louvain
# clustering
vip.clust <- LouvainClustering(pbmc.vip, dist.mat = pbmc.vip.dist, rmax = 350, 
                               rstep = 25)
saveRDS(vip.clust, file = "vip.clust.rds")

# Sort the clusters
opt.clust <- sort(vip.clust$clusterings[[which.max(vip.clust$sils)]])
saveRDS(opt.clust, file = "opt.clust.rds")

# Identify the Master Regulators (the proteins that are most
# deferentially active between clusters) using a Mann Whitney U-Test.
vip.mrs <- MasterRegulators(pbmc.vip, opt.clust)
saveRDS(vip.mrs, file = "vip.mrs.rds")

# Select the number of MRs that need to be displayed
num.mrs <- 10

# Order the cells based on optimal clusters
cell.order <- names(opt.clust)
saveRDS(cell.order, file = "cell.order.rds")

# Select the master regulator set
mr.set <- sapply(vip.mrs, function(x) {
  names(x$positive[1:num.mrs])
})
mr.set <- unique(unlist(as.list(mr.set)))
saveRDS(mr.set, file = "mr.set.rds")

# Cell type annotations gene expression labels
plot.mat.no.sort <- read.csv(file = "vip.clust.fibroblasts.overlay.with.gene_exp.csv")
plot.mat.no.sort <- plot.mat.no.sort[, 1]

# gene expression labels based on Pathological vs Non-pathological
# Fibroblasts
plot.mat.no.sort <- read.csv(file = "vip.clust.fibroblasts.overlay.with.gene_exp.int.csv")  #pathogenic vs non-pathogenic
# plot.mat.no.sort<- plot.mat.no.sort[,3]

plot.mat.no.sort <- plot.mat.no.sort[, 1]

# Sort by protein activity clusters
plot.mat <- pbmc.vip[mr.set, cell.order]  #sort by protein activity clusters
plot.mat <- pbmc.vip[mr.set, match(plot.mat.no.sort, colnames(pbmc.vip))]  # sort by pathological vs non-pathological

## Set colors and annotations
any(is.na(plot.mat))
plot.mat <- plot.mat[, is.finite(colSums(plot.mat))]

mat.breaks <- QuantileBreaks(as.matrix(plot.mat))

# Colors based on Cell types Crude classification by gene expression
# labels
gene.exp.labels <- read.csv(file = "vip.clust.fibroblasts.overlay.with.gene_exp.csv")

## Classification by gene expression labels Pathological FB vs Others
gene.exp.labels <- read.csv(file = "vip.clust.fibroblasts.overlay.with.gene_exp.int.csv")  #pathogenic vs non-pathogenic
gene.exp <- gene.exp.labels[, 5]

# Define color coding
clust.colors <- ClusterColors(length(unique(opt.clust)))
names(clust.colors) <- unique(opt.clust)
ge.colors <- ClusterColors(length(unique(gene.exp)))
names(ge.colors) <- unique(gene.exp)

# Color codes for Fibroblasts
ge.colors <- c(`Vascular smooth muscle` = "#D39200", `Intermediate pathological FB` = "#00BA38", 
               `Pathological FB` = "#F8766D", Pericyte = "#93AA00", `Airway smooth muscle` = "#00BA38", 
               `Adventitial FB` = "#00C19F", `Alveolar FB` = "#00B9E3", `Mesothelial FB` = "#619CFF", 
               `Other FB` = "#FF61C3")

# Color codes for Pathological vs Non-pathological Fibroblasts
ge.colors <- c(`Non-Pathological FB` = "#00BFC4", `Intermediate pathological FB` = "#00BA38", 
               `Pathological FB` = "#F8766D")
unique(gene.exp)

# Define the annotation dataframe
annot.df.1 <- data.frame(opt.clust[cell.order])
annot.df.1 <- data.frame(colnames(plot.mat))

# Just for Gene Expression based annotation
annot.df.1[, 2] <- gene.exp

# Consider this if both protein activity and gene expression based
# annotation need to be considered
annot.df.1[, 3] <- opt.clust[cell.order]

annot.df <- as.data.frame(annot.df.1[, 2:3])
annot.df <- as.data.frame(annot.df.1[, 2])
rownames(annot.df) <- colnames(plot.mat)
colnames(annot.df) <- c("Gene_Exp")  #, 'Prot_Act')
colnames(annot.df) <- c("Gene_Exp", "Prot_Act")

# fibroblasts
annot.df$Gene_Exp <- factor(annot.df$Gene_Exp, levels = c("Adventitial FB", 
                                                          "Alveolar FB", "Mesothelial FB", "Other FB1", "Other FB2", "Other FB3", 
                                                          "Pathological FB"))

annot.df$Gene_Exp <- factor(annot.df$Gene_Exp, levels = c("Pathological FB", 
                                                          "Vascular smooth muscle", "Pericyte", "Airway smooth muscle", "Adventitial FB", 
                                                          "Alveolar FB", "Mesothelial FB", "Intermediate pathological FB", "Other FB"))

# patho vs non-patho
annot.df$Gene_Exp <- factor(annot.df$Gene_Exp, levels = c("Pathological FB", 
                                                          "Non-Pathological FB", "Intermediate pathological FB"))

annot.color <- list(Gene_Exp = ge.colors)  #,

annot.color <- list(Gene_Exp = ge.colors, Prot_Act = clust.colors)

rownames(plot.mat) <- getSYMBOL(rownames(plot.mat), data = "org.Hs.eg")
test1 <- pheatmap(plot.mat, main = "Viper Clustering: Master Regulators\nFibroblasts", 
                  fontsize = 20, annotation_col = annot.df, annotation_colors = annot.color, 
                  cluster_cols = FALSE, show_colnames = FALSE, cluster_rows = FALSE, 
                  show_rownames = TRUE, fontsize_row = 7, fontsize_col = 12, breaks = mat.breaks, 
                  color = ColorLevels(length(mat.breaks) - 1, "vip"))
save_pheatmap_pdf(test1, "VIPER_MRs_Fibroblasts_Like_Updated_Updated_2.pdf")