# Author: Somnath Tagore, Ph.D. Title: Drug Target Identification
# Script Name: Drug_Target_Identification.R Last Updated: 02/05/2021

# Packages required for this analysis

library(pheatmap)

# Put here the variable name of the protein activity
fibroblasts.like.cells.raw.data.ges.vpmat <- readRDS("fibroblasts.like.cells.raw.data.ges.vpmat.rds")

# Loading the drug file
drugs <- readRDS("drugs.rds")

dim(drugs)

drugs[1:10, ]
plot.data <- fibroblasts.like.cells.raw.data.ges.vpmat
dim(plot.data)

# There are only the regulators that are druggable
plot.data = plot.data[intersect(as.character(unique(drugs$target)), rownames(plot.data)), 
]
dim(plot.data)
plot.data[1:5, 1:5]

# Load labels for clustering
gene.exp.labels <- read.csv(file = "vip.clust.fibroblasts.overlay.with.gene_exp.final.csv")  #pathogenic vs non-pathogenic

plot.data = plot.data[, match(gene.exp.labels[, 1], colnames(plot.data))]

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  unique(breaks)
}
mat_breaks <- quantile_breaks(plot.data[, 1], n = 100)

# Load labels based on protein-activity based clustering
vip.clus <- read.csv(file = "vip.clust.fibroblasts.csv")
annot.col = data.frame(cluster = gene.exp.labels[, 4])

rownames(annot.col) <- colnames(plot.data)

# Sort cells based on viper clustering
o = order(annot.col$cluster, decreasing = F)

plot.data = plot.data[, o]

# Integrate the scores
integrated_scores = apply(plot.data, 1, function(x) {
  sum(x)/sqrt(length(x))
})
plot.data = plot.data[order(integrated_scores, decreasing = T), ]

# total oncoprotein -log10 bonferroni adjusted p-value
plot.data <- pnorm(q = as.matrix(plot.data), lower.tail = FALSE)
for (i in 1:ncol(plot.data)) {
  padj.res <- plot.data[, i] * nrow(plot.data)
  padj.res[which(padj.res > 1)] = 1
  plot.data[, i] <- padj.res
}
plot.data <- -1 * log(plot.data, base = 10)
mat_breaks <- quantile_breaks(plot.data, n = 100)
mat_breaks = c(0, mat_breaks[which(mat_breaks > 2)])

# Heatmap colors
myColor <- colorRampPalette(c("white", "red"))(length(mat_breaks))
rownames(plot.data) <- getSYMBOL(rownames(plot.data), data = "org.Hs.eg")

# Unique colors
ge.colors.1 <- ClusterColors(length(unique(gene.exp.labels[, 4])))
names(ge.colors) <- unique(gene.exp)
unique(gene.exp)
ge.colors.1 <- list(ge.colors.1)

# Final heatmap
test <- pheatmap(plot.data, color = myColor, fontsize_row = 18, fontsize = 24, 
                 fontsize_col = 8, show_rownames = TRUE, show_colnames = F, cluster_cols = FALSE, 
                 cluster_rows = TRUE, annotation_col = annot.col, annotation_colors = ge.colors.1)

save_pheatmap_pdf(test, "Fibroblasts_OncoTarget_Updated_Updated.pdf")
