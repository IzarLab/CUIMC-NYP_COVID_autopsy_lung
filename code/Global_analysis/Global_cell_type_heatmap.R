#!/usr/bin/env Rscript

### title: Generate heatmaps of differential gene expression for top 5 marker genes
### in all cell type clusters, at either cell_type_main or intermediate level
### author: Yiping Wang date: 02/08/2021

filter_every_tenth_gene = FALSE
cell_type_levels_to_display = c("main", "intermediate")

for (i in 1:length(cell_type_level_to_display)) {
    # read in top 5 marker genes for each cell type group
    markers <- read.csv(paste0("markers_cell_type_", cell_type_levels_to_display[i], 
        ".csv"))
    markers <- markers %>% filter(p_val_adj < 0.05)
    topgenes <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% 
        arrange(cluster)
    
    # extract matrix of expression values from seurat data object, select rows that
    # correspond to marker genes
    DefaultAssay(data_lungs_all) = "RNA"
    datamat = data_lungs_all$RNA@scale.data[which(rownames(data_lungs_all) %in% topgenes$gene), 
        ]
    datamat = datamat[match(topgenes$gene, rownames(datamat)), ]
    datamat = datamat[, order(data_lungs_all$cell_type_main)]
    markertype = paste0(topgenes$cluster, " markers")
    celltype = data_lungs_all$cell_type_main
    celltype = celltype[order(data_lungs_all$cell_type_main)]
    
    # take only every cell in data matrix to display, to create smaller pdf file
    if (filter_every_tenth_gene) {
        celltype = celltype[seq(1, dim(datamat)[2], 10)]
        datamat = datamat[, seq(1, dim(datamat)[2], 10)]
    }
    
    annotation_col = data.frame(celltype = factor(celltype))
    rownames(annotation_col) = colnames(datamat)
    
    # some marker genes may be repeated within data matrix for repeat genes, add
    # number at end to distinguish repeats
    labels_row = rownames(datamat)
    annotation_row = data.frame(markertype = factor(markertype))
    uniquegenenames = unique(rownames(datamat))
    for (i in 1:length(uniquegenenames)) {
        if (sum(rownames(datamat) == uniquegenenames[i]) > 1) {
            matchidxs = which(rownames(datamat) == uniquegenenames[i])
            for (j in 1:length(matchidxs)) {
                rownames(datamat)[matchidxs[j]] = paste0(rownames(datamat)[matchidxs[j]], 
                  "_", j)
            }
        }
    }
    rownames(annotation_row) = rownames(datamat)
    
    # clip heatmap values to plus or minus 2.5, to ensure that color contrast is high
    datamat[datamat <= -2.5] = -2.5
    datamat[datamat >= 2.5] = 2.5
    
    pdf(paste0("markers_cell_type_", cell_type_levels_to_display[i], ".pdf"), height = 15, 
        width = 17)
    print(pheatmap(datamat, color = PurpleAndYellow(), annotation_col = annotation_col, 
        annotation_row = annotation_row, labels_row = labels_row, show_colnames = F, 
        cluster_rows = F, cluster_cols = F))
    dev.off()
}
