#!/usr/bin/env Rscript

### title: Create UMAPs of cell type classification and COVID-19 status across
### samples, as well as UMAPs highlighting cells from individual patients author:
### Yiping Wang date: 02/08/2021

consistentcolors = c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", "#005AC8", 
    "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", "#FA7850", 
    "#14D2DC", "#FA78FA")

# first two UMAPs plotted are used to generate Figures 1C and 1D

# UMAP of cell_type_intermediate, and COVID-19 vs. Control status, in one file
pdf(file = paste0("cell_type_intermediate_lungs_all_umap.pdf"), height = 7, width = 14)
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
tplot = DimPlot(data_lungs_all, reduction = "umap", label = T, group.by = "cell_type_intermediate", 
    repel = T, label.size = 2.5, shuffle = T)
AugmentPlot(tplot, dpi = 300)
print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
    override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
    vp = vplayout(1, 1))

tplot2 = DimPlot(data_lungs_all, reduction = "umap", label = F, group.by = "group", 
    repel = T, shuffle = T, cols = consistentcolors[1:2])
AugmentPlot(tplot2, dpi = 300)
print(tplot2 + ggtitle("patientgroup"), vp = vplayout(1, 2))
dev.off()

# same as above, but without cell_type_intermediate group labels, for ease of
# manual editing
pdf(file = paste0("cell_type_intermediate_lungs_all_umap_unlabeled.pdf"), height = 7, 
    width = 14)
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
tplot = DimPlot(data_lungs_all, reduction = "umap", label = F, group.by = "cell_type_intermediate", 
    repel = T, label.size = 2.5, shuffle = T)
# tplot[[1]]$layers[[1]]$aes_params$alpha = .5
AugmentPlot(tplot, dpi = 300)
print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
    override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
    vp = vplayout(1, 1))

tplot2 = DimPlot(data_lungs_all, reduction = "umap", label = F, group.by = "group", 
    repel = T, shuffle = T, cols = consistentcolors[1:2])
# tplot2[[1]]$layers[[1]]$aes_params$alpha = .2
AugmentPlot(tplot2, dpi = 300)
print(tplot2 + ggtitle("patientgroup"), vp = vplayout(1, 2))
dev.off()

# UMAP of cell_type_intermediate groups, in COVID-19 cells only
pdf(file = paste0("cell_type_intermediate_lungs_all_cov_umap.pdf"), height = 7, width = 7)
tplot = DimPlot(subset(data_lungs_all, group == "COVID-19"), reduction = "umap", 
    label = T, group.by = "cell_type_intermediate", repel = T, label.size = 2.5, 
    shuffle = T)
# tplot[[1]]$layers[[1]]$aes_params$alpha = .5
AugmentPlot(tplot, dpi = 300)
print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
    override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
    vp = vplayout(1, 1))
dev.off()

# UMAP of cell_type_intermediate groups, in Control cells only
pdf(file = paste0("cell_type_intermediate_lungs_all_ctr_umap.pdf"), height = 7, width = 7)
tplot = DimPlot(subset(data_lungs_all, group == "Control"), reduction = "umap", label = T, 
    group.by = "cell_type_intermediate", repel = T, label.size = 2.5, shuffle = T)
# tplot[[1]]$layers[[1]]$aes_params$alpha = .5
AugmentPlot(tplot, dpi = 300)
print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
    override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
    vp = vplayout(1, 1))
dev.off()

withlabels = c(TRUE, FALSE)
labelsuffix = c("", "_nolabels")

# highlight either epithelial, immune, or fibroblast cell groups in three
# different UMAP figures, all plotted on the same pdf file create this either
# with or without cell type group labels Immune cells is Figure 2A Epithelial
# cells is Figure 3A Fibroblasts is Figure 4A
for (z in 1:length(withlabels)) {
    DefaultAssay(data_lungs_all) = "RNA"
    pdf(file = paste0("cell_type_intermediate_lungs_all_epithelial_fibroblast_immune_subsets_umap", 
        labelsuffix[z], ".pdf"), height = 7, width = 21)
    pushViewport(viewport(layout = grid.layout(1, 3)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    thisuniquecelltypes = c("Airway epithelial cells", "AT1", "AT2", "Other epithelial cells")
    cellstohighlight = list()
    for (i in 1:length(thisuniquecelltypes)) {
        cellstohighlight = list.append(cellstohighlight, colnames(data_lungs_all)[data_lungs_all$cell_type_intermediate == 
            thisuniquecelltypes[i]])
    }
    sortuniquecelltypes = sort(unique(data_lungs_all$cell_type_intermediate))
    hues = seq(15, 375, length = length(sortuniquecelltypes) + 1)
    gg_color_hue = hcl(h = hues, l = 65, c = 100)[1:length(sortuniquecelltypes)]
    colstohighlight = c()
    for (i in 1:length(thisuniquecelltypes)) {
        colstohighlight = append(colstohighlight, gg_color_hue[match(thisuniquecelltypes[i], 
            sortuniquecelltypes)])
    }
    data_lungs_all_temp = data_lungs_all
    data_lungs_all_temp$cell_type_intermediate[!(data_lungs_all_temp$cell_type_intermediate %in% 
        thisuniquecelltypes)] = "Unselected"
    tplot = DimPlot(data_lungs_all_temp, reduction = "umap", label = withlabels[z], 
        group.by = "cell_type_intermediate", repel = T, label.size = 2.5, shuffle = T) + 
        scale_color_manual(labels = c(thisuniquecelltypes, "Unselected"), values = c(colstohighlight, 
            "gray"))
    AugmentPlot(tplot, dpi = 300)
    print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
        override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
        vp = vplayout(1, 1))
    
    thisuniquecelltypes = c("Fibroblasts", "Smooth muscle")
    colstohighlight = c()
    for (i in 1:length(thisuniquecelltypes)) {
        colstohighlight = append(colstohighlight, gg_color_hue[match(thisuniquecelltypes[i], 
            sortuniquecelltypes)])
    }
    data_lungs_all_temp = data_lungs_all
    data_lungs_all_temp$cell_type_intermediate[!(data_lungs_all_temp$cell_type_intermediate %in% 
        thisuniquecelltypes)] = "Unselected"
    tplot = DimPlot(data_lungs_all_temp, reduction = "umap", label = withlabels[z], 
        group.by = "cell_type_intermediate", repel = T, label.size = 2.5, shuffle = T) + 
        scale_color_manual(labels = c(thisuniquecelltypes, "Unselected"), values = c(colstohighlight, 
            "gray"))
    AugmentPlot(tplot, dpi = 300)
    print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
        override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
        vp = vplayout(1, 2))
    
    thisuniquecelltypes = c("B cells", "CD4+ T cells", "CD8+ T cells", "Cycling NK/T cells", 
        "Dendritic cells", "Macrophages", "Mast cells", "Monocytes", "NK cells", 
        "Plasma cells", "Tregs")
    colstohighlight = c()
    for (i in 1:length(thisuniquecelltypes)) {
        colstohighlight = append(colstohighlight, gg_color_hue[match(thisuniquecelltypes[i], 
            sortuniquecelltypes)])
    }
    data_lungs_all_temp = data_lungs_all
    data_lungs_all_temp$cell_type_intermediate[!(data_lungs_all_temp$cell_type_intermediate %in% 
        thisuniquecelltypes)] = "Unselected"
    tplot = DimPlot(data_lungs_all_temp, reduction = "umap", label = withlabels[z], 
        group.by = "cell_type_intermediate", repel = T, label.size = 2.5, shuffle = T) + 
        scale_color_manual(labels = c(thisuniquecelltypes, "Unselected"), values = c(colstohighlight, 
            "gray"))
    AugmentPlot(tplot, dpi = 300)
    print(tplot + ggtitle("cell_type_intermediate") + guides(col = guide_legend(nrow = 30, 
        override.aes = list(size = 5))) + theme(legend.text = element_text(size = 7)), 
        vp = vplayout(1, 3))
    dev.off()
}

# create series of UMAPs, where cells from one patient at a time are highlighted
# in blue
uniqueidents = unique(data_lungs_all$orig.ident)
pdf("cell_type_intermediate_lungs_all_individual_plots.pdf")
for (i in 1:length(uniqueidents)) {
    statusname = paste0(uniqueidents[i], "_status")
    statusnamebare = uniqueidents[i]
    eval(parse(text = paste0("data_lungs_all$", statusname, " = 0")))
    eval(parse(text = paste0("data_lungs_all$", statusname, "[data_lungs_all$orig.ident==uniqueidents[i]] = 1")))
    eval(parse(text = paste0("data_lungs_all$", statusnamebare, " = 0")))
    eval(parse(text = paste0("data_lungs_all$", statusnamebare, "[data_lungs_all$orig.ident==uniqueidents[i]] = 1")))
}
for (i in 1:floor(length(uniqueidents)/4)) {
    startidx = (i - 1) * 4 + 1
    endidx = min(i * 4, length(uniqueidents))
    aplot = FeaturePlot(data_lungs_all, features = uniqueidents[startidx:endidx], 
        min.cutoff = 0.5, max.cutoff = 1.5)
    AugmentPlot(aplot, dpi = 300)
    print(aplot)
}
dev.off()
