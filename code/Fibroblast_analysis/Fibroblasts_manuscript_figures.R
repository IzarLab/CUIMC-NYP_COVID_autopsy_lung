#!/usr/bin/env Rscript



#### script to generate the figures used in the
#### manuscript
library(Seurat)
library(ggplot2)
library(RColorBrewer)


# *************************************************************
# figures
seu <- readRDS("data/lungs_all/data_lungs_all_fibro.rds")
DimPlot(seu)
DimPlot(seu, group.by = "cell_type_fine")

# creating pathological and non-pathological labels
seu$pathological <- ifelse(seu$cell_type_fine %in% 
                             c("Pathological FB", "Intermediate pathological FB"), 
                           "Pathological", "non-Pathological")
DimPlot(seu, group.by = "pathological")

# computing pvalues
seu <- SetIdent(seu, value = "pathological")
FindMarkers(seu, features = c("CTHRC1", "COL1A1", "COL3A1"), 
            ident.1 = "Pathological")
#        p_val avg_logFC pct.1 pct.2 p_val_adj
# COL1A1     0 2.0490484 0.700 0.286         0
# COL3A1     0 1.6685641 0.811 0.429         0
# CTHRC1     0 0.9039856 0.237 0.048         0

# set the cell types as the active.ident
seu <- SetIdent(seu, value = "cell_type_fine")

# UMAPS color blind color scheme 
# c('#FAE6BE','#A0FA82', '#F0F032', '#0AB45A', '#FA7850',
# '#14D2DC', '#FA78FA', '#00A0FA', '#005AC8', '#8214A0')

print(DimPlot(seu, reduction = "umap", label = T, group.by = "group", 
              shuffle = T, raster = T, cols = c("#006E82", "#AA0A3C")))
# ggsave('manuscript/figs/umap_lung_fb_covCtrl_label.pdf')
print(DimPlot(seu, reduction = "umap", label = F, group.by = "group", 
              shuffle = T, raster = T, cols = c("#006E82", "#AA0A3C")))
# ggsave('manuscript/figs/umap_lung_fb_covCtrl_nolabel.pdf')

print(DimPlot(seu, reduction = "umap", label = T, group.by = "cell_type_fine", 
              shuffle = T, raster = T, cols = c("#8214A0", "#A0FA82", 
                                                "#00A0FA", "#0AB45A", "#FA78FA", "#14D2DC", 
                                                "#FA7850", "#00A0FA", "#F0F032")))
# ggsave('manuscript/figs/umap_lung_fb_cell_type_fine_label.pdf')
print(DimPlot(seu, reduction = "umap", label = F, group.by = "cell_type_fine", 
              shuffle = T, raster = T, cols = c("#8214A0", "#A0FA82", 
                                                "#00A0FA", "#0AB45A", "#FA78FA", "#14D2DC", 
                                                "#FA7850", "#00A0FA", "#F0F032")))
# ggsave('manuscript/figs/umap_lung_fb_cell_type_fine_nolabel.pdf',
# width = 8, height = 6)

print(DimPlot(seu, reduction = "umap", label = T, group.by = "pathological", 
              shuffle = T, raster = T, cols = c("#00A0FA", "#FA7850")))
# ggsave('manuscript/figs/umap_lung_fb_pathological_label.pdf')
print(DimPlot(seu, reduction = "umap", label = F, group.by = "pathological", 
              shuffle = T, raster = T, cols = c("#00A0FA", "#FA7850")))
# ggsave('manuscript/figs/umap_lung_fb_pathological_nolabel.pdf')


# DOTPLOTS pathological markers
DotPlot(seu, assay = "RNA", features = unlist(strsplit("COL1A1 COL3A1 LUM POSTN CTHRC1", 
                                                       split = " ")), cluster.idents = T, group.by = "cell_type_fine", 
        dot.scale = 7) + scale_color_viridis_c() + RotatedAxis()
# ggsave('manuscript/figs/dotplot_lung_fb_pathological.pdf',
# width = 6, height = 5)

# fibroblast subtype markers
temp <- factor(seu@meta.data$cell_type_fine, levels = c("Alveolar FB", 
                                                        "Adventitial FB", "Pathological FB", "Intermediate pathological FB", 
                                                        "Pericyte", "Airway smooth muscle", "Vascular smooth muscle", 
                                                        "Mesothelial FB", "Other FB"))
levels(temp)
seu@meta.data$cell_type_fine_ordered <- temp
DotPlot(seu, assay = "RNA", features = unlist(strsplit("LIMCH1 MACF1 GPC3 ROBO2 ABCA10 SCARA5 LEPR C3 COL1A1 LUM CTHRC1 ADAM12 POSTN CYP7B1 NCAM2 PDGFRB LAMC3 TRPC6 MYH11 ACTG2 ACTA2 RGS5 VWF TINAGL1 PKHD1L1 UPK3B CLDN1 ITGBL1 NALCN HECW1", 
                                                       split = " ")), cluster.idents = F, group.by = "cell_type_fine_ordered", 
        dot.scale = 7) + scale_color_viridis_c() + RotatedAxis()
# ggsave("manuscript/figs/dotplot_lung_fb_subtypeMarkers.pdf", 
       # width = 15, height = 5)

# fibroblasts in all cell types
seu_all <- readRDS("data/data_lungs_all_v3.rds")
DimPlot(seu_all, reduction = "umap", label = T)
DefaultAssay(seu_all) <- "RNA"
FeaturePlot(seu_all, features = c("ACTA2", "LAMA2", 
                                  "COL1A2", "COL1A1"), max.cutoff = "q95", order = T, 
            raster = T)  #+ scale_color_viridis_c()
# ggsave('manuscript/figs/umap_lung_all_fbMarkers.pdf')

# VIOLIN PLOTs pathological violins and pvals
VlnPlot(seu, assay = "RNA", features = c("CTHRC1", 
                                         "COL1A1", "COL3A1"), group.by = "pathological", 
        pt.size = 0, cols = c("#00A0FA", "#FA7850"))
# ggsave("manuscript/figs/violin_lung_fb_pathological_COL1A1_COL3A1_CTHRC1.pdf")  #, width = 15, height = 5)