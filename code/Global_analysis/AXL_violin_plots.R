#!/usr/bin/env Rscript

### title: Print out violin plots of AXL expression across cell_type_main,
### intermediate and fine classifications author: Yiping Wang date: 02/08/2021

# Set Identities in data for all samples, using either cell_type_main,
# intermediate, or fine then print out violin plots of AXL expression
uniquecelltypesmain = unique(data_lungs_all$cell_type_main)
for (i in 1:length(uniquecelltypesmain)) {
    data_lungs_all = SetIdent(data_lungs_all, cells = colnames(data_lungs_all)[data_lungs_all$cell_type_main == 
        uniquecelltypesmain[i]], value = uniquecelltypesmain[i])
}
pdf("Extended Data Myeloid 2 c.pdf")
aplot = VlnPlot(data_lungs_all, features = c("AXL"), log = T)
AugmentPlot(aplot, dpi = 300)
print(aplot)
dev.off()

uniquecelltypesintermediate = unique(data_lungs_all$cell_type_intermediate)
for (i in 1:length(uniquecelltypesintermediate)) {
    data_lungs_all = SetIdent(data_lungs_all, cells = colnames(data_lungs_all)[data_lungs_all$cell_type_intermediate == 
        uniquecelltypesintermediate[i]], value = uniquecelltypesintermediate[i])
}
pdf("AXL_allcelltypes_intermediate.pdf")
aplot = VlnPlot(data_lungs_all, features = c("AXL"), log = T)
AugmentPlot(aplot, dpi = 300)
print(aplot)
dev.off()

uniquecelltypesfine = unique(data_lungs_all$cell_type_fine)
for (i in 1:length(uniquecelltypesfine)) {
    data_lungs_all = SetIdent(data_lungs_all, cells = colnames(data_lungs_all)[data_lungs_all$cell_type_fine == 
        uniquecelltypesfine[i]], value = uniquecelltypesfine[i])
}
pdf("AXL_allcelltypes_fine.pdf", height = 7, width = 21)
aplot = VlnPlot(data_lungs_all, features = c("AXL"), log = T)
AugmentPlot(aplot, dpi = 300)
print(aplot)
dev.off()
