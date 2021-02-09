#!/usr/bin/env Rscript

### title: Add immune_status and cell_type_intermediate classification to seurat
### data object, print out tables to check intersection of cell type
### classifications among three different levels author: Yiping Wang date:
### 02/09/2021

library("rlist")

# By default, cell_type_intermediate is same as cell_type_fine, unless different
# grouping information is required
data_lungs_all$cell_type_intermediate = data_lungs_all$cell_type_fine
# By default, cells are assumed to be non-immune cells
data_lungs_all$immune_status = "Non-immune"
allcelltypesfine = sort(unique(data_lungs_all$cell_type_fine))

for (i in 1:length(allcelltypesfine)) {
    # If cell_type_fine is an immune cell type, change immune_status to immune
    if (str_detect(allcelltypesfine[i], "B cells") || allcelltypesfine[i] == "Plasma cells") {
        data_lungs_all$immune_status[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Immune"
    }
    if (str_detect(allcelltypesfine[i], "myeloid") || allcelltypesfine[i] == "Monocytes" || 
        allcelltypesfine[i] == "Macrophages" || allcelltypesfine[i] == "Dendritic cells" || 
        allcelltypesfine[i] == "Transitioning MDM" || allcelltypesfine[i] == "Monocyte-derived macrophages") {
        data_lungs_all$immune_status[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Immune"
    }
    if (allcelltypesfine[i] == "NK cells" || allcelltypesfine[i] == "CD8+ T cells" || 
        allcelltypesfine[i] == "CD4+ T cells" || allcelltypesfine[i] == "Tregs" || 
        allcelltypesfine[i] == "Cycling NK/T cells") {
        data_lungs_all$immune_status[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Immune"
    }
    if (allcelltypesfine[i] == "Mast cells") {
        data_lungs_all$immune_status[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Immune"
    }
    
    # Group some cell_type_fine classifications into larger groups in
    # cell_type_intermediate
    if (str_detect(allcelltypesfine[i], "B cells")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "B cells"
    }
    if (str_detect(allcelltypesfine[i], "Other myeloid") || allcelltypesfine[i] == 
        "Monocyte-derived macrophages" || allcelltypesfine[i] == "Transitioning MDM") {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Macrophages"
    }
    if (str_detect(allcelltypesfine[i], "Airway") && !(str_detect(allcelltypesfine[i], 
        "Airway smooth muscle"))) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Airway epithelial cells"
    }
    if (str_detect(allcelltypesfine[i], "Inflammed epithelial") || str_detect(allcelltypesfine[i], 
        "AT intermediate")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Other epithelial cells"
    }
    if (str_detect(allcelltypesfine[i], "ndothelial")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Endothelial cells"
    }
    if (str_detect(allcelltypesfine[i], "FB")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Fibroblasts"
    }
    if (str_detect(allcelltypesfine[i], "smooth muscle")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Smooth muscle"
    }
    if (str_detect(allcelltypesfine[i], "Pericytes")) {
        data_lungs_all$cell_type_intermediate[data_lungs_all$cell_type_fine == allcelltypesfine[i]] = "Fibroblasts"
    }
}

# To validate cell type groupings, create matrices containing the total number of
# cells in each cell_type_intermediate vs fine pairing, main vs. fine, and main
# vs. intermediate
allcelltypesintermediate = unique(data_lungs_all$cell_type_intermediate)
checkmatrix = matrix(0, length(allcelltypesfine), length(allcelltypesintermediate), 
    dimnames = list(allcelltypesfine, allcelltypesintermediate))
for (i in 1:length(allcelltypesfine)) {
    for (j in 1:length(allcelltypesintermediate)) {
        acelltypemain = unique(data_lungs_all$cell_type_main[data_lungs_all$cell_type_intermediate == 
            allcelltypesintermediate[j]])
        if (length(acelltypemain) == 1) {
            acelltypemain = acelltypemain[1]
        } else {
            # print(allcelltypesintermediate[j]) nonsense = nonsense+1
        }
        # print(paste0(allcelltypesintermediate[j],' ',acelltypemain))
        checkmatrix[i, j] = sum(data_lungs_all$cell_type_fine == allcelltypesfine[i] & 
            data_lungs_all$cell_type_intermediate == allcelltypesintermediate[j])
    }
}

write.table(checkmatrix, "check_cell_type_intermediate_vs_fine.csv", sep = ",", quote = F, 
    row.names = T, col.names = T)

allcelltypesmain = unique(data_lungs_all$cell_type_main)
checkmatrix2 = matrix(0, length(allcelltypesfine), length(allcelltypesmain), dimnames = list(allcelltypesfine, 
    allcelltypesmain))
for (i in 1:length(allcelltypesfine)) {
    for (j in 1:length(allcelltypesmain)) {
        checkmatrix2[i, j] = sum(data_lungs_all$cell_type_fine == allcelltypesfine[i] & 
            data_lungs_all$cell_type_main == allcelltypesmain[j])
    }
}

write.table(checkmatrix2, "check_cell_type_main_vs_fine.csv", sep = ",", quote = F, 
    row.names = T, col.names = T)

checkmatrix3 = matrix(0, length(allcelltypesintermediate), length(allcelltypesmain), 
    dimnames = list(allcelltypesintermediate, allcelltypesmain))
for (i in 1:length(allcelltypesintermediate)) {
    for (j in 1:length(allcelltypesmain)) {
        checkmatrix3[i, j] = sum(data_lungs_all$cell_type_intermediate == allcelltypesintermediate[i] & 
            data_lungs_all$cell_type_main == allcelltypesmain[j])
    }
}

write.table(checkmatrix3, "check_cell_type_main_vs_intermediate.csv", sep = ",", 
    quote = F, row.names = T, col.names = T)
