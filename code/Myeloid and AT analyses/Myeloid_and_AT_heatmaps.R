#!/usr/bin/env Rscript

### title: Generate heatmaps of differential gene expression for AT cells,
### fibroblasts, or myeloid cells, comparing either different groups of cells, or
### comparing COVID-19 samples vs. Control samples within the same cell type
### author: Yiping Wang date: 02/08/2021

# Creates Figures 3E, Extended Data Figure 4H, Extended Data Figure 5A,
# Extended Data Fibroblasts I

# list of options controls which kind of heatmaps to plot choose whether to plot
# heatmaps for AT cells, fibroblasts, or myeloid (default)
use_AT = F
use_fibroblasts = F
# if plotting heatmaps of AT cells, choose whether or not to use DATP annotation
# as a separate class in heatmap
use_predatp = F
# if plotting heatmaps of AT cells, choose whether canonical AT cell markers
# should always be included in plot
includefixedmarkers = F
# if plotting two cell groups versus each other, choose whether to plot only
# covid or only control cells
covid_only_in_versus = F
control_only_in_versus = F
# choose whether to select only every 10th cell for plotting in pdf, to reduce
# file size
downsample = F
# choose whether to compare two different cell groups against each other
isversus = F
# choose whether to compare multiple different cell groups against each other
ismultiplecomp = F

# select data to plot based on options, and add a celltypefield that determines
# which subgroups of cells to plot
if (use_AT) {
    data_cov_ctr = data_AT_int_cov_ctr
    
    data_cov_ctr$celltype = ""
    if (!use_predatp) {
        data_cov_ctr$celltype[data_cov_ctr$cell_type_datp == "AT1"] = "AT1"
        data_cov_ctr$celltype[data_cov_ctr$cell_type_datp == "AT2"] = "AT2"
        data_cov_ctr$celltype[data_cov_ctr$cell_type_datp == "AT intermediate"] = "AT intermediate"
        data_cov_ctr$celltype[data_cov_ctr$cell_type_datp == "DATP"] = "DATP"
    } else {
        data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "AT1"] = "AT1"
        data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "AT2"] = "AT2"
        data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "AT intermediate"] = "AT intermediate"
    }
} else if (use_fibroblasts) {
    data_cov_ctr = data_lungs_all
    
    data_cov_ctr$celltype = ""
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Alveolar FB"] = "All other FB"
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Mesothelial FB"] = "All other FB"
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Other FB"] = "All other FB"
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Adventitial FB"] = "All other FB"
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Intermediate pathological FB"] = "Intermediate pathological FB"
    data_cov_ctr$celltype[data_cov_ctr$cell_type_fine == "Pathological FB"] = "Pathological FB"
    data_cov_ctr = subset(data_cov_ctr, celltype != "")
} else {
    data_cov_ctr = data_myeloid_cov_ctr
    
    data_cov_ctr$celltype = ""
    data_cov_ctr$celltype[data_cov_ctr$AM_specific == "Alveolar macrophages"] = "Alveolar macrophages"
    data_cov_ctr$celltype[data_cov_ctr$AM_specific == "Monocyte-derived macrophages"] = "Monocyte-derived macrophages"
    data_cov_ctr$celltype[data_cov_ctr$AM_specific == "Monocytes"] = "Monocytes"
    data_cov_ctr$celltype[data_cov_ctr$AM_specific == "Transitioning MDM"] = "Transitioning MDM"
    data_cov_ctr = subset(data_cov_ctr, celltype != "")
}

DefaultAssay(data_cov_ctr) = "RNA"
data_cov_ctr = ScaleData(data_cov_ctr, features = rownames(data_cov_ctr))

# list of comparisons to make, uncomment the line that corresponds to figure that
# should be made celltypestotest = c('AT1') celltypestotest = c('AT1','AT2')
# celltypestotest = c('AT1_vs_AT2') celltypestotest =
# list(list(c('AT1'),c('AT2'))) celltypestotest = c('Alveolar
# macrophages_vs_Monocyte-derived macrophages')
celltypestotest = list(c("Alveolar macrophages"))
# celltypestotest = list(c('Alveolar macrophages'),c('Monocyte-derived
# macrophages'),c('Transitioning MDM'),c('Monocytes')) celltypestotest =
# list(c('Alveolar macrophages','Monocyte-derived macrophages','Transitioning
# MDM','Monocytes')) celltypestotest = list(c('Monocyte-derived
# macrophages','Transitioning MDM')) celltypestotest =
# list(list(c('Monocytes','Monocyte-derived macrophages','Transitioning
# MDM'),c('Alveolar macrophages'))) celltypestotest = list(list(c('Alveolar
# macrophages'),c('Monocyte-derived macrophages'),c('Monocytes'),c('Transitioning
# MDM'))) celltypestotest = list(list(c('Alveolar
# macrophages'),c('Monocyte-derived macrophages'),c('Monocytes'),c('Other
# Fibroblasts'))) celltypestotest = list(list(c('Intermediate pathological
# FB'),c('Other FB'),c('Pathological FB'))) celltypestotest =
# list(list(c('Intermediate pathological FB','Pathological FB'),c('All other
# FB')))
AT2givenmarkers = c("SFTPB", "SFTPC", "SFTPD", "ETV5")
AT1givenmarkers = c("AGER", "CLIC5", "CAV1", "PDPN")
for (i in 1:length(celltypestotest)) {
    Idents(data_cov_ctr) = ""
    suffix = ""
    DefaultAssay(data_cov_ctr) <- "RNA"
    
    # use FindMarkers to identify markers between groups of cells, depending on
    # options given above, and write out results to csv files
    if (isversus) {
        celltype1 = celltypestotest[[i]][[1]]
        celltype2 = celltypestotest[[i]][[2]]
        if (covid_only_in_versus) {
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype1) & data_cov_ctr$group == "COVID-19"], value = "celltype1")
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype2) & data_cov_ctr$group == "COVID-19"], value = "celltype2")
            data_cov_ctr_small = data_cov_ctr[, ((data_cov_ctr$celltype %in% celltype1) | 
                (data_cov_ctr$celltype %in% celltype2)) & data_cov_ctr$group == "COVID-19"]
        } else if (control_only_in_versus) {
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype1) & data_cov_ctr$group == "Control"], value = "celltype1")
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype2) & data_cov_ctr$group == "Control"], value = "celltype2")
            data_cov_ctr_small = data_cov_ctr[, ((data_cov_ctr$celltype %in% celltype1) | 
                (data_cov_ctr$celltype %in% celltype2)) & data_cov_ctr$group == "Control"]
        } else {
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype1)], value = "celltype1")
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltype2)], value = "celltype2")
            data_cov_ctr_small = data_cov_ctr[, (data_cov_ctr$celltype %in% celltype1) | 
                (data_cov_ctr$celltype %in% celltype2)]
        }
        
        topmarkers1 = FindMarkers(data_cov_ctr, ident.1 = "celltype1", ident.2 = "celltype2", 
            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
        topmarkers2 = FindMarkers(data_cov_ctr, ident.1 = "celltype2", ident.2 = "celltype1", 
            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
        topmarkers1$gene = rownames(topmarkers1)
        topmarkers2$gene = rownames(topmarkers2)
        topmarkers1$celltype = str_replace_all(paste(celltype1, sep = "_", collapse = "_"), 
            " ", "_")
        topmarkers2$celltype = str_replace_all(paste(celltype2, sep = "_", collapse = "_"), 
            " ", "_")
        
        csvsuffix = ""
        if (covid_only_in_versus) {
            csvsuffix = "_covid_only"
        }
        if (control_only_in_versus) {
            csvsuffix = "_control_only"
        }
        
        if (use_predatp) {
            cat(paste0("Up in ", str_replace_all(paste(celltype1, sep = "_", collapse = "_"), 
                " ", "_"), "\n"), file = paste0(str_replace_all(paste(celltype1, 
                sep = "_", collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                append = FALSE)
            write.table(topmarkers1, paste0(str_replace_all(paste(celltype1, sep = "_", 
                collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                row.names = F, col.names = T, sep = ",", quote = F, append = T)
            cat(paste0("Up in ", str_replace_all(paste(celltype2, sep = "_", collapse = "_"), 
                " ", "_"), "\n"), file = paste0(str_replace_all(paste(celltype1, 
                sep = "_", collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                append = TRUE)
            write.table(topmarkers2, paste0(str_replace_all(paste(celltype1, sep = "_", 
                collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                row.names = F, col.names = T, sep = ",", quote = F, append = T)
        } else {
            cat(paste0("Up in ", str_replace_all(paste(celltype1, sep = "_", collapse = "_"), 
                " ", "_"), "\n"), file = paste0(str_replace_all(paste(celltype1, 
                sep = "_", collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                append = FALSE)
            write.table(topmarkers1, paste0(str_replace_all(paste(celltype1, sep = "_", 
                collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                row.names = F, col.names = T, sep = ",", quote = F, append = T)
            cat(paste0("Up in ", str_replace_all(paste(celltype2, sep = "_", collapse = "_"), 
                " ", "_"), "\n"), file = paste0(str_replace_all(paste(celltype1, 
                sep = "_", collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                append = TRUE)
            write.table(topmarkers2, paste0(str_replace_all(paste(celltype1, sep = "_", 
                collapse = "_"), " ", "_"), "_vs_", str_replace_all(paste(celltype2, 
                sep = "_", collapse = "_"), " ", "_"), csvsuffix, "_predatp_markers.csv"), 
                row.names = F, col.names = T, sep = ",", quote = F, append = T)
        }
        
        if (includefixedmarkers) {
            topmarkers1 = c(rownames(topmarkers2)[1:25], rownames(topmarkers1)[1:25])
        } else {
            topmarkers1 = c(rownames(topmarkers1)[1:25], rownames(topmarkers2)[1:25])
        }
    } else if (ismultiplecomp) {
        for (j in 1:length(celltypestotest[[i]])) {
            celltypesjth = celltypestotest[[i]][[j]]
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
                celltypesjth)], value = "celltypesjth")
            data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[!(data_cov_ctr$celltype %in% 
                celltypesjth)], value = "celltypesother")
            data_cov_ctr_small = data_cov_ctr
            
            topmarkerstemp = FindMarkers(data_cov_ctr, ident.1 = "celltypesjth", 
                ident.2 = "celltypesother", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                assay = "RNA")
            topmarkerstemp$gene = rownames(topmarkerstemp)
            topmarkerstemp$celltype = str_replace_all(paste(celltypesjth, sep = "_", 
                collapse = "_"), " ", "_")
            
            if (j == 1) {
                topmarkers1 = rownames(topmarkerstemp)[1:25]
            } else {
                topmarkers1 = c(topmarkers1, rownames(topmarkerstemp)[1:25])
            }
        }
    } else {
        data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
            celltypestotest[[i]]) & data_cov_ctr$group == "Control"], value = "Control")
        data_cov_ctr = SetIdent(data_cov_ctr, cells = colnames(data_cov_ctr)[(data_cov_ctr$celltype %in% 
            celltypestotest[[i]]) & data_cov_ctr$group == "COVID-19"], value = "COVID-19")
        data_cov_ctr_small = data_cov_ctr[, data_cov_ctr$celltype %in% celltypestotest[[i]]]
        topmarkers1 = FindMarkers(data_cov_ctr, ident.1 = "Control", ident.2 = "COVID-19", 
            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
        topmarkers2 = FindMarkers(data_cov_ctr, ident.1 = "COVID-19", ident.2 = "Control", 
            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
        topmarkers1$gene = rownames(topmarkers1)
        topmarkers2$gene = rownames(topmarkers2)
        topmarkers1$condition = "Control"
        topmarkers2$condition = "COVID-19"
        
        csvsuffix = "_markers.csv"
        if (use_AT && use_predatp) {
            csvsuffix = "_predatp_markers.csv"
        }
        
        cat(paste0("Up in Control\n"), file = paste0(str_replace_all(paste(celltypestotest[[i]], 
            sep = "_", collapse = "_"), " ", "_"), csvsuffix), append = FALSE)
        write.table(topmarkers1, paste0(str_replace_all(paste(celltypestotest[[i]], 
            sep = "_", collapse = "_"), " ", "_"), csvsuffix), row.names = F, col.names = T, 
            sep = ",", quote = F, append = T)
        cat(paste0("Up in COVID-19\n"), file = paste0(str_replace_all(paste(celltypestotest[[i]], 
            sep = "_", collapse = "_"), " ", "_"), csvsuffix), append = TRUE)
        write.table(topmarkers2, paste0(str_replace_all(paste(celltypestotest[[i]], 
            sep = "_", collapse = "_"), " ", "_"), csvsuffix), row.names = F, col.names = T, 
            sep = ",", quote = F, append = T)
        topmarkers1 = c(rownames(topmarkers1)[1:25], rownames(topmarkers2)[1:25])
    }
    
    # include canonical AT markers in marker list if desired
    if (includefixedmarkers) {
        if ("CAV1" %in% topmarkers1) {
            AT1givenmarkers = c("AGER", "CLIC5", "PDPN")
        }
        topmarkers1 = topmarkers1[!(topmarkers1 %in% AT2givenmarkers)]
        topmarkers1 = topmarkers1[!(topmarkers1 %in% AT1givenmarkers)]
        data_cov_ctr_AT2givenidxs = match(AT2givenmarkers, rownames(data_cov_ctr_small))
        data_cov_ctr_AT1givenidxs = match(AT1givenmarkers, rownames(data_cov_ctr_small))
    }
    
    # find locations of marker genes within seurat data object
    data_cov_ctr_small1idxs = match(topmarkers1, rownames(data_cov_ctr_small))
    
    # create annotation status vector that records how markers were generated and for
    # which groups
    if (includefixedmarkers) {
        data_cov_ctr_small = data_cov_ctr_small[c(data_cov_ctr_AT2givenidxs, data_cov_ctr_AT1givenidxs, 
            data_cov_ctr_small1idxs), ]
        
        annotation_status = c(rep("AT2 original marker", length(data_cov_ctr_AT2givenidxs)), 
            rep("AT1 original marker", length(data_cov_ctr_AT1givenidxs)), rep("AT1/2 inferred markers", 
                length(data_cov_ctr_small1idxs)))
    } else {
        data_cov_ctr_small = data_cov_ctr_small[c(data_cov_ctr_small1idxs), ]
        if (isversus) {
            if (use_AT) {
                annotation_status = c(rep("AT1/2 inferred markers", length(data_cov_ctr_small1idxs)))
            } else {
                annotation_status = c(rep(paste0(paste(celltype1, sep = ", ", collapse = ", "), 
                  " versus ", paste(celltype2, sep = ", ", collapse = ", "), " inferred markers"), 
                  length(data_cov_ctr_small1idxs)))
            }
        } else if (ismultiplecomp) {
            annotation_status = c()
            for (j in 1:length(celltypestotest[[i]])) {
                celltypesjth = celltypestotest[[i]][[j]]
                annotation_status = c(annotation_status, rep(str_replace_all(paste(celltypesjth, 
                  sep = "_", collapse = "_"), " ", "_"), 25))
            }
        } else {
            annotation_status = c(rep("Control versus COVID-19 markers", length(data_cov_ctr_small1idxs)))
        }
        suffix = paste0(suffix, "nofixed_markers_")
        if (covid_only_in_versus && isversus) {
            suffix = paste0(suffix, "covid_only_")
        }
        if (control_only_in_versus && isversus) {
            suffix = paste0(suffix, "control_only_")
        }
    }
    
    # rearrange seurat data object so that gene names in $RNA@scale.data match up
    # with rest of object
    permuteval = match(rownames(data_cov_ctr_small), rownames(data_cov_ctr_small$RNA@scale.data))
    data_cov_ctr_small$RNA@scale.data = data_cov_ctr_small$RNA@scale.data[permuteval, 
        ]
    
    # reorder cell types so that it matches order of cell types to be compared, then
    # covid status within each cell type
    if (isversus) {
        reordercelltypes = c(celltype1, celltype2)
    } else {
        reordercelltypes = sort(unique(data_cov_ctr_small$celltype))
    }
    reorderarr = c()
    for (j in 1:length(reordercelltypes)) {
        reorderarrtemp = which(data_cov_ctr_small$celltype == reordercelltypes[j])
        reorderarrtemp = reorderarrtemp[order(data_cov_ctr_small$group[reorderarrtemp])]
        reorderarr = c(reorderarr, reorderarrtemp)
    }
    
    data_cov_ctr_small = data_cov_ctr_small[, reorderarr]
    if (ismultiplecomp) {
        datamat = data_cov_ctr$RNA@scale.data[which(rownames(data_cov_ctr) %in% topmarkers1), 
            ]
        datamat = datamat[match(topmarkers1, rownames(datamat)), ]
    } else {
        datamat = data_cov_ctr_small$RNA@scale.data
    }
    
    # extract signature strength for each cell, to be plotted on top of heatmap if
    # analyzing AT cells
    celltype = data_cov_ctr_small$celltype
    groupvec = data_cov_ctr_small$group
    if (use_AT) {
        DATP_short = data_cov_ctr_small$DATP_short1
        our_DATP_sig = data_cov_ctr_small$our_DATP_sig1
        AT1 = data_cov_ctr_small$AT11
        AT2 = data_cov_ctr_small$AT21
        primed_AT2 = data_cov_ctr_small$primed_AT21
        cycling_AT2 = data_cov_ctr_small$cycling_AT21
        UP_in_DATPs = data_cov_ctr_small$UP_in_DATPs1
    }
    
    # reorder data and column labels by cell type, then COVID-19 status, as above
    datamat = datamat[, reorderarr]
    celltype = celltype[reorderarr]
    groupvec = groupvec[reorderarr]
    if (use_AT) {
        DATP_short = DATP_short[reorderarr]
        our_DATP_sig = our_DATP_sig[reorderarr]
        AT1 = AT1[reorderarr]
        AT2 = AT2[reorderarr]
        primed_AT2 = primed_AT2[reorderarr]
        cycling_AT2 = cycling_AT2[reorderarr]
        UP_in_DATPs = UP_in_DATPs[reorderarr]
    }
    
    if (use_AT) {
        annotation_col = data.frame(group = factor(groupvec), celltype = factor(celltype), 
            DATP_sig = our_DATP_sig, primed_AT2 = primed_AT2, AT2 = AT2, AT1 = AT1)
    } else {
        annotation_col = data.frame(group = factor(groupvec), celltype = factor(celltype))
    }
    
    rownames(annotation_col) = colnames(datamat)
    
    # if multiple cell groups, some marker genes may be repeated within data matrix
    # for repeat genes, add number at end to distinguish repeats
    if (ismultiplecomp) {
        labels_row = rownames(datamat)
        annotation_row = data.frame(annotation_status = annotation_status)
        uniquegenenames = unique(rownames(datamat))
        for (z in 1:length(uniquegenenames)) {
            if (sum(rownames(datamat) == uniquegenenames[z]) > 1) {
                matchidxs = which(rownames(datamat) == uniquegenenames[z])
                for (j in 1:length(matchidxs)) {
                  rownames(datamat)[matchidxs[j]] = paste0(rownames(datamat)[matchidxs[j]], 
                    "_", j)
                }
            }
        }
        rownames(annotation_row) = rownames(datamat)
    } else {
        annotation_row = data.frame(annotation_status = annotation_status)
        
        rownames(annotation_row) = rownames(datamat)
    }
    
    # assign colors to each cell type that may be present in data
    logdatamat = log(datamat - 1.01 * min(min(datamat[datamat != 0])))
    consistentcolors = c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", "#005AC8", 
        "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", "#FA7850", 
        "#14D2DC", "#FA78FA")
    consistentcolors2 = c("#8214A0", "#0072B2", "#0AB45A", "#FA5078")
    celltypecolors = c()
    if (sum(celltype == "AT intermediate") != 0) {
        celltypecolors = append(celltypecolors, c(`AT intermediate` = consistentcolors2[1]))
    }
    if (sum(celltype == "AT1") != 0) {
        celltypecolors = append(celltypecolors, c(AT1 = consistentcolors2[2]))
    }
    if (sum(celltype == "AT2") != 0) {
        celltypecolors = append(celltypecolors, c(AT2 = consistentcolors2[3]))
    }
    if (sum(celltype == "DATP") != 0) {
        celltypecolors = append(celltypecolors, c(DATP = consistentcolors2[4]))
    }
    consistentcolors3 = c("#E73F74", "#F2B701", "#11A579", "#3969AC")
    if (sum(celltype == "Alveolar macrophages") != 0) {
        celltypecolors = append(celltypecolors, c(`Alveolar macrophages` = consistentcolors3[1]))
    }
    if (sum(celltype == "Monocyte-derived macrophages") != 0) {
        celltypecolors = append(celltypecolors, c(`Monocyte-derived macrophages` = consistentcolors3[2]))
    }
    if (sum(celltype == "Transitioning MDM") != 0) {
        celltypecolors = append(celltypecolors, c(`Transitioning MDM` = consistentcolors3[3]))
    }
    if (sum(celltype == "Monocytes") != 0) {
        celltypecolors = append(celltypecolors, c(Monocytes = consistentcolors3[4]))
    }
    if (sum(celltype == "Pathological FB") != 0) {
        celltypecolors = append(celltypecolors, c(`Pathological FB` = consistentcolors3[1]))
    }
    if (sum(celltype == "Intermediate pathological FB") != 0) {
        celltypecolors = append(celltypecolors, c(`Intermediate pathological FB` = consistentcolors3[2]))
    }
    if (sum(celltype == "All other FB") != 0) {
        celltypecolors = append(celltypecolors, c(`All other FB` = consistentcolors3[3]))
    }
    ann_colors = list(group = c(Control = consistentcolors[1], `COVID-19` = consistentcolors[2]), 
        celltype = celltypecolors)
    
    # use PurpleAndYellow palette to emulate Seurat heatmap color scheme
    colorpalette = PurpleAndYellow()
    realdatamax = max(max(datamat))
    realdatamin = min(min(datamat))
    colorpalette2 = colorpalette
    
    # clip heatmap values to plus or minus 2.5, to ensure that color contrast is high
    datamat[datamat <= -2.5] = -2.5
    datamat[datamat >= 2.5] = 2.5
    
    if (isversus) {
        printname = str_replace_all(paste0(paste(celltype1, sep = "_", collapse = "_"), 
            "_vs_", paste(celltype2, sep = "_", collapse = "_")), " ", "_")
    } else if (ismultiplecomp) {
        printname = str_replace_all(paste(celltypestotest[[i]], sep = "_", collapse = "_"), 
            " ", "_")
    } else {
        printname = str_replace_all(paste(celltypestotest[[i]], sep = "_", collapse = "_"), 
            " ", "_")
    }
    
    print(printname)
    # print heatmap to pdf file
    if (downsample) {
        suffix = paste0(suffix, "downsample_")
        annotation_col = annotation_col[seq(1, dim(datamat)[2], 10), ]
        datamat = datamat[, seq(1, dim(datamat)[2], 10)]
    }
    if (use_AT && use_predatp) {
        pdf(paste0(str_replace_all(printname, "-", "_"), "_predatp_cov_ctr_", suffix, 
            "heatmap.pdf"), width = 14, height = 14)
    } else if (ismultiplecomp) {
        pdf(paste0(str_replace_all(printname, "-", "_"), "_", suffix, "heatmap.pdf"), 
            width = 14, height = 14)
    } else {
        pdf(paste0(str_replace_all(printname, "-", "_"), "_cov_ctr_", suffix, "heatmap.pdf"), 
            width = 14, height = 14)
    }
    if (ismultiplecomp) {
        print(pheatmap(datamat, color = PurpleAndYellow(), annotation_col = annotation_col, 
            annotation_row = annotation_row, annotation_colors = ann_colors, labels_row = labels_row, 
            show_colnames = F, cluster_rows = F, cluster_cols = F, useRaster = T))
    } else {
        print(pheatmap(datamat, color = colorpalette2, annotation_col = annotation_col, 
            annotation_row = annotation_row, annotation_colors = ann_colors, show_colnames = F, 
            cluster_rows = F, cluster_cols = F, useRaster = T))
    }
    dev.off()
    
    # print violin plot of AXL expression in Alveolar macrophages, if options are
    # correct
    if (!(isversus) && paste(celltypestotest[[i]], sep = "_", collapse = "_") == 
        "Alveolar macrophages") {
        pdf("Extended Data Myeloid 2 b.pdf")
        subsetdata = subset(data_cov_ctr, AM_specific == "Alveolar macrophages")
        wilcoxpval = wilcox.test(subsetdata$RNA@scale.data[rownames(subsetdata) == 
            "AXL", subsetdata$group == "Control"], subsetdata$RNA@scale.data[rownames(subsetdata) == 
            "AXL", subsetdata$group == "COVID-19"])$p.val * 34546
        print(VlnPlot(subsetdata, features = "AXL", group.by = "group", pt.size = 0, 
            cols = consistentcolors[1:2], log = F) + ggtitle(paste0("AXL Control vs. COVID-19 Wilcoxon p-val: ", 
            toString(wilcoxpval))))
        dev.off()
    }
}
