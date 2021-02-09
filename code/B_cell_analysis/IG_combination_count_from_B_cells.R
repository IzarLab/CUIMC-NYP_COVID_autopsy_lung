#!/usr/bin/env Rscript

### title: Counting number of IG combinations in B cells authors: Jana Biermann,
### Yiping Wang date: 10/05/20

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)

# Run analysis for only covid samples, only control, and all samples
groups = c("cov", "ctr", "")
longgroupnames = c("COVID-19", "Control", "")

# construct result dataframe that stores light and heavy chain information for
# every cell, and which patient the combination is found in
res_patient = data.frame(lightchain = character(), heavychain = character(), orig.ident = character())
for (z in 1:length(groups)) {
    patient <- readRDS(paste0("/data/data_lungs_all_bcells.rds"))
    
    if (groups[z] != "") {
        patient = subset(patient, group == groups[z])
    }
    
    # Rerun Seurat workflow
    patient <- FindVariableFeatures(patient, selection.method = "vst", nfeatures = 2000)
    patient <- ScaleData(object = patient)
    patient <- RunPCA(object = patient)
    patient <- FindNeighbors(patient, dims = 1:20)
    patient <- FindClusters(patient, resolution = 1.8)
    patient <- RunUMAP(object = patient, dims = 1:20)
    
    ig <- grep("^IG", rownames(patient@assays$RNA@counts), value = T)
    
    # Create gene expression dataframe with only IG genes, delete empty rows and
    # columns
    tab <- as.data.frame(patient@assays$RNA@data[ig, ])
    nonzerocelltypes = patient$celltype_bped_main
    nonzerocelltypes = nonzerocelltypes[colSums(tab) > 0]
    tab <- tab[, colSums(tab) > 0]
    tab <- tab[rowSums(tab) > 0, ]
    
    # List heavy and light variable regions
    ighv <- grep("^IGHV", rownames(tab), value = T)
    iglv <- grep("^(IGLV|IGKV)", rownames(tab), value = T)
    
    # Set up empty dataframe for results
    res <- matrix(data = 0, nrow = length(iglv), ncol = length(ighv))
    colnames(res) <- ighv
    rownames(res) <- iglv
    res <- as.data.frame(res)
    
    # List constant regions, manually remove the ones that don't follow the right
    # pattern, and set up result dataframe for constant chains
    ig_con <- grep("^(IGHG|IGHA|IGHM|IGHE)", rownames(tab), value = T)
    ig_con <- ig_con[ig_con != "IGHGP"]
    res_con <- as.data.frame(matrix(data = NA, nrow = 1000, ncol = 5))
    colnames(res_con) <- c("barcode", "orig.ident", "heavy", "light", "constant")
    
    # Loop to identify cells with variable and constant chains expressed
    i <- 1
    nonzerocelltypes1 = c()
    for (c in 1:ncol(tab)) {
        # Remove values from previous loop
        rm("igl")
        rm("igh")
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
        
        # Only add cell to variable regions result dataframe if it has both heavy and
        # light variable chains expressed
        if (cell[igh] > 0 & cell[igl] > 0) {
            res[igl, igh] <- res[igl, igh] + 1
            nonzerocelltypes1 = append(nonzerocelltypes1, nonzerocelltypes[c])
            if (groups[z] != "") {
                anorig.ident = patient$orig.ident[match(colnames(tab)[c], colnames(tab))]
                res_patient_temp = data.frame(lightchain = igl, heavychain = igh, 
                  orig.ident = anorig.ident)
                res_patient = rbind(res_patient, res_patient_temp)
            }
        }
        
        # Only add cell to constant regions result dataframe if it has heavy and light
        # variable chains and a constant chain expressed
        if (cell[igh] > 0 & cell[igl] > 0 & cell[con] > 0) {
            res_con$barcode[i] <- colnames(tab)[c]
            res_con$heavy[i] <- igh
            res_con$light[i] <- igl
            res_con$constant[i] <- con
            res_con$orig.ident[i] <- patient$orig.ident[colnames(patient) == colnames(tab)[c]]
            i <- i + 1
        }
    }
    
    # Remove columns and rows that are empty
    res <- res[, colSums(res) > 0]
    res <- res[rowSums(res) > 0, ]
    
    # Add IG subtype info
    res_con$subtype <- substr(res_con$constant, 1, 4)
    if (groups[z] != "") {
        write.csv(res_con, paste0("table_IG_new_lungs_all_b_cell_int_all_chains_", 
            groups[z], ".csv"), row.names = F, quote = F)
    } else {
        write.csv(res_con, paste0("table_IG_new_lungs_all_b_cell_int_all_chains.csv"), 
            row.names = F, quote = F)
    }
    
    my_palette <- colorRampPalette(c("gray", "green"))(n = 100)
    # Save Figure 2h change height dimensions so light labels on y axis don't have
    # some dropout if figure too short
    if (groups[z] != "") {
        if (groups[z] == "cov") {
            pdf(paste0("IG_combination_frequency_", groups[z], ".pdf"), height = 12)
        } else {
            pdf(paste0("IG_combination_frequency_", groups[z], ".pdf"))
        }
    } else {
        pdf(paste0("IG_combination_frequency.pdf"), height = 12)
    }
    
    # Extended Data Figure B-cells d
    out = heatmap.2(as.matrix(res), trace = "none", scale = "none", margins = c(6, 
        6), cexRow = 0.7, cexCol = 0.7, hclustfun = function(d) {
        hclust(d, method = "average")
    }, xlab = "Heavy chains", ylab = "Light chains", col = rev(heat.colors(100)), 
        keysize = 1.5, lhei = c(1.7, 5), density.info = "none")
    dev.off()
    
    # write out light and heavy chains as plotted in combination frequency heatmap,
    # for use in IG_combination_plots.R
    if (groups[z] == "") {
        lightorder = rev(rownames(res)[out$rowInd])
        write.table(lightorder, "lightorder.txt", col.names = F, row.names = F, quote = F, 
            sep = "\t")
        heavyorder = colnames(res)[out$colInd]
        write.table(heavyorder, "heavyorder.txt", col.names = F, row.names = F, quote = F, 
            sep = "\t")
    }
    
    # write out overall frequency of light and heavy chains, in COVID-19 and Control
    # groups separately
    lightoveralldf = data.frame(lightchain = rownames(res), frequency = rowSums(res), 
        group = longgroupnames[z])
    heavyoveralldf = data.frame(heavychain = colnames(res), frequency = colSums(res), 
        group = longgroupnames[z])
    if (z == 1) {
        write.table(lightoveralldf, "lightoverallfreq.txt", row.names = F, col.names = T, 
            quote = F, sep = "\t")
        write.table(heavyoveralldf, "heavyoverallfreq.txt", row.names = F, col.names = T, 
            quote = F, sep = "\t")
    }
    if (z == 2) {
        write.table(lightoveralldf, "lightoverallfreq.txt", row.names = F, col.names = T, 
            quote = F, sep = "\t", append = T)
        write.table(heavyoveralldf, "heavyoverallfreq.txt", row.names = F, col.names = T, 
            quote = F, sep = "\t", append = T)
    }
}

consistentcolors = c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", "#005AC8", 
    "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", "#FA7850", 
    "#14D2DC", "#FA78FA")

# read in light chain overall frequency df make a new dataframe with light chain
# name linked to frequency in both COVID-19 and Control groups, summed
theme_set(theme_bw())
lightfreqdf = read.table("lightoverallfreq.txt", header = T, sep = "\t", quote = NULL)
uniquelightchain = unique(lightfreqdf$lightchain)
lightsortdf = data.frame(lightchain = uniquelightchain, frequency = 0)
for (i in 1:length(uniquelightchain)) {
    lightsortdf$frequency[lightsortdf$lightchain == uniquelightchain[i]] = sum(lightfreqdf$frequency[lightfreqdf$lightchain == 
        uniquelightchain[i]])
}

# same as above, for heavy chains
heavyfreqdf = read.table("heavyoverallfreq.txt", header = T, sep = "\t", quote = NULL)
uniqueheavychain = unique(heavyfreqdf$heavychain)
heavysortdf = data.frame(heavychain = uniqueheavychain, frequency = 0)
for (i in 1:length(uniqueheavychain)) {
    heavysortdf$frequency[heavysortdf$heavychain == uniqueheavychain[i]] = sum(heavyfreqdf$frequency[heavyfreqdf$heavychain == 
        uniqueheavychain[i]])
}

# count number of light chain occurrences in each patient
res_patient_tally_light = res_patient %>% group_by(lightchain, orig.ident) %>% tally()

# print barplot of light chain usage in each patient, sorted by overall light
# chain usage across all patients, using information in lightsortdf
ggplot(res_patient_tally_light, aes(fill = orig.ident, y = n, x = factor(lightchain, 
    level = lightsortdf$lightchain[rev(order(lightsortdf$frequency))]))) + geom_bar(position = "stack", 
    stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Light Chain") + ylab("Usage")
ggsave("Extended Data Figure B-cells 2i.pdf", width = 14)

# print light chain usage in COVID-19 and Control groups, sorted by overall usage
# of light chains across both groups
ggplot(lightfreqdf, aes(fill = group, y = frequency, x = factor(lightchain, level = lightsortdf$lightchain[rev(order(lightsortdf$frequency))]))) + 
    geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_fill_manual(values = consistentcolors[1:2], breaks = c("COVID-19", 
    "Control")) + xlab("Light Chain")
ggsave("Figure_2h.pdf", width = 14)

# similar analyses as above, for heavy chains
res_patient_tally_heavy = res_patient %>% group_by(heavychain, orig.ident) %>% tally()

ggplot(heavyfreqdf, aes(fill = group, y = frequency, x = factor(heavychain, level = heavysortdf$heavychain[rev(order(heavysortdf$frequency))]))) + 
    geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_fill_manual(values = consistentcolors[1:2], breaks = c("COVID-19", 
    "Control")) + xlab("Heavy Chain")
ggsave("Figure_2i.pdf", width = 14)

ggplot(res_patient_tally_heavy, aes(fill = orig.ident, y = n, x = factor(heavychain, 
    level = heavysortdf$heavychain[rev(order(heavysortdf$frequency))]))) + geom_bar(position = "stack", 
    stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Heavy Chain") + ylab("Usage")
ggsave("Extended Data Figure B-cells h.pdf", width = 14)
