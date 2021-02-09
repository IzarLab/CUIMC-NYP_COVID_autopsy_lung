#!/usr/bin/env Rscript

# script for 
# 1. identifying conserved ligand-receptor interaction across different samples which were obtained by running CellPhoneDB for each individual sample and 
# 2. computing log2FC and p-val of the ligand-receptor interaction between COVID and control
# 3. identifying the top interacting ligand-receptors between COVID and control
# 4. plotting log2FC and p-value of selected interactions between COVID and control


library(pheatmap)
library(ggplot2)


#*****************************************************************************************
#identify common ligand-receptor interactions across COIVD and Ctrl samples
#*****************************************************************************************

#*#identify common ligand-receptor interactions across different Ctrl samples
path1 <- c("cellphoneDB/Ctrl")
path2 <- c("out_cellTypeMain")
dirs_available <- list.dirs(path = path1,full.names = F, recursive = F)
dirs_available


# getting the common interactions and cell-cells
# patient <- dirs_available[1]
intr_pairs <- NULL
cell_cells <- NULL
for (patient in dirs_available) {
  print(patient)
  pvalues_path <-  file.path(path1,patient,path2,"pvalues.txt")
  print(pvalues_path)
  
  #reading the cpdb data
  temp_all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  print(dim(temp_all_pval))
  temp_all_pval[1:10,1:15]
  
  if(!is.null(intr_pairs)){
    intr_pairs = intersect(intr_pairs, temp_all_pval$interacting_pair)
    cell_cells = intersect(cell_cells, colnames(temp_all_pval)[12:ncol(temp_all_pval)])
  } else{
    intr_pairs = temp_all_pval$interacting_pair
    cell_cells = colnames(temp_all_pval)[12:ncol(temp_all_pval)]
  }
}
intr_pairs_ctrl <- intr_pairs
cell_cells_ctrl <- cell_cells
length(intr_pairs_ctrl)
head(intr_pairs_ctrl)
length(cell_cells_ctrl)
head(cell_cells_ctrl)



#*#identify common ligand-receptor interactions across different COVID samples

path1 <- c("cellphoneDB/Covid")
path2 <- c("out_cellTypeMain")
dirs_available <- list.dirs(path = path1,full.names = F, recursive = F)
dirs_available


# getting the common interactions and cell-cells
# patient <- dirs_available[1]
intr_pairs <- NULL
cell_cells <- NULL
for (patient in dirs_available) {
  print(patient)
  pvalues_path <-  file.path(path1,patient,path2,"pvalues.txt")
  print(pvalues_path)
  
  #reading the cpdb data
  temp_all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  print(dim(temp_all_pval))
  temp_all_pval[1:10,1:15]
  
  if(!is.null(intr_pairs)){
    intr_pairs = intersect(intr_pairs, temp_all_pval$interacting_pair)
    cell_cells = intersect(cell_cells, colnames(temp_all_pval)[12:ncol(temp_all_pval)])
  } else{
    intr_pairs = temp_all_pval$interacting_pair
    cell_cells = colnames(temp_all_pval)[12:ncol(temp_all_pval)]
  }
}
intr_pairs_covid <- intr_pairs
cell_cells_covid <- cell_cells
length(intr_pairs_covid)
head(intr_pairs_covid)
length(cell_cells_covid)
head(cell_cells_covid)


# get the common ligand-receptors and cell-cells between covid and ctrl
intr_pairs <- intersect(intr_pairs_ctrl,intr_pairs_covid)
cell_cells <- intersect(cell_cells_ctrl,cell_cells_covid)

length(intr_pairs)
length(cell_cells)
head(intr_pairs)
head(cell_cells)


#*****************************************************************************************
#************************#getting Ctrl interactions
#*****************************************************************************************
path1 <- c("cellphoneDB/Ctrl")
path2 <- c("out_cellTypeMain")
dirs_available <- list.dirs(path = path1,full.names = F, recursive = F)
dirs_available

# reading the common pvals and means
# patient <- dirs_available[1]
means_patients_val <- list()
for (patient in dirs_available) {
  print(patient)
  means_path <-  file.path(path1,patient,path2,"means.txt")
  print(means_path)
  
  if(file.exists(means_path)){
    temp_all_means = read.table(means_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
    print(dim(temp_all_means))
    temp_all_means[1:10,1:15]
    temp_pos1 <- match(intr_pairs,temp_all_means$interacting_pair)
    temp_pos2 <- match(cell_cells,colnames(temp_all_means))
    means_patients_val[[patient]] = temp_all_means[temp_pos1,temp_pos2]
    rownames(means_patients_val[[patient]]) <- temp_all_means$interacting_pair[temp_pos1]
    head(means_patients_val[[patient]])
    dim(means_patients_val[[patient]])
  } else {print(paste0("File not found: ",means_path))}
}

length(means_patients_val)
names(means_patients_val)
head(means_patients_val[[1]])

# read the means from all Control samples
# boxplot(means_patients_val[[1]])
temp <- as.matrix(do.call(cbind,means_patients_val))
dim(temp)
# boxplot(temp)
means_patients_val_ctrl <- array(temp,dim = c(dim(means_patients_val[[1]]), length(means_patients_val)))
dim(means_patients_val_ctrl)
CtrlSamples <- dim(means_patients_val_ctrl)[3]
# means_patients_val_ctrl[1:5,1:5,]
# means_patients_val_ctrl[1,1,1:7]

# compute the median of all means
median_means_patient_val_ctrl <- apply(means_patients_val_ctrl, c(1,2), median, na.rm=T)
dim(median_means_patient_val_ctrl)
# median_means_patient_val_ctrl
colnames(median_means_patient_val_ctrl) <- cell_cells
rownames(median_means_patient_val_ctrl) <- intr_pairs
median_means_patient_val_ctrl[1:10,1:10]
# boxplot(median_means_patient_val_ctrl)

# # dim(median_means_patient_val_ctrl)
# dim(all_means_ctrl)
# # median_means_patient_val_ctrl[1:5,1:5]
# all_means_ctrl[1:5,1:5]
# temp_pos1 <- match(intr_pairs, all_means_ctrl$interacting_pair)
# temp_pos2 <- match(cell_cells, colnames(all_means_ctrl))
# # all(median_means_patient_val_ctrl==all_means_ctrl[temp_pos1, temp_pos2])
# median_means_patient_val_ctrl <- all_means_ctrl[temp_pos1, temp_pos2]
# colnames(median_means_patient_val_ctrl) <- cell_cells
# rownames(median_means_patient_val_ctrl) <- intr_pairs
# remove(temp,means_patients_val ) 
# median_means_patient_val_ctrl[1:5,1:5]
# dim(median_means_patient_val_ctrl)


#*****************************************************************************************
#************************#getting Covid interactions 
#*****************************************************************************************
path1 <- c("cellphoneDB/Covid")
path2 <- c("out_cellTypeMain")
dirs_available <- list.dirs(path = path1,full.names = F, recursive = F)
dirs_available


# reading the common pvals and means
# patient <- dirs_available[1]
means_patients_val <- list()
for (patient in dirs_available) {
  print(patient)
  means_path <-  file.path(path1,patient,path2,"means.txt")
  print(means_path)
  
  if(file.exists(means_path)){
    temp_all_means = read.table(means_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
    print(dim(temp_all_means))
    temp_all_means[1:10,1:15]
    temp_pos1 <- match(intr_pairs,temp_all_means$interacting_pair)
    temp_pos2 <- match(cell_cells,colnames(temp_all_means))
    means_patients_val[[patient]] = temp_all_means[temp_pos1,temp_pos2]
    rownames(means_patients_val[[patient]]) <- temp_all_means$interacting_pair[temp_pos1]
    head(means_patients_val[[patient]])
    dim(means_patients_val[[patient]])
  } else {print(paste0("File not found: ",means_path))}
}

length(means_patients_val)
names(means_patients_val)
head(means_patients_val[[1]])

# get the means of all Control samples
# boxplot(means_patients_val[[1]])
temp <- as.matrix(do.call(cbind,means_patients_val))
dim(temp)
# boxplot(temp)
means_patients_val_covid <- array(temp,dim = c(dim(means_patients_val[[1]]), length(means_patients_val)))
# temp <-  means_patients_val_covid
# means_patients_val_covid <- temp
# testing
# means_patients_val_covid <- means_patients_val_covid[31:35,31:35,1:7]
# means_patients_val_covid
dim(means_patients_val_covid)
covidSamples <- dim(means_patients_val_covid)[3]
# means_patients_val_covid[1:5,1:5,]
# means_patients_val_covid[1,1,1:covidSamples]

# get the median
median_means_patient_val_covid <- apply(means_patients_val_covid, c(1,2), median, na.rm=T)
dim(median_means_patient_val_covid)
# dim(median_means_patient_val_covid)
# dim(all_means_covid)
# median_means_patient_val_covid[1:5,1:5]
# all_means_covid[1:5,1:5]
# temp_pos1 <- match(intr_pairs, all_means_covid$interacting_pair)
# temp_pos2 <- match(cell_cells, colnames(all_means_covid))
# # all(median_means_patient_val_covid==all_means_covid[temp_pos1, temp_pos2])
# median_means_patient_val_covid <- all_means_covid[temp_pos1, temp_pos2]
colnames(median_means_patient_val_covid) <- cell_cells
rownames(median_means_patient_val_covid) <- intr_pairs
median_means_patient_val_covid[1:5,1:5]
remove(temp,means_patients_val ) 
dim(median_means_patient_val_covid)

all(colnames(median_means_patient_val_covid)==colnames(median_means_patient_val_ctrl))
all(rownames(median_means_patient_val_covid)==rownames(median_means_patient_val_ctrl))

dim(means_patients_val_ctrl)
CtrlSamples
dim(means_patients_val_covid)
covidSamples


#********************************************** compute the pvalues using Wilcox test
nrows <- length(intr_pairs)
ncols <- length(cell_cells)
i=10
j=1
pval_covid_ctrl <- matrix(data = 1,nrow = nrows, ncol = ncols)
for (i in 1:nrows) {
  for (j in 1:ncols) {
    temp <- wilcox.test(x=means_patients_val_covid[i,j,1:covidSamples], y=means_patients_val_ctrl[i,j,1:CtrlSamples], alternative = "two.sided")
    pval_covid_ctrl[i,j] <- temp$p.value
    
  }
}

dim(pval_covid_ctrl)
colnames(pval_covid_ctrl) <- cell_cells
rownames(pval_covid_ctrl) <- intr_pairs
pval_covid_ctrl[1:5,1:5]
# means_patients_val_covid[31:35,31:35,]
# means_patients_val_ctrl[31:35,31:35,]
# pval_covid_ctrl[31:35,31:35]

# pval_covid_ctrl[10]
# pval_covid_ctrl[1:5,1:15]
# pval_covid_ctrl[1:15,1:5]

# checking
# pval_covid_ctrl[5:10,1:5] #[6,2] [10,1]
# means_patients_val_covid[5:10,1:5,]
# means_patients_val_ctrl[5:10,1:5,]

length(which(pval_covid_ctrl<0.05))
length(which(is.na(pval_covid_ctrl)))
temp <- apply(pval_covid_ctrl, 1, min, na.rm=T)
length(which(temp<0.05))
# converting NA to 1, since they came most probably because of the two lists for wilcox test being the same.
temp_pos <- which(is.na(pval_covid_ctrl))
pval_covid_ctrl[temp_pos]
pval_covid_ctrl[temp_pos] <- 1
pval_covid_ctrl[temp_pos]

#*************************************** compute the logFC of LR expression
dim(median_means_patient_val_ctrl)
dim(median_means_patient_val_covid)
# typeof(median_means_patient_val_ctrl)
# typeof(median_means_patient_val_covid)
# converting the median-means to numeric
median_ctrl <- apply(median_means_patient_val_ctrl, 2, as.numeric)
rownames(median_ctrl) <- rownames(median_means_patient_val_ctrl)
median_ctrl[1:5,1:5]
median_covid <- apply(median_means_patient_val_covid, 2, as.numeric) 
rownames(median_covid) <- rownames(median_means_patient_val_covid)
median_covid[1:5,1:5]

# temp <- median_ctrl[order(median_ctrl, decreasing = T)]
# plot(temp)
# temp <- median_covid[order(median_covid, decreasing = T)]
# plot(temp)

# saturating the lowest values to take log
length(which(median_ctrl<0.0001))
length(which(median_covid<0.0001))
median_ctrl[median_ctrl<0.001] <- 0.0001
median_covid[median_covid<0.001] <- 0.0001
median_ctrl <- log2(median_ctrl)
median_covid <- log2(median_covid)

log2FC_median <- median_covid - median_ctrl
temp <- log2FC_median[order(log2FC_median, decreasing = T)]
plot(temp)
length(which(temp>=2))
length(which(temp<=-2))



# saving files
log2FC_median[1:5,1:5]
pval_covid_ctrl[1:5,1:5]
all(rownames(log2FC_median)==rownames(log2FC_median))
all(colnames(log2FC_median)==colnames(log2FC_median))
# write.csv(log2FC_median, file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_logFC_MedianMeanMat.csv")
# write.csv(pval_covid_ctrl, file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pval.csv")


#************************************ FDR correction on all interactions
temp <- as.matrix(pval_covid_ctrl)
dim(temp) <- NULL #flatten a matrix by appending columns
pvalAdj_covid_ctrl <- p.adjust(p=temp,method = "fdr")
pvalAdj_covid_ctrl <- matrix(pvalAdj_covid_ctrl, nrow = nrow(pval_covid_ctrl), ncol = ncol(pval_covid_ctrl), byrow = F)
rownames(pvalAdj_covid_ctrl) <- rownames(pval_covid_ctrl)
colnames(pvalAdj_covid_ctrl) <- colnames(pval_covid_ctrl)
dim(pvalAdj_covid_ctrl)
pvalAdj_covid_ctrl[1:5,1:5]
pval_covid_ctrl[1:5,1:5]

# write.csv(pvalAdj_covid_ctrl, file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pvalAdj_FDR.csv")

#*****************************************************************************************
#************************#getting the top interactions 
#*****************************************************************************************
# median-meanExpr FC
log2FC_median <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_logFC_MedianMeanMat.csv", header = F)
log2FC_median[1:5,1:5]
temp1 <- log2FC_median[1,-1]
head(temp1)
temp2 <- log2FC_median[-1,1]
head(temp2)
log2FC_median <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_logFC_MedianMeanMat.csv",row.names=1)
log2FC_median[1:5,1:5]
colnames(log2FC_median) <- temp1
log2FC_median[1:5,1:5]
dim(log2FC_median)

# PVALUE
pval_covid_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pval.csv", header = F)
pval_covid_ctrl[1:5,1:5]
temp1 <- pval_covid_ctrl[1,-1]
head(temp1)
temp2 <- pval_covid_ctrl[-1,1]
head(temp2)
pval_covid_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pval.csv",row.names=1)
pval_covid_ctrl[1:5,1:5]
colnames(pval_covid_ctrl) <- temp1
pval_covid_ctrl[1:5,1:5]
dim(pval_covid_ctrl)

all(colnames(pval_covid_ctrl)==colnames(log2FC_median))
all(rownames(pval_covid_ctrl)==rownames(log2FC_median))

# PVALUE.ADJ
pvalAdj_covid_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pvalAdj_FDR.csv", header = F)
pvalAdj_covid_ctrl[1:5,1:5]
temp1 <- pvalAdj_covid_ctrl[1,-1]
head(temp1)
temp2 <- pvalAdj_covid_ctrl[-1,1]
head(temp2)
pvalAdj_covid_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_CovidCtrl_pvalAdj_FDR.csv",row.names=1)
pvalAdj_covid_ctrl[1:5,1:5]
colnames(pvalAdj_covid_ctrl) <- temp1
pvalAdj_covid_ctrl[1:5,1:5]
dim(pvalAdj_covid_ctrl)

all(colnames(pvalAdj_covid_ctrl)==colnames(log2FC_median))
all(rownames(pvalAdj_covid_ctrl)==rownames(log2FC_median))

intr_pairs <- rownames(log2FC_median)
cell_cells <- colnames(pval_covid_ctrl)


# getting top interactions
temp_pos1 <- which(log2FC_median>=2)
temp_pos1 <- c(temp_pos1, which(log2FC_median<=-2))
temp_pos2 <- which(pvalAdj_covid_ctrl<=0.1)  #0.05 cell type cutoff; 0.1, all inter cutoff

length(temp_pos1)
length(temp_pos2)
length(intersect(temp_pos1,temp_pos2))

temp_common <- intersect(temp_pos1,temp_pos2)

log2FC_median_req <- matrix(0, nrow = nrow(log2FC_median),ncol = ncol(log2FC_median))
log2FC_median_req[temp_common]= as.matrix(log2FC_median)[temp_common]
pvalAdj_covid_ctrl_req <- matrix(1, nrow = nrow(pvalAdj_covid_ctrl),ncol = ncol(pvalAdj_covid_ctrl))
pvalAdj_covid_ctrl_req[temp_common]= as.matrix(pvalAdj_covid_ctrl)[temp_common]

# finding the top cell-cells and intr-pairs
temp_mat <- matrix(0, nrow = nrow(log2FC_median),ncol = ncol(log2FC_median))
temp_mat[temp_common] <- 1
rownames(temp_mat) <- intr_pairs
colnames(temp_mat) <- cell_cells
temp_pos1 <- which(rowSums(temp_mat)>0)
temp_pos2 <- which(colSums(temp_mat)>0)
temp_inter <- data.frame(inter_pairs = intr_pairs[temp_pos1], nInter=rowSums(temp_mat)[temp_pos1])
temp_cells <- data.frame(cell_cells=cell_cells[temp_pos2], nInter=colSums(temp_mat)[temp_pos2])

head(temp_inter)
top_intr <- temp_inter[order(temp_inter$nInter, decreasing = T),]
head(temp_cells)
top_cells <- temp_cells[order(temp_cells$nInter, decreasing = T),]


#*****************************************************************************************
#************************#plot the top interactions
#*****************************************************************************************

log2FC_median[1:5,1:5]
pvalAdj_covid_ctrl[1:5,1:5]
head(top_cells)
head(top_intr)

# selected_rows <- temp_inter$inter_pairs
temp_pos <- grep("TGF",temp_inter$inter_pairs)
temp_inter$inter_pairs[temp_pos]
temp_pos <- grep("NOTCH",temp_inter$inter_pairs)
temp_inter$inter_pairs[temp_pos]
temp_inter$inter_pairs
selected_rows <- c("FGF1_TGFBR3","EGFR_TGFB1","TGFB1_aVb6 complex","TGFB2_TGFbeta receptor1","TGFB1_TGFbeta receptor2","TGFB3_TGFbeta receptor2","GDF9_TGFR_BMPR2","AXL_IL15RA","AXL_GAS6","ACVR1_BMPR2_BMP6","ACVR1_BMPR2_BMP8A","FGF10_FGFR1","FGF10_FGFR2","BMP2_SMO","BMPR1A_BMPR2_BMP2","BMR1A_ACR2A_BMP2")
temp_cells$cell_cells
temp_pos <- grep("Epithelial",temp_cells$cell_cells)
temp_cells$cell_cells[temp_pos]
temp_pos <- c(temp_pos,grep("Fibroblasts",temp_cells$cell_cells))
temp_cells$cell_cells[temp_pos]
selected_columns <- temp_cells$cell_cells[temp_pos]


sel_pval = pvalAdj_covid_ctrl[match(selected_rows, rownames(pvalAdj_covid_ctrl)), selected_columns]
head(sel_pval)
dim(sel_pval)
# temp_sel_pval = all_pval_ctrl[match(selected_rows, intr_pairs), selected_columns]
# head(temp_sel_pval)

sel_means = log2FC_median[match(selected_rows, rownames(log2FC_median)), selected_columns]
head(sel_means)
dim(sel_means)

all(colnames(sel_means)==colnames(sel_pval))
# match(colnames(sel_means),colnames(sel_pval))
# grep("Fibroblasts.1",colnames(sel_means))
# gsub("Fibroblasts.1","Fibroblasts",colnames(sel_means))
# colnames(sel_means) <- gsub("Fibroblasts.1","Fibroblasts",colnames(sel_means))
all(colnames(sel_means)==colnames(sel_pval))

all(rownames(sel_means)==rownames(sel_pval))

df_names = expand.grid(selected_rows, selected_columns)#create a data frame of all combinations of the supplied vectors
dim(df_names)
head(df_names)
pval = as.numeric(unlist(sel_pval))
plot(pval)
plot(-log10(pval))
pval[pval<0.0009] = 0.0009
plot.data = cbind(df_names,pval)
head(plot.data)
pr = as.numeric(unlist(sel_means))
plot(pr)
# pr[pr==0] = 1
# plot.data = cbind(plot.data,log2(pr))
plot.data = cbind(plot.data,pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)

plot.data$pair <- as.character(plot.data$pair)
plot.data$clusters <- as.character(plot.data$clusters)


#interacting pairs
unique(plot.data[,"pair"])

#correcting the ligand receptor interactions
# forCorrecting <- c("MET_HGF")
# correct <-       c("HGF_MET")
# correctionList <- data.frame(forCorrecting=forCorrecting, correct=correct, stringsAsFactors = F)
# 
# 
# # ligand_list <- reqGenes
# i <- 1
# plot.data_corrected <- plot.data
# for (i in 1:nrow(correctionList)) {
#   temp_pos <- grep(correctionList[i,"forCorrecting"], plot.data[,"pair"])
#   plot.data[temp_pos,]
#   j=temp_pos[1]
#   for (j in temp_pos) {
#     plot.data[j,]
#     # temp <- plot.data[j,]
#     # temp
#     correctionList[i,]
#     plot.data_corrected[j,"pair"] <- correctionList[i,"correct"] #correcting the interaction
#     #correcting the clusters
#     # temp_clusters <- plot.data[j,"clusters"]
#     temp_clusters <- strsplit(plot.data[j,"clusters"], split = "\\|")[[1]]
#     plot.data_corrected[j,"clusters"] <- paste0(temp_clusters[2],"|",temp_clusters[1])
#     # plot.data[j,]
#     # plot.data_corrected[j,]
#   }
#   plot.data_corrected[temp_pos,]
# }
plot.data_corrected <- plot.data
head(plot.data_corrected)
unique(plot.data_corrected[,"clusters"])

# #keeping only xxx as the source of the ligands
# temp_pos <- grep("xxx\\|", plot.data_corrected[,"clusters"])
# plot.data_corrected[temp_pos,"clusters"]
# plot.data_corrected[-temp_pos,"clusters"]
# plot.data_corrected_req <- plot.data_corrected[temp_pos,]
# plot.data_corrected_req[,"clusters"]
# 
# #keeping only HGF_MET interactions
# temp_pos <- grep("HGF_MET", plot.data_corrected_req[,"pair"])
# plot.data_corrected_req[temp_pos,]
# plot.data_corrected_req[-temp_pos,]
# plot.data_corrected_req <- plot.data_corrected_req[temp_pos,]
# plot.data_corrected_req[,"clusters"]
plot.data_corrected_req <- plot.data_corrected


my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
ggplot(plot.data_corrected_req,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  # scale_size(breaks = c(0,1,2,3)) +
  scale_color_gradientn(colors=my_palette) + #, limits=c(-2,2)
  # lims(colour = c(-2,2)) +
  labs(size="-log10(pvalue.Adj)", color="Log2FC (COVID vs Ctrl)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
filename = 'cellphoneDB/analysis/dotplot_covidVsctrl_topSel.pdf'
width = 15
height = 10
# ggsave(filename,width = width, height = height, limitsize=F)#width = width, height = height,,  device = cairo_pdf


#*****************************************************************************************
#************************#plot the selected interactions
#*****************************************************************************************

log2FC_median[1:5,1:5]
pvalAdj_covid_ctrl[1:5,1:5]
head(temp_cells)
head(temp_inter)

reqGenes <- read.table("data/cpdb_genesReq.txt", quote="", comment.char="")
head(reqGenes)
dim(reqGenes)
reqGenes <- reqGenes[,"V1"] #getting to a character vector
# reqGenes <- c(reqGenes, "HGF","CD47")
reqGenes

#selecting the rows
temp_pos <- unlist(sapply(reqGenes, function(x, cpdb_list){
  c(grep(paste0(x,"[[:punct:]]"), cpdb_list), grep(paste0("[[:punct:]]",x,"$"), cpdb_list), grep(paste0("[[:punct:]]",x,"[[:punct:]]"), cpdb_list)) #this will avoid homolog names eg: NGF, NGFR
},rownames(log2FC_median)))
# all_means_ctrl[,"interacting_pair"]))
# all_means_ctrl[unique(temp_pos),]
temp_pos
length(temp_pos)
length(unique(temp_pos))
temp_pos_req <- unique(temp_pos)
rownames(log2FC_median)[temp_pos_req]
log2FC_median[temp_pos_req,]
selected_rows = rownames(log2FC_median)[temp_pos_req]
# selected_rows <- temp_inter$inter_pairs

selected_columns <- temp_cells$cell_cells


sel_pval = pvalAdj_covid_ctrl[match(selected_rows, rownames(pvalAdj_covid_ctrl)), selected_columns]
head(sel_pval)
dim(sel_pval)
# temp_sel_pval = all_pval_ctrl[match(selected_rows, intr_pairs), selected_columns]
# head(temp_sel_pval)

sel_means = log2FC_median[match(selected_rows, rownames(log2FC_median)), selected_columns]
head(sel_means)
dim(sel_means)

all(colnames(sel_means)==colnames(sel_pval))
all(rownames(sel_means)==rownames(sel_pval))

df_names = expand.grid(selected_rows, selected_columns)#create a data frame of all combinations of the supplied vectors
dim(df_names)
head(df_names)
# df_names[which(df_names[,"Var1"]=="MET_HGF"),]
pval = as.numeric(unlist(sel_pval))
plot(pval)
plot(-log10(pval))
pval[pval<0.0009] = 0.0009
plot.data = cbind(df_names,pval)
head(plot.data)
pr = as.numeric(unlist(sel_means))
plot(pr)
# pr[pr==0] = 1
# plot.data = cbind(plot.data,log2(pr))
plot.data = cbind(plot.data,pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)

plot.data$pair <- as.character(plot.data$pair)
plot.data$clusters <- as.character(plot.data$clusters)


#interacting pairs
unique(plot.data[,"pair"])

#correcting the ligand receptor interactions
# forCorrecting <- c("MET_HGF")
# correct <-       c("HGF_MET")
# correctionList <- data.frame(forCorrecting=forCorrecting, correct=correct, stringsAsFactors = F)
# 
# 
# # ligand_list <- reqGenes
# i <- 1
# plot.data_corrected <- plot.data
# for (i in 1:nrow(correctionList)) {
#   temp_pos <- grep(correctionList[i,"forCorrecting"], plot.data[,"pair"])
#   plot.data[temp_pos,]
#   j=temp_pos[1]
#   for (j in temp_pos) {
#     plot.data[j,]
#     # temp <- plot.data[j,]
#     # temp
#     correctionList[i,]
#     plot.data_corrected[j,"pair"] <- correctionList[i,"correct"] #correcting the interaction
#     #correcting the clusters
#     # temp_clusters <- plot.data[j,"clusters"]
#     temp_clusters <- strsplit(plot.data[j,"clusters"], split = "\\|")[[1]]
#     plot.data_corrected[j,"clusters"] <- paste0(temp_clusters[2],"|",temp_clusters[1])
#     # plot.data[j,]
#     # plot.data_corrected[j,]
#   }
#   plot.data_corrected[temp_pos,]
# }
plot.data_corrected <- plot.data
head(plot.data_corrected)
unique(plot.data_corrected[,"clusters"])

# #keeping only xxx as the source of the ligands
# temp_pos <- grep("xxx\\|", plot.data_corrected[,"clusters"])
# plot.data_corrected[temp_pos,"clusters"]
# plot.data_corrected[-temp_pos,"clusters"]
# plot.data_corrected_req <- plot.data_corrected[temp_pos,]
# plot.data_corrected_req[,"clusters"]
# 
# #keeping only HGF_MET interactions
# temp_pos <- grep("HGF_MET", plot.data_corrected_req[,"pair"])
# plot.data_corrected_req[temp_pos,]
# plot.data_corrected_req[-temp_pos,]
# plot.data_corrected_req <- plot.data_corrected_req[temp_pos,]
# plot.data_corrected_req[,"clusters"]
plot.data_corrected_req <- plot.data_corrected


my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
ggplot(plot.data_corrected_req,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  # scale_size(breaks = c(0,1,2,3)) +
  scale_color_gradientn(colors=my_palette) + #, limits=c(-2,2)
  # lims(colour = c(-2,2)) +
  labs(size="-log10(pvalue.Adj)", color="Log2FC (COVID vs Ctrl)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
filename = 'cellphoneDB/analysis/dotplot_covidVsctrl_selected.pdf'
width = 15
height = 10
ggsave(filename,width = width, height = height, limitsize=F)#width = width, height = height,,  device = cairo_pdf

