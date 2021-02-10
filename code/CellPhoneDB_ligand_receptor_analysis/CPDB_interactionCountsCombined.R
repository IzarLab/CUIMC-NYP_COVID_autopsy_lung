#!/usr/bin/env Rscript


# script for analyzing and plotting the interaction
# counts between cell types of different samples
# obtained by CellPhoneDB run in each sample


library(pheatmap)
library(igraph)
library(ggpubr)
library(ggplot2)


# *****************************************************************************************
# getting Ctrl interaction-counts
# *****************************************************************************************
# get the count_network from each individual sample
# and combined to a single file
path1 <- c("cellphoneDB/Ctrl")
path2 <- c("out_cellTypeMain/count_network.txt")
dirs_available <- list.dirs(path = path1, full.names = F, 
                            recursive = F)
dirs_available

# patient <- dirs_available[1]
counts_mat <- NULL
for (patient in dirs_available) {
  print(patient)
  temp_file <- file.path(path1, patient, path2)
  print(temp_file)
  if (file.exists(temp_file)) {
    count_network <- read.delim(temp_file, row.names = NULL, 
                                stringsAsFactors = FALSE)
    # head(count_network) count_network <-
    # count_network
    count_network$source_target <- paste0(count_network$SOURCE, 
                                          "_", count_network$TARGET)  #combining the source and target names for sorting
    count_network$source_target <- gsub(" ", "-", 
                                        count_network$source_target)  #removing spaces in the names for sorting and comparison
    count_network <- count_network[order(count_network$source_target), 
    ]  #sorting to make the order in different patients match
    head(count_network)
    temp_counts <- count_network[, c("count", "source_target")]  #getting only the required columns
    colnames(temp_counts) <- paste0(patient, "-", 
                                    colnames(temp_counts))
    # head(temp_counts)
    if (is.null(counts_mat)) {
      counts_mat <- cbind(count_network[, c("SOURCE", 
                                            "TARGET")], temp_counts)
    } else {
      if (all(counts_mat[, ncol(counts_mat)] == 
              temp_counts[, 2])) {
        counts_mat <- cbind(counts_mat, temp_counts)
      } else {
        print(paste0("Error: interaction mismatch in ", 
                     patient))
      }
      
    }
    # head(counts_mat)
  } else {
    print(paste0("File not found: ", temp_file))
  }
  
}
head(counts_mat)
# keeping only the counts and other required
# columns
counts_mat <- counts_mat[, c(1, 2, grep("count", colnames(counts_mat)))]
head(counts_mat)
# write.csv(counts_mat,file =
# 'cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Ctrl.csv')

# *****************************************************************************************
# getting Covid interaction-counts
# *****************************************************************************************
# get the count_network from each individual sample
# and combined to a single file
path1 <- c("cellphoneDB/Covid")
path2 <- c("out_cellTypeMain/count_network.txt")
dirs_available <- list.dirs(path = path1, full.names = F, 
                            recursive = F)
dirs_available

# patient <- 'L17cov' #'L22cov' 'L17cov'
# dirs_available[1]
counts_mat <- NULL
for (patient in dirs_available) {
  print(patient)
  temp_file <- file.path(path1, patient, path2)
  print(temp_file)
  if (file.exists(temp_file)) {
    count_network <- read.delim(temp_file, row.names = NULL, 
                                stringsAsFactors = FALSE)
    head(count_network)
    # count_network <- count_network
    count_network$source_target <- paste0(count_network$SOURCE, 
                                          "_", count_network$TARGET)  #combining the source and traget names for sorting
    count_network$source_target <- gsub(" ", "-", 
                                        count_network$source_target)  #removing spaces in the names for sorting and comparison
    count_network <- count_network[order(count_network$source_target), 
    ]  #sorting to make the order in different patients match
    head(count_network)  #
    temp_counts <- count_network[, c("count", "source_target")]  #getting only the required columns
    colnames(temp_counts) <- paste0(patient, "-", 
                                    colnames(temp_counts))
    head(temp_counts)  #
    if (is.null(counts_mat)) {
      counts_mat <- cbind(count_network[, c("SOURCE", 
                                            "TARGET")], temp_counts)
    } else {
      if (all(counts_mat[, ncol(counts_mat)] == 
              temp_counts[, 2])) {
        counts_mat <- cbind(counts_mat, temp_counts)
      } else {
        print(paste0("Warning: interaction mismatch in ", 
                     patient))
        # temp_cell_cell <-
        # setdiff(counts_mat[,ncol(counts_mat)],temp_counts[,2])
        # temp_nullCounts <-
        # data.frame(count=rep(0,length(temp_cell_cell)),
        # temp_cell_cell) colnames(temp_nullCounts) <-
        # colnames(temp_counts) head(temp_nullCounts)
        temp_newCounts <- counts_mat[, c(ncol(counts_mat) - 
                                           1, ncol(counts_mat))]
        colnames(temp_newCounts) <- colnames(temp_counts)
        temp_newCounts[, 1] <- rep(0, nrow(temp_newCounts))
        head(temp_newCounts)
        temp_pos <- match(temp_counts[, 2], 
                          temp_newCounts[, 2])
        temp_newCounts[temp_pos, ] <- temp_counts[, 
        ]
        # # check the missing interactions temp <-
        # setdiff(seq(1,nrow(temp_newCounts)),temp_pos)
        # temp_newCounts[temp,2]
        if (all(counts_mat[, ncol(counts_mat)] == 
                temp_newCounts[, 2])) {
          print(paste0("Update: interaction mismatch corrected with zero counts in ", 
                       patient))
          counts_mat <- cbind(counts_mat, temp_newCounts)
        } else {
          print(paste0("Error: interaction mismatch could not be corrected in ", 
                       patient))
        }
      }
      
    }
    head(counts_mat)  #
    # 
  } else {
    print(paste0("File not found: ", temp_file))
  }
  
}
head(counts_mat)
# keeping only the counts and other required
# columns
counts_mat <- counts_mat[, c(1, 2, grep("count", colnames(counts_mat)))]
head(counts_mat)
# write.csv(counts_mat,file =
# 'cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Covid.csv')


# *****************************************************************************************
# find the median counts for COVID and Control samples separately
# *****************************************************************************************
# *CONTROL
counts_mat_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Ctrl.csv", 
                            row.names = 1)
head(counts_mat_ctrl)

temp <- apply(counts_mat_ctrl[, 3:ncol(counts_mat_ctrl)], 
              MARGIN = 1, FUN = median)
counts_mat_ctrl_med <- cbind(counts_mat_ctrl[, 1:2], 
                             temp)
colnames(counts_mat_ctrl_med)[3] <- "count"
head(counts_mat_ctrl_med)

count_matrix_ctrl = matrix(counts_mat_ctrl_med[, "count"], 
                           nrow = length(unique(counts_mat_ctrl_med[, "SOURCE"])), 
                           ncol = length(unique(counts_mat_ctrl_med[, "TARGET"])), 
                           byrow = T)
rownames(count_matrix_ctrl) = unique(counts_mat_ctrl_med[, 
                                                         "SOURCE"])
colnames(count_matrix_ctrl) = unique(counts_mat_ctrl_med[, 
                                                         "TARGET"])
count_matrix_ctrl
# write.csv(count_matrix_ctrl,file =
# 'cellphoneDB/analysis/ligandReceptors_Counts_median_Ctrl.csv')


# COVID
counts_mat_covid <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Covid.csv", 
                             row.names = 1)
head(counts_mat_covid)
temp <- apply(counts_mat_covid[, 3:ncol(counts_mat_covid)], 
              MARGIN = 1, FUN = median)
counts_mat_covid_med <- cbind(counts_mat_covid[, 1:2], 
                              temp)
colnames(counts_mat_covid_med)[3] <- "count"
head(counts_mat_covid_med)

count_matrix_covid = matrix(counts_mat_covid_med[, 
                                                 "count"], nrow = length(unique(counts_mat_covid_med[, 
                                                                                                     "SOURCE"])), ncol = length(unique(counts_mat_covid_med[, 
                                                                                                                                                            "TARGET"])), byrow = T)
rownames(count_matrix_covid) = unique(counts_mat_covid_med[, 
                                                           "SOURCE"])
colnames(count_matrix_covid) = unique(counts_mat_covid_med[, 
                                                           "TARGET"])
count_matrix_covid
# write.csv(count_matrix_covid,file =
# 'cellphoneDB/analysis/ligandReceptors_Counts_median_Covid.csv')


# *****************************************************************************************
# #creating cell-cell interactome
# *****************************************************************************************

count_matrix_ctrl <- as.matrix(read.csv(file = "cellphoneDB/analysis/ligandReceptors_Counts_median_Ctrl.csv", 
                                        row.names = 1))
head(count_matrix_ctrl)

count_matrix_covid <- as.matrix(read.csv(file = "cellphoneDB/analysis/ligandReceptors_Counts_median_Covid.csv", 
                                         row.names = 1))
head(count_matrix_covid)


# get cellTypeMain freq
prop_cellTypeMain <- read.csv("data/ligandReceptors_cellType_proportions_main.csv")  #boxplot_proportions_main.csv
head(prop_cellTypeMain, 10)
length(unique(prop_cellTypeMain$orig.ident))
table(prop_cellTypeMain$cell_type_main)
dim(prop_cellTypeMain)
# get control proportions
temp_pos <- grep("C", prop_cellTypeMain$orig.ident)
table(prop_cellTypeMain$orig.ident[temp_pos])
table(prop_cellTypeMain$cell_type_main[temp_pos])
length(unique(prop_cellTypeMain$orig.ident[temp_pos]))
freq_ctrl <- matrix(prop_cellTypeMain$freq[temp_pos], 
                    nrow = 9, ncol = 7, byrow = F)
rownames(freq_ctrl) <- prop_cellTypeMain$cell_type_main[1:9]
unique(prop_cellTypeMain$orig.ident[temp_pos])
colnames(freq_ctrl) <- unique(prop_cellTypeMain$orig.ident[temp_pos])
head(freq_ctrl)
colSums(freq_ctrl)

# get covid proportions
temp_pos <- grep("L", prop_cellTypeMain$orig.ident)
table(prop_cellTypeMain$orig.ident[temp_pos])
table(prop_cellTypeMain$cell_type_main[temp_pos])
length(unique(prop_cellTypeMain$orig.ident[temp_pos]))
freq_covid <- matrix(prop_cellTypeMain$freq[temp_pos], 
                     nrow = 9, ncol = 20, byrow = F)
rownames(freq_covid) <- prop_cellTypeMain$cell_type_main[1:9]
unique(prop_cellTypeMain$orig.ident[temp_pos])
colnames(freq_covid) <- unique(prop_cellTypeMain$orig.ident[temp_pos])
head(freq_covid)
colSums(freq_covid)

# *********************** create a igraph object
# from counts matrix 
#****************** control
net <- graph_from_adjacency_matrix(count_matrix_ctrl, 
                                   mode = "undirected", weighted = T, diag = T)
net
# layout
layOut = layout_in_circle(net)  #get the circular layout
# edge weight
E(net)$width <- E(net)$weight/5
# # node size
nodeSize <- apply(freq_ctrl, MARGIN = 1, FUN = median)
V(net)$size <- nodeSize * 100

# saving the figures
pdf("results/counts_interactome_ctrl.pdf")
plot(net, layout = layOut, vertex.shape = "circle", 
     vertex.label.dist = -2, edge.color = "#529EFF", 
     vertex.color = "#A3AA00")  #, edge.curved=TRUE)
title("Control")
dev.off()
dev.off()


## ****************** covid
net <- graph_from_adjacency_matrix(count_matrix_covid, 
                                   mode = "undirected", weighted = T, diag = T)
net
# layout
layOut = layout_in_circle(net)  #get the circular layout
# edge weight
E(net)$width <- E(net)$weight/5
# # node size
nodeSize <- apply(freq_covid, MARGIN = 1, FUN = median)
V(net)$size <- nodeSize * 100

pdf("results/counts_interactome_covid.pdf")
plot(net, layout = layOut, vertex.shape = "circle", 
     vertex.label.dist = -2, edge.color = "#f37735", 
     vertex.color = "#A3AA00")
title("Covid")
dev.off()
dev.off()



# *****************************************************************************************
# #analysis of interactions between major cell types
# *****************************************************************************************
# *CONTROL
counts_mat_ctrl <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Ctrl.csv", 
                            row.names = 1)
head(counts_mat_ctrl)
dim(counts_mat_ctrl)
# temp_pos <- which(counts_mat_ctrl$SOURCE %in%
# c('Fibroblasts','Epithelial cells','Endothelial
# cells'))
counts_mat_ctrl_req <- counts_mat_ctrl  #[temp_pos,]
counts_mat_ctrl_req[1:5, 1:5]

temp <- counts_mat_ctrl_req[, 3:ncol(counts_mat_ctrl)]
rownames(temp) <- paste0(counts_mat_ctrl_req[, 1], 
                         "|", counts_mat_ctrl_req[, 2])
temp[1:5, 1:5]  #<- as.data.frame(temp)

# converting to long-form for plotting
temp_plot <- stack(temp)
temp_plot$cell_cell <- rep(rownames(temp), ncol(temp))
head(temp_plot)

ggplot(temp_plot, aes(x = cell_cell, y = values)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1))



# COVID
counts_mat_covid <- read.csv(file = "cellphoneDB/analysis/ligandReceptors_combinedCountsMat_Covid.csv", 
                             row.names = 1)
head(counts_mat_covid)

dim(counts_mat_covid)
# temp_pos <- which(counts_mat_covid$SOURCE %in%
# c('Fibroblasts','Epithelial cells','Endothelial
# cells'))
counts_mat_covid_req <- counts_mat_covid  #[temp_pos,]
counts_mat_covid_req[1:5, 1:5]

temp <- counts_mat_covid_req[, 3:ncol(counts_mat_covid)]
rownames(temp) <- paste0(counts_mat_covid_req[, 1], 
                         "|", counts_mat_covid_req[, 2])
temp[1:5, 1:5]  #<- as.data.frame(temp)

temp1 <- stack(temp)
temp1$cell_cell <- rep(rownames(temp), ncol(temp))
head(temp1)
ggplot(temp1, aes(x = cell_cell, y = values)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

temp_plot <- rbind(temp_plot, temp1)
table(temp_plot$ind)

temp_plot$group <- substr(temp_plot$ind, start = 1, 
                          stop = 1)
temp_plot$group <- gsub("C", "Control", temp_plot$group)
temp_plot$group <- gsub("L", "COVID", temp_plot$group)
# plotting for every possible cell-cell combination
ggplot(temp_plot, aes(x = cell_cell, y = values, fill = group)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1)) + scale_fill_manual(values = c("#005AC8", 
                                                                                               "#AA0A3C")) + ylab("Number of interactions") + 
  stat_compare_means(method = "t.test", aes(label = ..p.signif.., 
                                            group = group))


# plotting for all interactions from each cell type
temp <- c("Fibroblasts", "Epithelial cells", "Endothelial cells", 
          "Mast cell-like", "Myeloid", "Neuron-like cells", 
          "B cells", "APC-like", "T cells")
i = temp[2]
temp_plot$cell <- NA
for (i in temp) {
  # print(i)
  temp_pos <- grep(paste0(i, "\\|"), temp_plot$cell_cell)
  table(temp_plot$cell_cell[temp_pos])
  temp_plot$cell[temp_pos] <- i
}
ggplot(temp_plot, aes(x = cell, y = values, fill = group)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1)) + scale_fill_manual(values = c("#006E82", 
                                                                                               "#AA0A3C")) + ylab("Number of interactions") + 
  geom_dotplot(binaxis = "y", binwidth = 0.5, dotsize = 0.5, 
               stackdir = "center", position = position_dodge(1)) + 
  stat_compare_means(method = "t.test", aes(label = ..p.format.., 
                                            group = group))  #abel = ..p.signif..

ggplot(temp_plot, aes(x = cell, y = values, color = group)) + 
  geom_boxplot(lwd = 1, fatten = 1) + theme(axis.text.x = element_text(angle = 90, 
                                                                       hjust = 1)) + scale_color_manual(values = c("#006E82", 
                                                                                                                   "#AA0A3C")) + ylab("Number of interactions") + 
  geom_dotplot(binaxis = "y", binwidth = 0.5, dotsize = 1, 
               stackdir = "centerwhole", position = position_dodge(0.75)) + 
  stat_compare_means(method = "t.test", size = 3, 
                     label.y = 150, aes(label = paste0("p=", ..p.format..), 
                                        group = group))  #abel = ..p.signif..
ggsave(filename = "manuscript/figs/CellPhoneDB/boxplot_dotplot_counts_ctrlCovid.pdf")