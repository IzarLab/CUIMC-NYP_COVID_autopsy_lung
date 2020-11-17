library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)
library(grid)
library(dplyr)
library(ggsignif)
library(ggpubr)
library(pheatmap)
library(SingleCellExperiment)
library(SingleR)

working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"
covid_sigs<-read.csv(paste0(working_directory,'/covid_signatures_reduced.csv'),na.strings = '')
for(c in 1:ncol(covid_sigs)){
  covid_sigs[,c] = toupper(covid_sigs[,c])
}
