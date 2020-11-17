#!/usr/bin/env Rscript
print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(purrr)
library(cowplot)
library(stringr)

working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

pat_list<-c('lungs_02_cov', 'lungs_03_cov', 'lungs_05_cov', 'lungs_06_cov', 'lungs_08_cov', 
            'lungs_11_cov', 'lungs_12_cov', 'lungs_14_ctr', 'lungs_15_ctr')

# load objects
for(pat in pat_list){
  assign(pat,readRDS(paste0(pat,'/data_',pat,'_cb.rds')))
}

# create object.list
object.list<-NULL
for(obj in grep('lungs_',ls(),value=T)){
  object.list<-c(object.list,eval(parse(text = obj)))
}

# find anchors
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:20)#, k.filter = 70)

# integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:20)#, k.weight = 70)

# normal workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu)
seu <- RunUMAP(object = seu, dims = 1:20)

DefaultAssay(seu) <- "RNA"
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- "integrated"

#markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) %>% write.csv('data/lungs_all/markers_lungs_all.csv',row.names = F)

saveRDS(seu, file = paste0(working_directory,'/lungs_all_yiping.rds'))
