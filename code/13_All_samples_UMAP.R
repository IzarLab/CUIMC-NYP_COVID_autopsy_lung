working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

#Manually broad celltype categories, based on visual inspection of UMAP of Seurat clusters
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(0,1,5,9,13,16,17,22),'Myeloid','NA')
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(10,12,18,23,25),'Epithelial cells',combined@meta.data$cell_type_main)
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(7,21),'Endothelial cells',combined@meta.data$cell_type_main)
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(4,8,11,15,20),'Fibroblasts',combined@meta.data$cell_type_main)
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(19),'Mast cells',combined@meta.data$cell_type_main)
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(6,24),'B cells',combined@meta.data$cell_type_main)
combined@meta.data$cell_type_main<-ifelse(combined@meta.data$integrated_snn_res.0.8 %in% c(2,3,14),'T cells',combined@meta.data$cell_type_main)

pdf(paste0(workingdirectory,"/clustermap.pdf"))
print(DimPlot(combined, reduction = "umap",label = F,group.by = 'cell_type_main'))
print(DimPlot(combined, reduction = "umap",group.by = 'ident',label=T,repel=T,label.size=2.5))
dev.off()
