working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"
source("8_Loadsoftware_libraries.R")

combined = readRDS(paste0(workingdirectory,"/lungs_all_yiping.rds"))
userevisedtuft = TRUE

source("13_All_samples_UMAP.R")

combined = combined[,!is.na(combined$celltype_bped_main)]

source("11_Add_gene_signature_module_scores.R")

embed = Embeddings(combined,reduction="umap")
pdf(file=paste0(workingdirectory,"/Figure_1cd"),height=7,width=14)
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
tplot = DimPlot(combined, reduction = "umap",label = T,group.by = 'overallclassification',repel = T,label.size = 2.5)
tplot[[1]]$layers[[1]]$aes_params$alpha = .5
print(tplot + 
  ggtitle('overallclassification') +
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=7)),
  vp=vplayout(1,1))

tplot2 = DimPlot(combined, reduction = "umap",label = F,group.by = 'group',repel=T)
tplot2[[1]]$layers[[1]]$aes_params$alpha = .2
print(tplot2+ggtitle('patientgroup'),vp=vplayout(1,2))
dev.off()

saveRDS(combined,workingdirectory("/combinedwithtwoizarcontrols.rds"))