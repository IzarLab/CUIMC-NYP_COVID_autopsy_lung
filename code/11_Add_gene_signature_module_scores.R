working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

#immune cell types are taken from SingleR's bped database into overallclassification
#other cells are initially classified as other
immunetypes_bped_main = c("Macrophages","B-cells","CD8+ T-cells","NK cells","CD4+ T-cells","Monocytes","DC","Eosinophils","Neutrophils")
combined$immunestatus = (combined$celltype_bped_main %in% immunetypes_bped_main)
combined$overallclassification = "Other"
combined$overallclassification[combined$immunestatus] = combined$celltype_bped_main[combined$immunestatus]
combined$overallclassification[combined$celltype_bped_fine=="Plasma cells"] = "Plasma cells"

#signatures for covid analysis are read in and used to assign module scores
#for AT1 and AT2, only cells previously assigned to the Epithelial cells cluster are allowed to have a positive module score
# list of all assigned module scores is stored in columnscoreslist
covid_sigs<-read.csv('covid_signatures_reduced.csv',na.strings = '')
columnscoreslist = list()
for(c in 1:ncol(covid_sigs)){
  covid_sigs[,c] = toupper(covid_sigs[,c])
  combined = AddModuleScore(object = combined, features = list(na.omit(covid_sigs[,c])), name = colnames(covid_sigs)[c], assay = 'RNA', search = T)
  if (names(covid_sigs)[c]=="AT1")
  {
    combined$AT11[combined$cell_type_main!="Epithelial cells"] = 0
  }
  if (names(covid_sigs)[c]=="AT2")
  {
    combined$AT21[combined$cell_type_main!="Epithelial cells"] = 0
  }
  columnscores = eval(parse(text=paste0("combined$",colnames(covid_sigs)[c],"1")))
  columnscoreslist[[c]] = columnscores
}

source("12_Overall_celltype_classification.R")

#redefine a tuft cell signature, based on cells whose original tuft cell signature scores were > 1.0, and were located in cluster 20 of the seurat clustering
#rerun cell type assignment in overallclassification, using the new tuft signature
if (userevisedtuft)
{

  #combined = SetIdent(combined, cells=colnames(combined)[combined$tuft_cell1>1.0 & (combined$overallclassification %in% c("tuft_cell","tuft_1_montoro","tuft_2_montoro")) & combined$cell_type_main=="Fibroblasts"], value="high_tuft_cell_score")
  combined = SetIdent(combined, cells=colnames(combined)[combined$tuft_cell1>1.0 & (combined$overallclassification %in% c("tuft_cell","tuft_1_montoro","tuft_2_montoro")) & combined@meta.data$integrated_snn_res.0.8==20], value="high_tuft_cell_score")
  high_tuft_markers = FindMarkers(combined,ident.1="high_tuft_cell_score",ident.2=NULL)
  high_tuft_markers_short = rownames(high_tuft_markers)[1:30]
  combinedold = combined

  immunetypes_bped_main = c("Macrophages","B-cells","CD8+ T-cells","NK cells","CD4+ T-cells","Monocytes","DC","Eosinophils","Neutrophils")
  combined$immunestatus = (combined$celltype_bped_main %in% immunetypes_bped_main)
  combined$overallclassification = "Other"
  combined$overallclassification[combined$immunestatus] = combined$celltype_bped_main[combined$immunestatus]
  combined$overallclassification[combined$celltype_bped_fine=="Plasma cells"] = "Plasma cells"

  covid_sigs<-read.csv('covid_signatures_reduced.csv',na.strings = '')
  columnscoreslist = list()
  covid_sigs$tuft_cell = NULL
  covid_sigs$tuft_1_montoro = NULL
  covid_sigs$tuft_2_montoro = NULL
  covid_sigs$high_tuft_markers_short = c(high_tuft_markers_short,rep(NA,dim(covid_sigs)[1]-length(high_tuft_markers_short)))

  for(c in 1:ncol(covid_sigs)){
    covid_sigs[,c] = toupper(covid_sigs[,c])
    combined = AddModuleScore(object = combined, features = list(na.omit(covid_sigs[,c])), name = colnames(covid_sigs)[c], assay = 'RNA', search = T)
    if (names(covid_sigs)[c]=="AT1")
    {
      combined$AT11[combined$cell_type_main!="Epithelial cells"] = 0
    }
    if (names(covid_sigs)[c]=="AT2")
    {
      combined$AT21[combined$cell_type_main!="Epithelial cells"] = 0
    }
    if (userevisedtuft && names(covid_sigs)[c]=="high_tuft_markers_short")
    {
      combined$high_tuft_markers_short1[combined$cell_type_main!="Fibroblasts"] = 0
      combined$high_tuft_markers_short1[!(combinedold$overallclassification %in% c("tuft_cell","tuft_1_montoro","tuft_2_montoro"))] = 0
      combined$high_tuft_markers_short1[combined$high_tuft_markers_short1<0.7] = 0
    }
    columnscores = eval(parse(text=paste0("combined$",colnames(covid_sigs)[c],"1")))
    columnscoreslist[[c]] = columnscores#[!combined$immunestatus]
  }
  source("12_Overall_celltype_classification.R")
}
