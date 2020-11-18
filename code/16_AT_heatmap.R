working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library(pheatmap)

combined = readRDS(paste0(working_directory,"/combinedwithtwoizarcontrols.rds"))
DefaultAssay(combined) = "RNA"

#read in DATP signatures
DATP_sigs<-read.csv('DATP_sigs.csv',na.strings = '')
for(c in 1:ncol(DATP_sigs)){
  DATP_sigs[,c] = toupper(DATP_sigs[,c])
  combined = AddModuleScore(object = combined, features = list(na.omit(DATP_sigs[,c])), name = colnames(DATP_sigs)[c], assay = 'RNA', search = T)
}

#make new ATstatus field, and assign values based on overallclassification
combined$ATstatus = "None"
combined$ATstatus[combined$overallclassification=="AT1"] = "AT1"
combined$ATstatus[combined$overallclassification=="AT2"] = "AT2"
#find markers for AT1 vs. AT2 status
combined = SetIdent(combined, cells=colnames(combined)[combined$ATstatus=="AT1"],value="AT1")
combined = SetIdent(combined, cells=colnames(combined)[combined$ATstatus=="AT2"],value="AT2")
atmarkers = FindMarkers(combined,ident.1="AT1",ident.2="AT2")

covid_sigs_all = c(na.omit(covid_sigs$AT1),na.omit(covid_sigs$AT2))
topatmarkers = rownames(atmarkers)[1:30]
write.table(topatmarkers,paste0(workingdirectory,"Extended_Data_Table_5.csv"),row.names=T,col.names=T,sep=",",quote=F)
topatmarkers = topatmarkers[!(topatmarkers %in% covid_sigs_all)]
covid_sigs_all = c(covid_sigs_all,topatmarkers)

#subset data based on cells having AT classification, and located in the Epithelial cells cluster
combinedsmall = combined[,combined$ATstatus!="None" & combined$cell_type_main=="Epithelial cells"]
DefaultAssay(combinedsmall) <- "RNA"

#make four new classes for AT1 and 2 status, combined with whether UP_in_DATPs signature is above or below the median
combinedsmall$AT_DATP_status = "None"
median1 = median(combinedsmall$UP_in_DATPs1[combinedsmall$ATstatus=="AT1"])
median2 = median(combinedsmall$UP_in_DATPs1[combinedsmall$ATstatus=="AT2"])
combinedsmall$AT_DATP_status[combinedsmall$ATstatus=="AT1" & combinedsmall$UP_in_DATPs1<median1] = "AT1_low_UP_in_DATPs"
combinedsmall$AT_DATP_status[combinedsmall$ATstatus=="AT1" & combinedsmall$UP_in_DATPs1>median1] = "AT1_high_UP_in_DATPs"
combinedsmall$AT_DATP_status[combinedsmall$ATstatus=="AT2" & combinedsmall$UP_in_DATPs1<median2] = "AT2_low_UP_in_DATPs"
combinedsmall$AT_DATP_status[combinedsmall$ATstatus=="AT2" & combinedsmall$UP_in_DATPs1>median2] = "AT2_high_UP_in_DATPs"

#print violin plots of CAV1 expression in AT1, ETV5 in AT2
#calculate wilcoxon p-vals of differential expression between covid and control in these groups
pdf(paste0(workingdirectory,"/Figure_3de.pdf"))
print(VlnPlot(combinedsmall[,combinedsmall$ATstatus=="AT1"],features="CAV1",pt.size=0.1,split.by="group",combine=F,split.plot=T,log=T))
print(VlnPlot(combinedsmall[,combinedsmall$ATstatus=="AT2"],features="ETV5",pt.size=0.1,split.by="group",combine=F,split.plot=T,log=T))
dev.off()
print(wilcox.test(as.numeric(combinedsmall$RNA[rownames(combinedsmall)=="CAV1",combinedsmall$ATstatus=="AT1" & combinedsmall$group=="cov"]),as.numeric(combinedsmall$RNA[rownames(combinedsmall)=="CAV1",combinedsmall$ATstatus=="AT1" & combinedsmall$group=="ctr"]))$p.val)
print(wilcox.test(as.numeric(combinedsmall$RNA[rownames(combinedsmall)=="ETV5",combinedsmall$ATstatus=="AT2" & combinedsmall$group=="cov"]),as.numeric(combinedsmall$RNA[rownames(combinedsmall)=="ETV5",combinedsmall$ATstatus=="AT2" & combinedsmall$group=="ctr"]))$p.val)

#select genes in combinedsmall that are only in the original AT1 signature, AT2 signature, or new AT1 vs 2 signature, in that order
combinedsmall1idxs = match(na.omit(covid_sigs$AT1),rownames(combinedsmall))
combinedsmall2idxs = match(na.omit(covid_sigs$AT2),rownames(combinedsmall))
combinedsmall3idxs = match(topatmarkers,rownames(combinedsmall))
covid_sigs_status = c(rep("AT1 original marker",length(combinedsmall1idxs)), rep("AT2 original marker",length(combinedsmall2idxs)))
covid_sigs_status = c(covid_sigs_status,rep("AT1/2 Cluster Markers",length(topatmarkers)))
combinedsmall = combinedsmall[c(combinedsmall1idxs,combinedsmall2idxs,combinedsmall3idxs),]
permuteval = match(rownames(combinedsmall),rownames(combinedsmall$RNA@scale.data))
combinedsmall$RNA@scale.data = combinedsmall$RNA@scale.data[permuteval,]

#score combinedsmall for cell cycle markers, assign them as cycling if G2M or S scores are above a threshold, otherwise non-cycling
combinedsmall <- CellCycleScoring(combinedsmall, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
pdf(paste0(workingdirectory,"/AT_cellcycle_scatter.pdf"))
print(FeatureScatter(combinedsmall,"G2M.Score","S.Score"))
dev.off()
combinedsmall$cycling = "Non-cycling"
combinedsmall$cycling[combinedsmall$G2M.Score>0.08 | combinedsmall$S.Score>0.08] = "Cycling"


#p-vals and median differences for UP_in_DATPs, DOWN_in_DATPs, cycling scores, and primed scores in covid vs. ctr
up_pval = wilcox.test(combinedsmall$UP_in_DATPs1[combinedsmall$group=="cov"],combinedsmall$UP_in_DATPs1[combinedsmall$group=="ctr"])$p.value
up_mediandiff = median(combinedsmall$UP_in_DATPs1[combinedsmall$group=="cov"]) - median(combinedsmall$UP_in_DATPs1[combinedsmall$group=="ctr"])
down_pval = wilcox.test(combinedsmall$DOWN_in_DATPs1[combinedsmall$group=="cov"],combinedsmall$DOWN_in_DATPs1[combinedsmall$group=="ctr"])$p.value
down_mediandiff = median(combinedsmall$DOWN_in_DATPs1[combinedsmall$group=="cov"]) - median(combinedsmall$DOWN_in_DATPs1[combinedsmall$group=="ctr"])
cycling_pval = wilcox.test(combinedsmall$cycling_AT21[combinedsmall$group=="cov"],combinedsmall$cycling_AT21[combinedsmall$group=="ctr"])$p.value
cycling_mediandiff = median(combinedsmall$cycling_AT21[combinedsmall$group=="cov"]) - median(combinedsmall$cycling_AT21[combinedsmall$group=="ctr"])
primed_pval = wilcox.test(combinedsmall$primed_AT21[combinedsmall$group=="cov"],combinedsmall$primed_AT21[combinedsmall$group=="ctr"])$p.value
primed_mediandiff = median(combinedsmall$primed_AT21[combinedsmall$group=="cov"]) - median(combinedsmall$primed_AT21[combinedsmall$group=="ctr"])

pdf(paste0(workingdirectory,"/Figure_3c.pdf"),width=14,height=14)

#make heatmap of AT1 vs. AT2 differential genes, with row and column annotations
reorder = order(combinedsmall$ATstatus,combinedsmall$group)#,combinedsmall$UP_in_DATPs1,combinedsmall$DOWN_in_DATPs1)
datamat = combinedsmall$RNA@scale.data
ATstatus = combinedsmall$ATstatus
groupvec = combinedsmall$group
primed_AT2 = combinedsmall$primed_AT21
cycling_AT2 = combinedsmall$cycling_AT21
reducedAT2_in_DATPs = combinedsmall$reducedAT2_in_DATPs1
UP_in_DATPs = combinedsmall$UP_in_DATPs1
DOWN_in_DATPs = combinedsmall$DOWN_in_DATPs1
reducedAT1_in_DATPs = combinedsmall$reducedAT1_in_DATPs1
KIarr = combined$RNA@scale.data[rownames(combined)=="MKI67",combined$ATstatus!="None"]
phasearr = combinedsmall$Phase
G2_Mscore = combinedsmall$G2M.Score
Sscore = combinedsmall$S.Score
cycling = combinedsmall$cycling

datamat = datamat[,reorder]
ATstatus = ATstatus[reorder]
groupvec = groupvec[reorder]
primed_AT2 = primed_AT2[reorder]
cycling_AT2 = cycling_AT2[reorder]
reducedAT2_in_DATPs = reducedAT2_in_DATPs[reorder]
UP_in_DATPs = UP_in_DATPs[reorder]
DOWN_in_DATPs = DOWN_in_DATPs[reorder]
reducedAT1_in_DATPs = reducedAT1_in_DATPs[reorder]
KIarr = KIarr[reorder]
phasearr = phasearr[reorder]
G2_Mscore = G2_Mscore[reorder]
Sscore = Sscore[reorder]
cycling = cycling[reorder]

annotation_col = data.frame(ATstatus = factor(ATstatus), group = factor(groupvec), cycling_status = cycling, primed_AT2_score = primed_AT2, cycling_AT2_score = cycling_AT2, UP_in_DATPs_score = UP_in_DATPs)#, DOWN_in_DATPs_score = DOWN_in_DATPs,phase = phasearr, KI67expr = KIarr,G2_Mscore = G2_Mscore, Sscore = Sscore, reducedAT1_in_DATPsscore = reducedAT1_in_DATPs, reducedAT2_in_DATPsscore = reducedAT2_in_DATPs)
rownames(annotation_col) = colnames(datamat)
annotation_row = data.frame(annotation_status = covid_sigs_status)
rownames(annotation_row) = rownames(datamat)
logdatamat = log(datamat-1.01*min(min(datamat)))
print(pheatmap(logdatamat, color = PurpleAndYellow(), annotation_col = annotation_col, annotation_row = annotation_row, show_colnames = F, cluster_rows = F, cluster_cols = F))#, legend_breaks = c(-4,-3,-2,-1,0,1,2,max(logdatamat)), legend_labels = c("-4","-3","-2","-1","0","1","2","Log Expression\n")))
dev.off()
DefaultAssay(combined) = "integrated"