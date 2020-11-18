working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"
combined = readRDS(paste0(working_directory,"/combinedwithtwoizarcontrols.rds"))

#print violin plot of DCLK1 expression in Tuft vs. Non-tuft cells
#perform wilcoxon test for DCLK1 expression in Tuft vs. Non-tuft cells
pdf(paste0(workingdirectory,"/Extended_Data_Figure_6a.pdf"))
combined$tuft_status = "Non-tuft"
combined$tuft_status[combined$overallclassification=="high_tuft_markers_short"] = "Tuft"
print(VlnPlot(combined,features="DCLK1",group.by="tuft_status",pt.size=0.1,split.by="group",combine=F,split.plot=T,log=F))
dev.off()
print(wilcox.test(as.numeric(combined$RNA[rownames(combined)=="DCLK1",combined$tuft_status=="Tuft"]),as.numeric(combined$RNA[rownames(combined)=="DCLK1",combined$tuft_status=="Non-tuft"]))$p.val)

#print violin plots of tuft-1 and tuft-2 signature strengths in tuft cells
combined = SetIdent(combined, cells=colnames(combined)[combined$overallclassification=="high_tuft_markers_short"],value="high_tuft_markers_short")
#combined = AddModuleScore(object = combined, features = list(na.omit(covid_sigs$tuft_1_montoro)), name = "tuft_1_montoro", assay = 'RNA', search = T)
#combined = AddModuleScore(object = combined, features = list(na.omit(covid_sigs$tuft_2_montoro)), name = "tuft_2_montoro", assay = 'RNA', search = T)
pdf(paste0(workingdirectory,"/Extended_Data_Figure_6bc.pdf"))
plotobj = VlnPlot(combined[,combined$tuft_status=="Tuft"],features=c("tuft_1_montoro1","tuft_2_montoro1"),pt.size=0.1,combine=T,split.by = "group",split.plot=T,log=F,same.y.lims=T) + scale_y_continuous(limits = c(-0.1,0.2))
print(plotobj)
dev.off()
#perform simple binomial test for frequency of higher tuft-1 signatures over tuft-2 signatures, compared to higher tuft-2 signatures over tuft-1 signatures
sum1 = sum(combined[,combined$tuft_status=="Tuft"]$tuft_2_montoro1>combined[,combined$tuft_status=="Tuft"]$tuft_1_montoro1)
sum2 = sum(combined[,combined$tuft_status=="Tuft"]$tuft_1_montoro1>combined[,combined$tuft_status=="Tuft"]$tuft_2_montoro1)
print(binom.test(sum1,sum1+sum2,0.5)$p.val)
