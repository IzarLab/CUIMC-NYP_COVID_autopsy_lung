working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library("textTinyR")

combined = readRDS(paste0(working_directory,"/combinedwithtwoizarcontrols.rds"))
DefaultAssay(combined) = "RNA"

tissuestotest = sort(unique(combined$overallclassification))
genestotest = c("ACE2","TMPRSS2","CTSL","BSG","IL18","IL1A","IL25","IFNA1A","IFNA1","IFNA2","IFNA2A","IFNA2B","IL29","IL28A","IL28B","IL1B")

#for each gene and tissue in genes and tissuestotest, calculate it's wilcoxon p-value for differential expression between covid and ctr samples, and store in dotplotdf
#if p-val<0.1, store above information as well as log fold change in dotplotdf2
#if p-val<.01, store above information in dotplotdf3
dotplotdf = matrix(nrow=length(genestotest),ncol=length(tissuestotest))
dotplotdf2 = data.frame(tissue=character(), gene=character(), pval=double(), logfoldchange=double())
dotplotdf3 = data.frame(tissue=character(), gene=character())
rownames(dotplotdf) = genestotest
colnames(dotplotdf) = tissuestotest
datamat = combined$RNA@scale.data
for (i in 1:length(tissuestotest))
{
  for (j in 1:length(genestotest))
  {
    datacol1 = combined$RNA@data[rownames(combined$RNA)==genestotest[j],combined$overallclassification==tissuestotest[i] & combined$group=="cov"]
    datacol2 = combined$RNA@data[rownames(combined$RNA)==genestotest[j],combined$overallclassification==tissuestotest[i] & combined$group=="ctr"]
    if (length(datacol1)>0 && length(datacol2)>0)
    {
      arr1 = datacol1
      arr2 = datacol2
      if (length(arr1)>2 && length(arr2)>2 && (min(arr1)!=max(arr1) || min(arr2)!=max(arr2)) && sum(arr1)!=0 && sum(arr2)!=0)
      {
	print(length(arr1))
	print(length(arr2))
	arr1save = arr1
	arr2save = arr2
	wilc_pval = wilcox.test(arr1, arr2)$p.value

        dotplotdf[j,i] = wilc_pval
	if (wilc_pval<.1)
	{
	  tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), gene=str_replace_all(genestotest[j],"\\."," "), pval=wilc_pval, logfoldchange=log(mean(arr1)/mean(arr2),2))
	  dotplotdf2 = rbind(dotplotdf2,tempdf)
	  if (wilc_pval<0.01)
	  {
	    tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), gene=str_replace_all(genestotest[j],"\\."," "))
	    dotplotdf3 = rbind(dotplotdf3,tempdf)
	  }
	}
      }
    }
  }
}
#write csv file of wilcoxon p-values
write.table(dotplotdf,"specific_tissues_test_pval.csv",row.names=T,col.names=T,sep=",",quote=F)

theme_set(theme_bw())
#make dotplots of differential expression for four sets of genes, with size standing for wilcoxon p-value, and color for log fold change
dotplotdf2_2 = dotplotdf2[dotplotdf2$gene %in% c("IL18","IL1A"),]
dotplotdf3_2 = dotplotdf3[dotplotdf3$gene %in% c("IL18","IL1A"),]

#alternative color scale for greater clarity: scale_colour_gradient2(low="red", high="blue", mid="black", midpoint = 0, name=waiver())
#to add asterisks next to gene-tissue combinations with p-val < .01, add this code: geom_point(data = dotplotdf3_1, color="black", size = 2, fill=NA, shape=8, position = position_nudge(x = 0.2, y = 0.2))

ggplot(dotplotdf2_2, aes(x = tissue, y = gene, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(dotplotdf2_2$logfoldchange),0,max(dotplotdf2_2$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005","0.25","0.5",">=.075")) + ylab("Gene") + xlab("Cell Type")
ggsave(paste0(workingdirectory,"/Figure_3h.pdf"),width=11,height=4)

# write out csv file with log fold change and p-val for IL18 and IL1A
write.table(dotplotdf2_2,paste0(workingdirectory,"/IL18_IL1A_table.csv"),row.names=T,col.names=T,sep=",",quote=F)

#perform similar analyses and create similar figures as above, for three interesting pathway signatures
misc_sigs<-read.csv('misc_signatures.csv',na.strings = '')
for(c in 1:ncol(misc_sigs)){
  misc_sigs[,c] = toupper(misc_sigs[,c])
  combined = AddModuleScore(object = combined, features = list(na.omit(misc_sigs[,c])), name = colnames(misc_sigs)[c], assay = 'RNA', search = T)
}
misc_sigs_arr = c("Type.I.interferon.abbreviated1","inflammasome.receptors1","chemotaxis1")

dotplotdf2 = data.frame(tissue=character(), sig=character(), pval=double(), logfoldchange=double())
dotplotdf3 = data.frame(tissue=character(), sig=character())
dotplotdf = matrix(nrow=length(misc_sigs_arr),ncol=length(tissuestotest))
rownames(dotplotdf) = misc_sigs_arr
colnames(dotplotdf) = tissuestotest
for (i in 1:length(tissuestotest))
{
  for (j in 1:length(misc_sigs_arr))
  {
    misc_sigs_genes = na.omit(misc_sigs[,j])
    datacol = eval(parse(text=paste0("combined$",misc_sigs_arr[j])))
    datacol1 = combined$RNA@data[rownames(combined$RNA) %in% misc_sigs_genes,combined$overallclassification==tissuestotest[i] & combined$group=="cov"]
    datacol2 = combined$RNA@data[rownames(combined$RNA) %in% misc_sigs_genes,combined$overallclassification==tissuestotest[i] & combined$group=="ctr"]
    if (!is.null(dim(datacol1)) && !is.null(dim(datacol2)))
    {
      arr1 = sparse_Sums(datacol1,rowSums=F)
      arr2 = sparse_Sums(datacol2,rowSums=F)
      if (length(arr1)>2 && length(arr2)>2 && (min(arr1)!=max(arr1) || min(arr2)!=max(arr2)) && sum(arr1)!=0 && sum(arr2)!=0)
      {
	arr1save = arr1
	arr2save = arr2
	wilc_pval = wilcox.test(arr1, arr2)$p.value
	print(tissuestotest[i])
	print(misc_sigs_arr[j])
	print(wilc_pval)
	dotplotdf[j,i] = wilc_pval
	if (wilc_pval<.10)
	{
	  tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), sig=str_replace_all(str_replace_all(misc_sigs_arr[j],"\\."," "),"1",""), pval=wilc_pval, logfoldchange=log(mean(arr1)/mean(arr2),2))
	  dotplotdf2 = rbind(dotplotdf2,tempdf)
	  if (wilc_pval<0.01)
	  {
	    tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), sig=str_replace_all(str_replace_all(misc_sigs_arr[j],"\\."," "),"1$",""))
	    dotplotdf3 = rbind(dotplotdf3,tempdf)
	  }
	}
      }
    }
  }
}
write.table(dotplotdf,"specific_tissues_test_sigs_pval.csv",row.names=T,col.names=T,sep=",",quote=F)

dotplotdf2$sig[dotplotdf2$sig=="Type I interferon abbreviated"] = "Type I interferon abbreviated"
dotplotdf3$sig[dotplotdf3$sig=="Type I interferon abbreviated"] = "Type I interferon abbreviated"
theme_set(theme_bw())

ggplot(dotplotdf2, aes(x = tissue, y = sig, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(dotplotdf2$logfoldchange),0,max(dotplotdf2$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005",".025",".05",">=.075")) + ylab("Pathway") + xlab("Cell Type")
ggsave(paste0(workingdirectory,"/Figure_3i.pdf"),width=11,height=5)
DefaultAssay(combined) = "integrated"

write.table(dotplotdf2,"misc_signatures_table.csv",row.names=T,col.names=T,sep=",",quote=F)
