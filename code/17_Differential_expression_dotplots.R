working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library("textTinyR")

combined = readRDS(paste0(working_directory,"/combinedwithtwoizarcontrols.rds"))
DefaultAssay(combined) = "RNA"

tissuestotest = sort(unique(combined$overallclassification))
genestotest = c("ACE2","TMPRSS2","CTSL","BSG","IL18","IL1A","IL25","IFNA1A","IFNA1","IFNA2","IFNA2A","IFNA2B","IL29","IL28A","IL28B","IL1B")

#for each gene and tissue in genes and tissuestotest, calculate it's wilcoxon p-value for differential expression between covid and ctr samples, and store in filldf
#if p-val<0.1, store above information as well as log fold change in filldf2
#if p-val<.01, store above information in filldf3
filldf = matrix(nrow=length(genestotest),ncol=length(tissuestotest))
filldf2 = data.frame(tissue=character(), gene=character(), pval=double(), logfoldchange=double())
filldf3 = data.frame(tissue=character(), gene=character())
rownames(filldf) = genestotest
colnames(filldf) = tissuestotest
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

        filldf[j,i] = wilc_pval
	if (wilc_pval<.1)
	{
	  tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), gene=str_replace_all(genestotest[j],"\\."," "), pval=wilc_pval, logfoldchange=log(mean(arr1)/mean(arr2),2))
	  filldf2 = rbind(filldf2,tempdf)
	  if (wilc_pval<0.01)
	  {
	    tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), gene=str_replace_all(genestotest[j],"\\."," "))
	    filldf3 = rbind(filldf3,tempdf)
	  }
	}
      }
    }
  }
}
#write csv file of wilcoxon p-values
write.table(filldf,"specific_tissues_test_pval.csv",row.names=T,col.names=T,sep=",",quote=F)

theme_set(theme_bw())
#make dotplots of differential expression for four sets of genes, with size standing for wilcoxon p-value, and color for log fold change
filldf2_1 = filldf2[filldf2$gene %in% c("TMPRSS2","CTSL","BSG"),]
filldf2_2 = filldf2[filldf2$gene %in% c("IL18","IL1A"),]
filldf2_3 = filldf2[filldf2$gene %in% c("IL18"),]
filldf2_4 = filldf2[filldf2$gene %in% c("IL18","IL1B"),]
filldf3_1 = filldf3[filldf3$gene %in% c("TMPRSS2","CTSL","BSG"),]
filldf3_2 = filldf3[filldf3$gene %in% c("IL18","IL1A"),]
filldf3_3 = filldf3[filldf3$gene %in% c("IL18"),]
filldf3_4 = filldf3[filldf3$gene %in% c("IL18","IL1B"),]

#alternative color scale for greater clarity: scale_colour_gradient2(low="red", high="blue", mid="black", midpoint = 0, name=waiver())
#to add asterisks next to gene-tissue combinations with p-val < .01, add this code: geom_point(data = filldf3_1, color="black", size = 2, fill=NA, shape=8, position = position_nudge(x = 0.2, y = 0.2))

ggplot(filldf2_1, aes(x = tissue, y = gene, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(filldf2_1$logfoldchange),0,max(filldf2_1$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005","0.25","0.5",">=.075")) + ylab("Gene") + xlab("Cell Type")
ggsave("extra_dotplot_1.pdf",width=11,height=6)

ggplot(filldf2_2, aes(x = tissue, y = gene, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(filldf2_2$logfoldchange),0,max(filldf2_2$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005","0.25","0.5",">=.075")) + ylab("Gene") + xlab("Cell Type")
ggsave("Figure_3h.pdf",width=11,height=4)

ggplot(filldf2_3, aes(x = tissue, y = gene, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(filldf2_3$logfoldchange),0,max(filldf2_3$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005","0.25","0.5",">=.075")) + ylab("Gene") + xlab("Cell Type")
ggsave("extra_dotplot_3.pdf",width=11,height=4)

ggplot(filldf2_4, aes(x = tissue, y = gene, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(filldf2_4$logfoldchange),0,max(filldf2_4$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005","0.25","0.5",">=.075")) + ylab("Gene") + xlab("Cell Type")
ggsave("extra_dotplot_4.pdf",width=11,height=4)

# write out csv file with log fold change and p-val for IL18 and IL1A
write.table(filldf2_2,"IL18_IL1A_table.csv",row.names=T,col.names=T,sep=",",quote=F)

#perform similar analyses and create similar figures as above, for three interesting pathway signatures
misc_sigs<-read.csv('misc_signatures.csv',na.strings = '')
for(c in 1:ncol(misc_sigs)){
  misc_sigs[,c] = toupper(misc_sigs[,c])
  combined = AddModuleScore(object = combined, features = list(na.omit(misc_sigs[,c])), name = colnames(misc_sigs)[c], assay = 'RNA', search = T)
}
misc_sigs_arr = c("Type.I.interferon.abbreviated1","inflammasome.receptors1","chemotaxis1")

filldf2 = data.frame(tissue=character(), sig=character(), pval=double(), logfoldchange=double())
filldf3 = data.frame(tissue=character(), sig=character())
filldf = matrix(nrow=length(misc_sigs_arr),ncol=length(tissuestotest))
rownames(filldf) = misc_sigs_arr
colnames(filldf) = tissuestotest
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
	filldf[j,i] = wilc_pval
	if (wilc_pval<.10)
	{
	  tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), sig=str_replace_all(str_replace_all(misc_sigs_arr[j],"\\."," "),"1",""), pval=wilc_pval, logfoldchange=log(mean(arr1)/mean(arr2),2))
	  filldf2 = rbind(filldf2,tempdf)
	  if (wilc_pval<0.01)
	  {
	    tempdf = data.frame(tissue=str_replace_all(tissuestotest[i],"_"," "), sig=str_replace_all(str_replace_all(misc_sigs_arr[j],"\\."," "),"1$",""))
	    filldf3 = rbind(filldf3,tempdf)
	  }
	}
      }
    }
  }
}
write.table(filldf,"specific_tissues_test_sigs_pval.csv",row.names=T,col.names=T,sep=",",quote=F)

filldf2$sig[filldf2$sig=="Type I interferon abbreviated"] = "Type I interferon abbreviated"
filldf3$sig[filldf3$sig=="Type I interferon abbreviated"] = "Type I interferon abbreviated"
theme_set(theme_bw())

ggplot(filldf2, aes(x = tissue, y = sig, color = logfoldchange, size = pval)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", "purple", "red"), breaks = c(min(filldf2$logfoldchange),0,max(filldf2$logfoldchange))) + scale_size(trans="reverse", breaks=c(.005,.025,.05,.075),labels=c("<=.005",".025",".05",">=.075")) + ylab("Pathway") + xlab("Cell Type")
ggsave("Figure_3i.pdf",width=11,height=5)
DefaultAssay(combined) = "integrated"

write.table(filldf2,"misc_signatures_table.csv",row.names=T,col.names=T,sep=",",quote=F)
