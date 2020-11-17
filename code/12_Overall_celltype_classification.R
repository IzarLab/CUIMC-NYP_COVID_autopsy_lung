working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

#realclashlist is a list of boolean vectors, one vector for each cell type
#each vector has one boolean value for every cell in the dataset
#TRUE values in the ith boolean vector indicates that classification as cell type i has some sort of significant overlap with another cell type at a cell
#significant clashes are determined if either another cell type has higher median scores than cell type i, or if cell type i has higher scores but the difference is not significant according to a wilcoxon p-value
mediandiffmatrix = matrix(0,nrow=length(columnscoreslist),ncol=length(columnscoreslist),dimnames=list(colnames(covid_sigs),colnames(covid_sigs)))
clashpercentmatrix = matrix(0,nrow=length(columnscoreslist),ncol=length(columnscoreslist),dimnames=list(colnames(covid_sigs),colnames(covid_sigs)))
realclashlist = list()
for (i in 1:length(columnscoreslist))
{
  print(paste0(colnames(covid_sigs)[i]," ",i))
  realclash = rep(F,length(columnscoreslist[[i]]))
  for (j in 1:length(columnscoreslist))
  {
    if (i!=j && colnames(covid_sigs)[j]!="tuft_1_montoro" && colnames(covid_sigs)[j]!="tuft_2_montoro")
    {
      totali = sum(columnscoreslist[[i]]>0)
      totaloverlap = sum(columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0)
      mediandiff = median(columnscoreslist[[i]][columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0]) - median(columnscoreslist[[j]][columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0])
      if (totaloverlap>1 || (totaloverlap>0 && mediandiff!=0))
      {
	p_val = wilcox.test(columnscoreslist[[i]][columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0],columnscoreslist[[j]][columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0])$p.value
      }
      else
      {
	p_val = 0
      }
      mediandiffmatrix[i,j] = mediandiff
      clashpercentmatrix[i,j] = totaloverlap/totali
      if (p_val>.05 || (totaloverlap>0 && mediandiff<0))
      {
	realclash[columnscoreslist[[i]]>0 & columnscoreslist[[j]]>0] = T
      }
    }
  }
  realclashlist[[i]] = realclash
}

#for cells which have FALSE values for a given vector in realclashlist, give them the corresponding cell type classification in overallclassification
#an exception is for tuft_1_montoro and tuft_2_montoro. Any cell that is originally classified as tuft_cell, and with positive montoro signatures as well, as classified as tuft_1 or 2_montoro
for(c in 1:ncol(covid_sigs))
{
  if (colnames(covid_sigs)[c]!="tuft_1_montoro" && colnames(covid_sigs)[c]!="tuft_2_montoro")
  {
    combined$overallclassification[columnscoreslist[[c]]>0 & combined$overallclassification=="Other" & !realclashlist[[c]]] = colnames(covid_sigs)[c]
  }
}
for(c in 1:ncol(covid_sigs))
{
  if (colnames(covid_sigs)[c]=="tuft_1_montoro" || colnames(covid_sigs)[c]=="tuft_2_montoro")
  {
    combined$overallclassification[columnscoreslist[[c]]>0 & combined$overallclassification=="tuft_cell"] = colnames(covid_sigs)[c]
  }
}