working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library(ggplot2)

combined = readRDS(paste0(workingdirectory,"/combinedwithtwoizarcontrols.rds"))
uniqueIDs = unique(combined$ID)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#read in data into IGtable variable, filtering out rows with NAs for both heavy and light chains
IGtable = read.table("data/lungs_all/b_cells/Extended_Data_Table_4.csv",header=T,sep=",",quote="")
IGtable = IGtable[!is.na(IGtable$light) & !is.na(IGtable$heavy),]
#make list of heavy and light chains that matches axes of Figure
lightorder = rev(c("IGKV1D-33","IGKV1-27","IGKV1-9","IGKV2-28","IGLVI-42","IGLVI-70","IGLV3-10","IGLV1-41","IGKV1D-13","IGKV2D-30","IGKV1-16","IGKV2D-29","IGLV3-16","IGKV1-17","IGLV9-49","IGKV1D-12","IGKV1D-39","IGLV3-9","IGLV7-43","IGKV1-6","IGLV1-44","IGLV2-11","IGLV2-8","IGLV1-40","IGLV2-23","IGLV1-51","IGLV7-46","IGLV2-14","IGLV1-47","IGKV3-11","IGLV3-1","IGLV4-69","IGKV3-15","IGKV2-30","IGKV1-39","IGLV3-21","IGKV4-1","IGKV3-20","IGLV3-19","IGLV10-54","IGLV3-25","IGLV3-27","IGKV1-5"))
heavyorder = c("IGHV3-15","IGHV3-7","IGHV4-34","IGHV3-23","IGHV1-69D","IGHV3-30","IGHV6-1","IGHV3-21","IGHV1-46","IGHV1-18","IGHV1-24","IGHV3-33","IGHV3-48","IGHV5-51","IGHV7-4-1","IGHV3-74","IGHV3-11","IGHV2-26","IGHV4-28","IGHV4-4","IGHV1-3","IGHV2-5","IGHV3-43","IGHV4-61","IGHV3-13","IGHV5-10-1","IGHV4-59","IGHV3-66","IGHV3-49","IGHV3-20","IGHV3-29","IGHV2-70","IGHV1-2","IGHV4-31","IGHV2-70D","IGHV4-39","IGHV1-69")
uniquelight = lightorder
uniqueheavy = heavyorder
#store indices of where IGtable light and heavy chains match uniquelight and uniqueheavy
IGtable$lightidxs = match(IGtable$light,uniquelight)
IGtable$heavyidxs = match(IGtable$heavy,uniqueheavy)

#make list of all light and heavy chain combinations that occur in IGtable, store in combs
#make list of combinations that are occur more than once, store in dupcombs
combs = c()
dupcombs = c()
for (i in 1:length(IGtable$light))
{
  if (sum(combs==paste0(IGtable$light[i]," ",IGtable$heavy[i]))!=0 & sum(dupcombs==paste0(IGtable$light[i]," ",IGtable$heavy[i]))==0)
  {
    dupcombs = append(dupcombs, paste0(IGtable$light[i]," ",IGtable$heavy[i]))
  }
  combs = append(combs, paste0(IGtable$light[i]," ",IGtable$heavy[i]))
}

#for duplicated combinations, it is possible that each entry may have different values in the constant field
#find all unique values in constant field for each duplicated combination
#in the cell for each duplicated combination, subdivide the cell along the light chain axis, so that each subcell will represent a unique constant value
#store the endpoints of each subcell in lightidxs and lightidxsceil fields
duptable = data.frame(light=character(),heavy=character(),constant=character(),lightidxs=integer(),lightidxsceil=integer(),heavyidxs=integer())
for (i in 1:length(dupcombs))
{
  words = strsplit(dupcombs[i]," ")[[1]]
  smalltable = na.omit(IGtable[IGtable$light==words[1] & IGtable$heavy==words[2],])
  smalltableconstant = unique(smalltable$constant)
  if (length(smalltableconstant)>1)
  {
    for (j in 1:length(smalltableconstant))
    {
      duptable = rbind(duptable,data.frame(light=words[1],heavy=words[2],constant=smalltableconstant[j],lightidxs=match(words[1],uniquelight)+(j-1)/length(smalltableconstant),lightidxsceil=match(words[1],uniquelight)+j/length(smalltableconstant),heavyidxs=match(words[2],uniqueheavy)))
    }
  }
}

#create heatmap of all light ahd heavy chain combinations, with subcells for multiple constant values appearing in one combination
ggplot(IGtable, aes(x=heavyidxs+0.5,y=lightidxs+0.5,fill=constant,width=1,height=1))+geom_tile()+geom_rect(data=duptable, aes(xmin=heavyidxs,xmax=heavyidxs+1,ymin=lightidxs,ymax=lightidxsceil,fill=constant), colour="black", size=0.3)+scale_x_continuous("Heavy Chain",breaks=1:(length(uniqueheavy)+1),labels=append(uniqueheavy,""))+scale_y_continuous("Light Chain",breaks=1:(length(uniquelight)+1),labels=append(uniquelight,"")) + scale_fill_manual(breaks = c("IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"),values=cbPalette[1:6]) + theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.grid = element_line(size = .3, colour="black"), panel.ontop = T)
ggsave("Figure_3f.pdf")

#count number of times each IG chain combination occurs in IGtable, store in sharedIGtable, along with lightidxs and heavyidxs
sharedIGtable = data.frame(heavy=character(),light=character(),occurrences=integer(),heavyidxs=integer(),lightidxs=integer())
for (i in 1:length(uniqueIDs))
{
  IGtablesmall = IGtable[IGtable$barcode %in% names(combined$ID[combined$ID==uniqueIDs[i]]),]
  if (dim(IGtablesmall)[1]!=0)
  {
    checktable = data.frame(heavy=character(),light=character())
    for (j in 1:dim(IGtablesmall)[1])
    {
      if (!is.na(IGtablesmall$heavy[j]) & !is.na(IGtablesmall$light[j]))
      {
        if (sum(checktable$heavy==IGtablesmall$heavy[j] & checktable$light==IGtablesmall$light[j])==0)
	{
	  smalltable = data.frame(heavy=IGtablesmall$heavy[j],light=IGtablesmall$light[j])
          checktable = rbind(checktable,smalltable)

	  matchingarr = (sharedIGtable$heavy==IGtablesmall$heavy[j] & sharedIGtable$light==IGtablesmall$light[j])
	  if (sum(matchingarr)==0)
	  {
	    smalltable = data.frame(heavy=IGtablesmall$heavy[j],light=IGtablesmall$light[j],occurrences=1,heavyidxs=match(IGtablesmall$heavy[j],uniqueheavy),lightidxs=match(IGtablesmall$light[j],uniquelight))
	    sharedIGtable = rbind(sharedIGtable,smalltable)
	  }
	  else {
	    sharedIGtable$occurrences[matchingarr] = sharedIGtable$occurrences[matchingarr]+1
	  }
	}
      }
    }
  }
}

#print heatmap of number of times each IG chain combination occurs in IGtable
sharedIGtable$occurrences = as.character(sharedIGtable$occurrences)
ggplot(sharedIGtable, aes(x=heavyidxs+0.5,y=lightidxs+0.5,fill=occurrences,width=1,height=1))+geom_tile()+scale_x_continuous("Heavy Chain",breaks=1:(length(uniqueheavy)+1),labels=append(uniqueheavy,""))+scale_y_continuous("Light Chain",breaks=1:(length(uniquelight)+1),labels=append(uniquelight,"")) + scale_fill_manual(breaks = c("1","2","3","4"),values=cbPalette[1:4]) + theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.grid = element_line(size = .3, colour="black"), panel.ontop = T)
ggsave("IGheatmap_shared_across_samples.pdf")

#create heatmaps of IG combinations that occur in each individual sample
for (i in 1:length(uniqueIDs))
{
  if (uniqueIDs[i] %in% c("14","15"))
  {
    suffix = "ctr"
  }
  else
  {
    suffix = "cov"
  }
  IGtablesmall = IGtable[IGtable$barcode %in% names(combined$ID[combined$ID==uniqueIDs[i]]),]

  if (dim(IGtablesmall)[1]!=0)
  {
    combssmall = c()
    dupcombssmall = c()
    duptablesmall = data.frame(light=character(),heavy=character(),constant=character(),lightidxs=integer(),lightidxsceil=integer(),heavyidxs=integer())
    for (j in 1:length(IGtablesmall$light))
    {
      if (sum(combssmall==paste0(IGtablesmall$light[j]," ",IGtablesmall$heavy[j]))!=0 & sum(dupcombssmall==paste0(IGtablesmall$light[j]," ",IGtablesmall$heavy[j]))==0)
      {
	dupcombssmall = append(dupcombssmall, paste0(IGtablesmall$light[j]," ",IGtablesmall$heavy[j]))
      }
      combssmall = append(combssmall, paste0(IGtablesmall$light[j]," ",IGtablesmall$heavy[j]))
    }
    if (length(dupcombssmall)!=0)
    {
      for (k in 1:length(dupcombssmall))
      {
	words = strsplit(dupcombssmall[k]," ")[[1]]
	smalltable = na.omit(IGtablesmall[IGtablesmall$light==words[1] & IGtablesmall$heavy==words[2],])
	smalltableconstant = unique(smalltable$constant)
	if (length(smalltableconstant)>1)
	{
	  for (l in 1:length(smalltableconstant))
	  {
	    duptablesmall = rbind(duptablesmall,data.frame(light=words[1],heavy=words[2],constant=smalltableconstant[l],lightidxs=match(words[1],uniquelight)+(l-1)/length(smalltableconstant),lightidxsceil=match(words[1],uniquelight)+l/length(smalltableconstant),heavyidxs=match(words[2],uniqueheavy)))
	  }
	}
      }
    }

    p = ggplot(IGtablesmall, aes(x=heavyidxs+0.5,y=lightidxs+0.5,fill=constant,width=1,height=1))+geom_tile()
    if (dim(duptablesmall)[1]>0)
    {
      p = p + geom_rect(data=duptablesmall, aes(xmin=heavyidxs,xmax=heavyidxs+1,ymin=lightidxs,ymax=lightidxsceil,fill=constant), colour="black", size=0.3)
    }
    p = p+ scale_x_continuous("Heavy Chain",breaks=1:(length(uniqueheavy)+1),labels=append(uniqueheavy,""),limits=c(1,length(uniqueheavy)))+
    scale_y_continuous("Light Chain",breaks=1:(length(uniquelight)+1),labels=append(uniquelight,""),limits=c(1,length(uniquelight)))+
    scale_fill_manual(breaks = c("IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"),values=cbPalette[1:6]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.grid = element_line(size = .3, colour="black"), panel.ontop = T)# + xlim(0,100) + ylim(0,100)
    ggsave(paste0("IGheatmap_",uniqueIDs[i],"_",suffix,".pdf"))
  }
}