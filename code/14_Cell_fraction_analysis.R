working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library(ggplot2)

combined = readRDS(paste0(workingdirectory,"/combinedwithtwoizarcontrols.rds"))

combined$overallclassification[combined$overallclassification=="Not classified"] = "Other"

#merge Macrophages and Monocytes into one groups
theme_set(theme_bw())
dftosum = data.frame(orig.ident = combined$orig.ident, group=combined$group, overallclassification = combined$overallclassification)
dftosum$overallclassification[dftosum$overallclassification=="Macrophages"] = "Macrophages/Monocytes"
dftosum$overallclassification[dftosum$overallclassification=="Monocytes"] = "Macrophages/Monocytes"

#calculate frequencies of each cell type in each sample, store in dfsummed
#calculate mean and sd of cell type frequencies in covid and control, store in dfsummed2
dfsummed = dftosum %>% group_by(orig.ident,overallclassification,group) %>% tally()
dfsummed = dfsummed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
dfsummed2 = dfsummed %>% group_by(overallclassification,group) %>% summarise(meanfreq = mean(freq), sdfreq = sd(freq))

#plot either boxplots, or grouped barplots with error bars, of cell type frequency, in covid and control
ggboxplot(dfsummed[dfsummed$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils"),], x="overallclassification", y="freq", color="group", add="jitter") + ylim(0,0.8) + stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("cellcount_immune_two_izar_control_boxplot.png",width=8,height=7)
ggboxplot(dfsummed[!(dfsummed$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils")),], x="overallclassification", y="freq", color="group", add="jitter") + ylim(0,0.8) + stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("cellcount_nonimmune_two_izar_control_boxplot.png",width=13,height=7)
ggplot(dfsummed2,aes(fill=group, y=meanfreq, x=overallclassification)) + geom_bar(position="dodge",stat="identity") + geom_errorbar(aes(ymin=meanfreq-sdfreq, ymax=meanfreq+sdfreq), width=.2, position=position_dodge(.9)) + theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("cellcount_two_izar_control_barplot.png",width=14,height=7)

dfsummed2$immunestatus = "Non-Immune"
dfsummed2$immunestatus[(dfsummed2$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes"))] = "Immune"
library(plotly)
library(orca)

#create new dataframe dfsummed2, with overall frequency of all cell types across covid or control samples
dfsummedoverall = data.frame(celltype=rownames(prop.table(table(combined$overallclassification[combined$group=="cov"]))),frequency=prop.table(table(combined$overallclassification[combined$group=="cov"])),group="cov")
dfsummedoverall2 = data.frame(celltype=rownames(prop.table(table(combined$overallclassification[combined$group=="ctr"]))),frequency=prop.table(table(combined$overallclassification[combined$group=="ctr"])),group="ctr")
dfsummedoverall3 = rbind(dfsummedoverall,dfsummedoverall2)
dfsummedoverall3$overallclassification=dfsummedoverall3$celltype
dfsummedoverall3$meanfreq=dfsummedoverall3$frequency.Freq
dfsummedoverall3$immunestatus = "Non-Immune"
dfsummedoverall3$immunestatus[(dfsummed2$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes"))] = "Immune"
dfsummed2 = dfsummedoverall3

#plot pie charts of immune cell type proportions and non-immune proportions, in either covid samples, control samples, or both
#also plot pie charts of total immune vs. non-immune proportions across covid and control samples
titlesarr = c("Covid Immune Vs. Non-Immune","Control Immune Vs. Non-Immune","Covid Immune","Control Immune","Covid Non-Immune","Control Non-Immune","Covid","Control")
suffixesarr = c("cov_immune_vs_nonimmune","ctr_immune_vs_nonimmune","cov_immune","ctr_immune","cov_nonimmune","ctr_nonimmune","cov","ctr")
figurenamesarr = c("Extended_Data_Figure_1d","Extended_Data_Figure_1c","Figure_1f","Figure_1e","Figure_1h","Figure_1g","Extended_Data_Figure_1b","Extended_Data_Figure_1a")
for (z in 1:length(suffixesarr))
{
  if (z %in% c(1,3,5,7))
  {
    selectarr = (dfsummed2$group=="cov")
  } else {
    selectarr = (dfsummed2$group=="ctr")
  }
  if (z %in% c(3,4))
  {
    selectarr = (selectarr & (dfsummed2$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes")))
  }
  if (z %in% c(5,6))
  {
    selectarr = (selectarr & !(dfsummed2$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes")))
  }
  if (z %in% c(1,2))
  {
    p <- plot_ly(dfsummed2[selectarr,], labels = ~immunestatus, values = ~meanfreq, type = 'pie',textposition = 'outside',textinfo = 'label+percent', sort=F)# %>%
  } else {
    p <- plot_ly(dfsummed2[selectarr,], labels = ~overallclassification, values = ~meanfreq, type = 'pie',textposition = 'outside',textinfo = 'label+percent', sort=F)# %>%
  }
  p <- p %>%  layout(title = 'Cell Frequency',
	   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, type = "category", categoryorder = "array", categoryarray = unique(dfsummed2$overallclassification)),
	   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, type = "category", categoryorder = "array", categoryarray = unique(dfsummed2$overallclassification)),
	   showlegend=F)

  m <- list(
    l = 200,
    r = 200,
    b = 400,
    t = 600,
    pad = 4
  )
  p <- p %>% layout(autosize = F, width = 1500, height = 1500, margin = m)
  t <- list(
    family = "sans serif",
    size = 16,
    color = 'blue')
  p <- p %>% layout(title=paste0(titlesarr[z]," Cell Fractions"),font=t)
  orca(p, paste0(figurenamesarr[z],".pdf"),format="pdf")
}