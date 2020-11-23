working_directory = "~/CUIMC-NYP_COVID_autopsy_lung/code"

library(ggplot2)

combined = readRDS(paste0(workingdirectory,"/combinedwithtwoizarcontrols.rds"))

combined$overallclassification[combined$overallclassification=="Not classified"] = "Other"

#merge Macrophages and Monocytes into one groups
theme_set(theme_bw())
df_tobesummed = data.frame(orig.ident = combined$orig.ident, group=combined$group, overallclassification = combined$overallclassification)
df_tobesummed$overallclassification[df_tobesummed$overallclassification=="Macrophages"] = "Macrophages/Monocytes"
df_tobesummed$overallclassification[df_tobesummed$overallclassification=="Monocytes"] = "Macrophages/Monocytes"

#calculate frequencies of each cell type in each sample, store in df_summed
#calculate mean and sd of cell type frequencies in covid and control, store in dfsummed2
df_summed = df_tobesummed %>% group_by(orig.ident,overallclassification,group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
df_summed_across_samples = df_summed %>% group_by(overallclassification,group) %>% summarise(meanfreq = mean(freq), sdfreq = sd(freq))

#plot either boxplots, or grouped barplots with error bars, of cell type frequency, in covid and control
ggboxplot(df_summed[df_summed$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils"),], x="overallclassification", y="freq", color="group", add="jitter") + ylim(0,0.8) + stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(workingdirectory,"/cellcount_immune_two_izar_control_boxplot.png"),width=8,height=7)

ggboxplot(df_summed[!(df_summed$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils")),], x="overallclassification", y="freq", color="group", add="jitter") + ylim(0,0.8) + stat_compare_means(aes(group=group), label="p.signif", method="wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(workingdirectory,"/cellcount_nonimmune_two_izar_control_boxplot.png"),width=13,height=7)

ggplot(df_summed_across_samples,aes(fill=group, y=meanfreq, x=overallclassification)) + geom_bar(position="dodge",stat="identity") + geom_errorbar(aes(ymin=meanfreq-sdfreq, ymax=meanfreq+sdfreq), width=.2, position=position_dodge(.9)) + theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(workingdirectory,"cellcount_two_izar_control_barplot.png"),width=14,height=7)

df_summed_across_samples$immunestatus = "Non-Immune"
df_summed_across_samples$immunestatus[(dfsummed2$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes"))] = "Immune"
library(plotly)
library(orca)

#create new dataframe dfsummed2, with overall frequency of all cell types across covid or control samples
df_summed_overall_cov = data.frame(celltype=rownames(prop.table(table(combined$overallclassification[combined$group=="cov"]))),frequency=prop.table(table(combined$overallclassification[combined$group=="cov"])),group="cov")
df_summed_overall_ctr = data.frame(celltype=rownames(prop.table(table(combined$overallclassification[combined$group=="ctr"]))),frequency=prop.table(table(combined$overallclassification[combined$group=="ctr"])),group="ctr")
df_summed_overall = rbind(df_summed_overall_cov,df_summed_overall_ctr)
df_summed_overall$overallclassification=df_summed_overall$celltype
df_summed_overall$meanfreq=df_summed_overall3frequency.Freq
df_summed_overall$immunestatus = "Non-Immune"
df_summed_overall$immunestatus[(df_summed_overall$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes"))] = "Immune"

#plot pie charts of immune cell type proportions and non-immune proportions, in either covid samples, control samples, or both
#also plot pie charts of total immune vs. non-immune proportions across covid and control samples
titlesarr = c("Covid Immune Vs. Non-Immune","Control Immune Vs. Non-Immune","Covid Immune","Control Immune","Covid Non-Immune","Control Non-Immune","Covid","Control")
suffixesarr = c("cov_immune_vs_nonimmune","ctr_immune_vs_nonimmune","cov_immune","ctr_immune","cov_nonimmune","ctr_nonimmune","cov","ctr")
figurenamesarr = c("Extended_Data_Figure_1d","Extended_Data_Figure_1c","Figure_1f","Figure_1e","Figure_1h","Figure_1g","Extended_Data_Figure_1b","Extended_Data_Figure_1a")
for (z in 1:length(suffixesarr))
{
  if (z %in% c(1,3,5,7))
  {
    selectarr = (df_summed_overall$group=="cov")
  } else {
    selectarr = (df_summed_overall$group=="ctr")
  }
  if (z %in% c(3,4))
  {
    selectarr = (selectarr & (df_summed_overall$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes")))
  }
  if (z %in% c(5,6))
  {
    selectarr = (selectarr & !(df_summed_overall$overallclassification %in% c("B-cells","Plasma cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Macrophages/Monocytes","NK cells","Neutrophils","Macrophages","Monocytes")))
  }
  if (z %in% c(1,2))
  {
    p <- plot_ly(df_summed_overall[selectarr,], labels = ~immunestatus, values = ~meanfreq, type = 'pie',textposition = 'outside',textinfo = 'label+percent', sort=F)# %>%
  } else {
    p <- plot_ly(df_summed_overall[selectarr,], labels = ~overallclassification, values = ~meanfreq, type = 'pie',textposition = 'outside',textinfo = 'label+percent', sort=F)# %>%
  }
  p <- p %>%  layout(title = 'Cell Frequency',
	   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, type = "category", categoryorder = "array", categoryarray = unique(df_summed_overall$overallclassification)),
	   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, type = "category", categoryorder = "array", categoryarray = unique(df_summed_overall$overallclassification)),
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
  orca(p, paste0(workingdirectory,"/",figurenamesarr[z],".pdf"),format="pdf")
}