#!/usr/bin/env Rscript

### title: Print out boxplots of cell type frequency within myeloid compartment,
### for covid and ctr samples author: Yiping Wang date: 02/08/2021

consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
    "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
    "#FA7850", "#14D2DC", "#FA78FA")

df_tobesummed = data.frame(orig.ident = data_myeloid_cov_ctr$orig.ident, group = data_myeloid_cov_ctr$group, 
    AM_specific = data_myeloid_cov_ctr$AM_specific)
df_summed = df_tobesummed %>% group_by(orig.ident, AM_specific, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

ggboxplot(df_summed, x = "AM_specific", y = "freq", color = "group", add = "jitter") + 
    ylim(0, 0.8) + stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("Myeloid_cell_type_boxplot.pdf"), width = 11, height = 7)
