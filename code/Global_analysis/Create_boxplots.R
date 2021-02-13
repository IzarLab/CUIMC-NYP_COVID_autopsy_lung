#!/usr/bin/env Rscript

### title: Generate boxplots of cell type frequencies, grouped by either disease
### status, or sex of samples author: Yiping Wang date: 02/08/2021

consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
    "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
    "#FA7850", "#14D2DC", "#FA78FA")

# calculate frequencies of cell_type_main classes in each sample
df_tobesummed = data.frame(orig.ident = data_lungs_all$orig.ident, group = data_lungs_all$group, 
    cell_type_main = data_lungs_all$cell_type_main, cell_type_fine = data_lungs_all$cell_type_fine, 
    cell_type_intermediate = data_lungs_all$cell_type_intermediate, immune_status = data_lungs_all$immune_status)
df_summed = df_tobesummed %>% group_by(orig.ident, cell_type_intermediate, cell_type_main, 
    immune_status, group) %>% tally()
df_summed = df_summed %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# calculate frequencies of cell_type_fine classes in each sample
df_tobesummed_fine = data.frame(orig.ident = data_lungs_all$orig.ident, group = data_lungs_all$group, 
    cell_type_fine = data_lungs_all$cell_type_fine, immune_status = data_lungs_all$immune_status)
df_summed_fine = df_tobesummed_fine %>% group_by(orig.ident, cell_type_fine, immune_status, 
    group) %>% tally()
df_summed_fine = df_summed_fine %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

write.table(df_summed_fine, "boxplot_proportions_fine.csv", sep = ",", row.names = F, 
    col.names = T, quote = F)

# calculate frequencies of cell_type_main classes in each sample
df_tobesummed_main = data.frame(orig.ident = data_lungs_all$orig.ident, group = data_lungs_all$group, 
    cell_type_main = data_lungs_all$cell_type_main, immune_status = data_lungs_all$immune_status)
df_summed_main = df_tobesummed_main %>% group_by(orig.ident, cell_type_main, immune_status, 
    group) %>% tally()
df_summed_main = df_summed_main %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# calculate frequencies of cell_type_intermediate classes in each sample
df_tobesummed_intermediate = data.frame(orig.ident = data_lungs_all$orig.ident, group = data_lungs_all$group, 
    cell_type_intermediate = data_lungs_all$cell_type_intermediate, immune_status = data_lungs_all$immune_status)
df_summed_intermediate = df_tobesummed_intermediate %>% group_by(orig.ident, cell_type_intermediate, 
    immune_status, group) %>% tally()
df_summed_intermediate = df_summed_intermediate %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# calculate frequencies of cell_type_intermediate classes only within immune
# compartment
df_tobesummed_intermediate_immune = data.frame(orig.ident = data_lungs_all$orig.ident, 
    group = data_lungs_all$group, cell_type_intermediate = data_lungs_all$cell_type_intermediate, 
    immune_status = data_lungs_all$immune_status)
df_tobesummed_intermediate_immune = df_tobesummed_intermediate_immune[df_tobesummed_intermediate_immune$immune_status == 
    "Immune", ]
df_summed_intermediate_immune = df_tobesummed_intermediate_immune %>% group_by(orig.ident, 
    cell_type_intermediate, immune_status, group) %>% tally()
df_summed_intermediate_immune = df_summed_intermediate_immune %>% group_by(orig.ident) %>% 
    mutate(freq = n/sum(n))

# calculate frequencies of cell_type_intermediate classes only within nonimmune
# compartment
df_tobesummed_intermediate_nonimmune = data.frame(orig.ident = data_lungs_all$orig.ident, 
    group = data_lungs_all$group, cell_type_intermediate = data_lungs_all$cell_type_intermediate, 
    immune_status = data_lungs_all$immune_status)
df_tobesummed_intermediate_nonimmune = df_tobesummed_intermediate_nonimmune[df_tobesummed_intermediate_nonimmune$immune_status == 
    "Non-immune", ]
df_summed_intermediate_nonimmune = df_tobesummed_intermediate_nonimmune %>% group_by(orig.ident, 
    cell_type_intermediate, immune_status, group) %>% tally()
df_summed_intermediate_nonimmune = df_summed_intermediate_nonimmune %>% group_by(orig.ident) %>% 
    mutate(freq = n/sum(n))

# make boxplots of cell type frequencies in COVID-19 and Control, using cell type
# frequencies defined above
ggboxplot(df_summed_main, x = "cell_type_main", y = "freq", color = "group", add = "jitter") + 
    ylim(0, 0.8) + stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("Figure_1E.pdf"), width = 11, height = 7)
ggboxplot(subset(df_summed_main, immune_status == "Immune"), x = "cell_type_main", 
    y = "freq", color = "group", add = "jitter") + ylim(0, 0.8) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("cell_type_main_lungs_all_immune_boxplot.pdf"), width = 8, height = 7)
ggboxplot(subset(df_summed_main, immune_status == "Non-immune"), x = "cell_type_main", 
    y = "freq", color = "group", add = "jitter") + ylim(0, 0.8) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test", size = 2) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("cell_type_main_lungs_all_nonimmune_boxplot.pdf"), width = 8, height = 7)

# Extended Data Figure 2B
ggboxplot(df_summed_intermediate_nonimmune, x = "cell_type_intermediate", y = "freq", 
    color = "group", add = "jitter") + ylim(0, 1) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test", size = 2, label.y = 0.9) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("Extended_Data_Figure_2B.pdf"), width = 6, 
    height = 5)
# Extended Data Figure 2C
ggboxplot(df_summed_intermediate_immune, x = "cell_type_intermediate", y = "freq", 
    color = "group", add = "jitter") + ylim(0, 1) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test", size = 2, label.y = 0.9) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("Extended_Data_Figure_2C.pdf"), width = 6, 
    height = 5)

ggboxplot(df_summed_fine_fibroblast, x = "cell_type_fine", y = "freq", color = "group", 
    add = "jitter") + ylim(0, 0.8) + stat_compare_means(aes(group = group), label = "p.format", 
    method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("cell_type_fine_lungs_all_fibroblast_boxplot.pdf"), width = 13, height = 7)
ggboxplot(df_summed_fine_fibroblast_norm_allcells, x = "cell_type_fine", y = "freq", 
    color = "group", add = "jitter") + ylim(0, 0.3) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("cell_type_fine_lungs_all_fibroblast_norm_allcells_boxplot.pdf"), width = 13, 
    height = 7)

# Extended Data Figure 2A
ggboxplot(df_summed_intermediate, x = "cell_type_intermediate", y = "freq", color = "group", 
    add = "jitter") + ylim(0, 0.6) + stat_compare_means(aes(group = group), label = "p.format", 
    method = "wilcox.test", size = 2, label.y = 0.55) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
ggsave(paste0("Extended_Data_Figure_2A.pdf"), width = 12, 
    height = 5)

# make boxplot of all cell_type_intermediate frequencies in COVID-19 and Control, spread
# across two rows
allcelltypesintermediate = sort(unique(df_summed_intermediate$cell_type_intermediate))
df_summed_intermediate$class = 1
df_summed_intermediate$class[df_summed_intermediate$cell_type_intermediate %in% allcelltypesintermediate[11:20]] = 2
p1 = ggboxplot(df_summed_intermediate[df_summed_intermediate$class == 1, ], x = "cell_type_intermediate", 
    y = "freq", color = "group", add = "jitter") + xlab("") + ylim(0, 0.8) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test", size = 3, label.y = 0.7) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, size = 7)) + scale_colour_manual(values = consistentcolors[1:2])
p2 = ggboxplot(df_summed_intermediate[df_summed_intermediate$class == 2, ], x = "cell_type_intermediate", 
    y = "freq", color = "group", add = "jitter") + guides(color = FALSE) + ylim(0, 
    0.8) + stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test", 
    size = 3, label.y = 0.7) + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
    size = 7)) + scale_colour_manual(values = consistentcolors[1:2])

pdf("cell_type_intermediate_lungs_all_boxplot.pdf", width = 12, height = 10)
print(plot_grid(p1, p2, labels = "", nrow = 2, align = "hv", axis = "tblr"))
dev.off()

# calculate ratio of AT2 to AT1 cells in each patient, compare ratio in COVID-19 vs.
# Control samples in a boxplot
at2freq = df_summed_intermediate[df_summed_intermediate$cell_type_intermediate == 
    "AT2", ]$freq
at1freq = df_summed_intermediate[df_summed_intermediate$cell_type_intermediate == 
    "AT1", ]$freq
group = df_summed_intermediate[df_summed_intermediate$cell_type_intermediate == "AT1", 
    ]$group
for (i in 1:length(at2freq)) {
    atratio = at2freq/at1freq
}
df_at = data.frame(atratio = atratio, group = group, a = "")
ggboxplot(df_at, x = "a", y = "atratio", color = "group", add = "jitter") + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2]) + ylab("AT2/AT1 Ratio") + 
    xlab("")
ggsave(paste0("Figure_3K.pdf"), width = 4, 
    height = 7)

# calculate frequencies of fibroblast classes, merging Intermediate pathological
# FB with Pathological FB
df_tobesummed_fine_fibroblast = data.frame(orig.ident = data_lungs_all$orig.ident, 
    group = data_lungs_all$group, cell_type_fine = data_lungs_all$cell_type_fine)
df_tobesummed_fine_fibroblast$cell_type_fine[df_tobesummed_fine_fibroblast$cell_type_fine == 
    "Intermediate pathological FB"] = "Pathological FB"
df_tobesummed_fine_fibroblast = df_tobesummed_fine_fibroblast[df_tobesummed_fine_fibroblast$cell_type_fine %in% 
    c("Adventitial FB", "Alveolar FB", "Mesothelial FB", "Other FB1", "Other FB2", 
        "Other FB3", "Other FB", "Intermediate pathological FB", "Pathological FB"), 
    ]
df_summed_fine_fibroblast = df_tobesummed_fine_fibroblast %>% group_by(orig.ident, 
    cell_type_fine, group) %>% tally()
df_summed_fine_fibroblast = df_summed_fine_fibroblast %>% group_by(orig.ident) %>% 
    mutate(freq = n/sum(n))

# calculate frequencies of pathological fibroblasts in each patient, compare ratio in COVID-19 vs.
# Control samples in a boxplot
pathological_FB_freq = df_summed_intermediate[df_summed_intermediate$cell_type_intermediate == 
    "Pathological FB", ]$freq
group = df_summed_intermediate[df_summed_intermediate$cell_type_intermediate == "Pathological FB", 
    ]$group
df_pfb = data.frame(pathological_FB_freq = pathological_FB_freq, group = group, a = "")
ggboxplot(df_at, x = "a", y = "pathological_FB_freq", color = "group", add = "jitter") + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2]) + ylab("Fraction among Fibroblasts") + 
    xlab("Pathological Fibroblasts")
ggsave(paste0("Figure_4G.pdf"), width = 4, 
    height = 7)

# calculate frequencies of fibroblast classes, normalizing by total frequency of
# all cells
df_tobesummed_fine_fibroblast_norm_allcells = data.frame(orig.ident = data_lungs_all$orig.ident, 
    group = data_lungs_all$group, cell_type_fine = data_lungs_all$cell_type_fine)
df_summed_fine_fibroblast_norm_allcells = df_tobesummed_fine_fibroblast_norm_allcells %>% 
    group_by(orig.ident, cell_type_fine, group) %>% tally()
df_summed_fine_fibroblast_norm_allcells = df_summed_fine_fibroblast_norm_allcells %>% 
    group_by(orig.ident) %>% mutate(freq = n/sum(n))
df_summed_fine_fibroblast_norm_allcells = subset(df_summed_fine_fibroblast_norm_allcells, 
    cell_type_fine %in% c("Adventitial FB", "Alveolar FB", "Mesothelial FB", "Other FB1", 
        "Other FB2", "Other FB3", "Other FB", "Intermediate pathological FB", "Pathological FB"))

df_summed_fine_fibroblast_norm_allcells = subset(df_summed_fine_fibroblast_norm_allcells, 
    fibroblast_type != "")

# plot COVID-19 and Control frequencies of fibroblast classes in boxplot
ggboxplot(df_summed_fine_fibroblast_norm_allcells, x = "fibroblast_type", y = "freq", 
    color = "group", add = "jitter") + ylim(c(0, 0.4)) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2]) + ylab("Fraction of Cells")
ggsave(paste0("Extended_Data_Figure_12D.pdf"), width = 4, height = 7)

# calculate frequencies of macrophage classes, normalizing by total frequency of
# immune cells
df_tobesummed_fine_macrophage_norm_allimmunecells = data.frame(orig.ident = data_lungs_all$orig.ident, 
    group = data_lungs_all$group, cell_type_fine = data_lungs_all$cell_type_fine, 
    immune_status = data_lungs_all$immune_status)
df_tobesummed_fine_macrophage_norm_allimmunecells = subset(df_tobesummed_fine_macrophage_norm_allimmunecells, 
    immune_status == "Immune")
df_summed_fine_macrophage_norm_allimmunecells = df_tobesummed_fine_macrophage_norm_allimmunecells %>% 
    group_by(orig.ident, cell_type_fine, group) %>% tally()
df_summed_fine_macrophage_norm_allimmunecells = df_summed_fine_macrophage_norm_allimmunecells %>% 
    group_by(orig.ident) %>% mutate(freq = n/sum(n))
df_summed_fine_macrophage_norm_allimmunecells = subset(df_summed_fine_macrophage_norm_allimmunecells, 
    cell_type_fine %in% c("Alveolar macrophages", "Monocyte-derived macrophages", 
        "Monocytes", "Transitioning MDM"))

# plot COVID-19 and Control frequencies of macrophage classes in boxplot
ggboxplot(df_summed_fine_macrophage_norm_allimmunecells, x = "cell_type_fine", y = "freq", 
    color = "group", add = "jitter") + ylim(c(0, 0.4)) + stat_compare_means(aes(group = group), 
    label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2]) + ylab("Fraction of Cells")
ggsave(paste0("Extended_Data_Figure_4G.pdf"), width = 4, height = 7)

# for either just COVID-19 samples, just Control samples, or for all samples,
# perform similar analyses as above, but this time, comparing cell frequencies
# divided by sex, rather than disease status makes Extended Data Figure 3A and B
boxplotgroups = c("COVID-19", "Control", "")
for (i in 1:length(boxplotgroups)) {
    if (boxplotgroups[i] != "") {
        data_lungs_all_temp = subset(data_lungs_all, group == boxplotgroups[i])
        df_tobesummed_sex = data.frame(orig.ident = data_lungs_all_temp$orig.ident, 
            sex = data_lungs_all_temp$sex, cell_type_main = data_lungs_all_temp$cell_type_main, 
            cell_type_fine = data_lungs_all_temp$cell_type_fine, cell_type_intermediate = data_lungs_all_temp$cell_type_intermediate, 
            immune_status = data_lungs_all_temp$immune_status)
        df_tobesummed_intermediate_sex = data.frame(orig.ident = data_lungs_all_temp$orig.ident, 
            sex = data_lungs_all_temp$sex, cell_type_intermediate = data_lungs_all_temp$cell_type_intermediate, 
            immune_status = data_lungs_all_temp$immune_status)
        df_tobesummed_main_sex = data.frame(orig.ident = data_lungs_all_temp$orig.ident, 
            sex = data_lungs_all_temp$sex, cell_type_main = data_lungs_all_temp$cell_type_main, 
            immune_status = data_lungs_all_temp$immune_status)
        if (boxplotgroups[i] == "Control") {
            suffix = "_ctr"
        } else {
            suffix = "_cov"
        }
    } else {
        df_tobesummed_sex = data.frame(orig.ident = data_lungs_all$orig.ident, sex = data_lungs_all$sex, 
            cell_type_main = data_lungs_all$cell_type_main, cell_type_fine = data_lungs_all$cell_type_fine, 
            cell_type_intermediate = data_lungs_all$cell_type_intermediate, immune_status = data_lungs_all$immune_status)
        df_tobesummed_intermediate_sex = data.frame(orig.ident = data_lungs_all$orig.ident, 
            sex = data_lungs_all$sex, cell_type_intermediate = data_lungs_all$cell_type_intermediate, 
            immune_status = data_lungs_all$immune_status)
        df_tobesummed_main_sex = data.frame(orig.ident = data_lungs_all$orig.ident, 
            sex = data_lungs_all$sex, cell_type_main = data_lungs_all$cell_type_main, 
            immune_status = data_lungs_all$immune_status)
        suffix = ""
    }
    df_summed_sex = df_tobesummed_sex %>% group_by(orig.ident, cell_type_intermediate, 
        cell_type_main, immune_status, sex) %>% tally()
    df_summed_sex = df_summed_sex %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
    
    df_summed_main_sex = df_tobesummed_main_sex %>% group_by(orig.ident, cell_type_main, 
        immune_status, sex) %>% tally()
    df_summed_main_sex = df_summed_main_sex %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))
    df_summed_intermediate_sex = df_tobesummed_intermediate_sex %>% group_by(orig.ident, 
        cell_type_intermediate, immune_status, sex) %>% tally()
    df_summed_intermediate_sex = df_summed_intermediate_sex %>% group_by(orig.ident) %>% 
        mutate(freq = n/sum(n))
    
    ggboxplot(df_summed_main_sex, x = "cell_type_main", y = "freq", color = "sex", 
        add = "jitter") + ylim(0, 0.8) + stat_compare_means(aes(group = sex), label = "p.format", 
        method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_colour_manual(values = consistentcolors[1:2])
    ggsave(paste0("cell_type_main_lungs_all_sex", suffix, "_boxplot.pdf"), width = 11, 
        height = 7)
    
    # This plot can be either Extended Data Sex Differences a or b
    ggboxplot(df_summed_intermediate_sex, x = "cell_type_intermediate", y = "freq", 
        color = "sex", add = "jitter") + ylim(0, 0.6) + stat_compare_means(aes(group = sex), 
        label = "p.format", method = "wilcox.test", size = 2, label.y = 0.55) + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2])
    ggsave(paste0("cell_type_intermediate_lungs_all_sex", suffix, "_onerow_boxplot.pdf"), 
        width = 12, height = 5)
    
    allcelltypesintermediate = sort(unique(df_summed_intermediate_sex$cell_type_intermediate))
    df_summed_intermediate_sex$class = 1
    df_summed_intermediate_sex$class[df_summed_intermediate_sex$cell_type_intermediate %in% 
        allcelltypesintermediate[11:20]] = 2
    p1 = ggboxplot(df_summed_intermediate_sex[df_summed_intermediate_sex$class == 
        1, ], x = "cell_type_intermediate", y = "freq", color = "sex", add = "jitter") + 
        xlab("") + ylim(0, 0.8) + stat_compare_means(aes(group = sex), label = "p.format", 
        method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
        size = 7)) + scale_colour_manual(values = consistentcolors[1:2])
    p2 = ggboxplot(df_summed_intermediate_sex[df_summed_intermediate_sex$class == 
        2, ], x = "cell_type_intermediate", y = "freq", color = "sex", add = "jitter") + 
        guides(color = FALSE) + ylim(0, 0.8) + stat_compare_means(aes(group = sex), 
        label = "p.format", method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1, size = 7)) + scale_colour_manual(values = consistentcolors[1:2])
    pdf(paste0("cell_type_intermediate_lungs_all_sex", suffix, "_boxplot.pdf"), width = 14, 
        height = 10)
    print(plot_grid(p1, p2, labels = "", nrow = 2, align = "hv", axis = "tblr"))
    dev.off()
}
