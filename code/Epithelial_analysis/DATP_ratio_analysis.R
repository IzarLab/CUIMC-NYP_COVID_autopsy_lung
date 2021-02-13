#!/usr/bin/env Rscript

### title: Calculate ratio of DATP annotated cells among all AT cells, make boxplot
### of this ratio in COVID-19 and Control samples, and correlate this ratio with
### symptom-death interval author: Yiping Wang date: 02/08/2021

df_tobesummed_fine_AT = data.frame(orig.ident = data_AT_int_cov_ctr$orig.ident, group = data_AT_int_cov_ctr$group, 
    cell_type_fine = data_AT_int_cov_ctr$cell_type_fine)
df_summed_fine_AT = df_tobesummed_fine_AT %>% group_by(orig.ident, cell_type_fine, 
    group) %>% tally()
df_summed_fine_AT = df_summed_fine_AT %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

df_tobesummed_ourdatp = data.frame(orig.ident = data_AT_int_cov_ctr$orig.ident, group = data_AT_int_cov_ctr$group, 
    cell_type_ourdatp = data_AT_int_cov_ctr$cell_type_ourdatp)
df_summed_ourdatp = df_tobesummed_ourdatp %>% group_by(orig.ident, cell_type_ourdatp, 
    group) %>% tally()
df_summed_ourdatp = df_summed_ourdatp %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# for each patient sample, calculate ratio of number of cells with DATP
# annotation, versus total number of AT cells for that patient
datp_ratio_df = data.frame(datp_ratio = double, orig.ident = character(), group = character())
uniqueidents = unique(df_summed_ourdatp$orig.ident)
for (i in 1:length(uniqueidents)) {
    if (sum(df_summed_ourdatp$orig.ident == uniqueidents[i] & df_summed_ourdatp$cell_type_ourdatp == 
        "DATP") != 0) {
        tempdf = data.frame(datp_ratio = df_summed_ourdatp$n[df_summed_ourdatp$orig.ident == 
            uniqueidents[i] & df_summed_ourdatp$cell_type_ourdatp == "DATP"]/sum(df_summed_fine_AT$n[df_summed_fine_AT$orig.ident == 
            uniqueidents[i]]), orig.ident = uniqueidents[i], group = unique(df_summed_ourdatp$group[df_summed_ourdatp$orig.ident == 
            uniqueidents[i]]))
        datp_ratio_df = rbind(datp_ratio_df, tempdf)
    } else {
        tempdf = data.frame(datp_ratio = 0, orig.ident = uniqueidents[i], group = unique(df_summed_ourdatp$group[df_summed_ourdatp$orig.ident == 
            uniqueidents[i]]))
        datp_ratio_df = rbind(datp_ratio_df, tempdf)
    }
}

# plot boxplot of DATP ratios in COVID-19 and Control
datp_ratio_df$a = ""
consistentcolors = colors <- c("#006E82", "#AA0A3C", "#8214A0", "#00A0FA", "#FA5078", 
    "#005AC8", "#CC79A7", "#FAE6BE", "#0072B2", "#A0FA82", "#F0F032", "#0AB45A", 
    "#FA7850", "#14D2DC", "#FA78FA")
ggboxplot(datp_ratio_df, x = "a", y = "datp_ratio", color = "group", add = "jitter") + 
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = consistentcolors[1:2]) + 
    xlab("") + ylab("Fraction of DATPs Among AT cells")
ggsave(paste0("Figure_3k.pdf"), width = 4, height = 7)

# for each sample in datp_ratio_df, extract symptoms-death interval and enter
# into datp_ratio_df
datp_ratio_df$interval_death_symptoms_onset_days = 0
uniqueidents = unique(data_lungs_all$orig.ident)
for (i in 1:length(uniqueidents)) {
    onetime = unique(data_lungs_all$interval_death_symptoms_onset_days[data_lungs_all$orig.ident == 
        uniqueidents[i]])
    print(onetime)
    datp_ratio_df$interval_death_symptoms_onset_days[datp_ratio_df$orig.ident == 
        uniqueidents[i]] = onetime
}

datp_ratio_df = subset(datp_ratio_df, !(is.na(interval_death_symptoms_onset_days)))

# Correlate DATP ratio with symptom-death interval
pearsoncorr = (cor.test(datp_ratio_df$datp_ratio, datp_ratio_df$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(datp_ratio_df$datp_ratio, datp_ratio_df$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(datp_ratio_df, aes(x = interval_death_symptoms_onset_days, y = datp_ratio)) + 
    geom_point() + geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, 
    "= 0.0958, p-val = 0.197"))) + xlab("Days from symptom onset to death") + ylab("DATP fraction among AT cells")
ggsave("datp_ratio_vs_symptom_death_interval.pdf", width = 7, height = 7)
