#!/usr/bin/env Rscript

### title: Correlate AT and Pathological FB cell types frequencies with
### symptom-death intervals author: Yiping Wang date: 02/08/2021

# obtain frequencyes of cell_type_intermediate in all patients
df_tobesummed_intermediate = data.frame(orig.ident = data_lungs_all$orig.ident, cell_type_intermediate = data_lungs_all$cell_type_intermediate)
df_summed_intermediate = df_tobesummed_intermediate %>% group_by(orig.ident, cell_type_intermediate) %>% 
    tally()
df_summed_intermediate = df_summed_intermediate %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

# extract symptom-death interval from data_lungs_all and match with
# cell_type_intermediate frequencies, and also calculate ratio of AT2 cells to
# AT1 cells
df_summed_intermediate$interval_death_symptoms_onset_days = 0
df_summed_intermediate$atratio = 0
uniqueidents = unique(data_lungs_all$orig.ident)
for (i in 1:length(uniqueidents)) {
    onetime = unique(data_lungs_all$interval_death_symptoms_onset_days[data_lungs_all$orig.ident == 
        uniqueidents[i]])
    print(onetime)
    df_summed_intermediate$interval_death_symptoms_onset_days[df_summed_intermediate$orig.ident == 
        uniqueidents[i]] = onetime
    df_summed_intermediate$atratio[df_summed_intermediate$orig.ident == uniqueidents[i]] = df_summed_intermediate$freq[df_summed_intermediate$orig.ident == 
        uniqueidents[i] & df_summed_intermediate$cell_type_intermediate == "AT2"]/df_summed_intermediate$freq[df_summed_intermediate$orig.ident == 
        uniqueidents[i] & df_summed_intermediate$cell_type_intermediate == "AT1"]
}

df_summed_intermediate = subset(df_summed_intermediate, !(is.na(interval_death_symptoms_onset_days)))

# correlate AT1 frequency with symptom-death interval
scatterplotdf = subset(df_summed_intermediate, cell_type_intermediate == "AT1")
pearsoncorr = (cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = freq)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.00633, p-val = 0.746"))) + 
    xlab("Days from symptom onset to death") + ylab("AT1 frequency")
ggsave("AT1_freq_vs_symptom_death_interval.pdf", width = 7, height = 7)

# correlate AT2 frequency with symptom-death interval
scatterplotdf = subset(df_summed_intermediate, cell_type_intermediate == "AT2")
pearsoncorr = (cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = freq)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.120, p-val = 0.146"))) + 
    xlab("Days from symptom onset to death") + ylab("AT2 frequency")
ggsave("AT2_freq_vs_symptom_death_interval.pdf", width = 7, height = 7)

# correlate AT2/AT1 ratio with symptom-death interval
scatterplotdf = subset(df_summed_intermediate, cell_type_intermediate == "AT1")
pearsoncorr = (cor.test(scatterplotdf$atratio, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$atratio, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = atratio)) + 
    geom_point() + geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, 
    "= 0.421, p-val = 0.00266"))) + xlab("Days from symptom onset to death") + ylab("AT2/1 ratio")
ggsave("Extended Data Figure Epithelial Cells r.pdf", width = 7, height = 7)

# perform same analyses as above, but with cell frequencies normalized only
# within the non-immune compartment
df_tobesummed_intermediate = data.frame(orig.ident = data_lungs_all$orig.ident, cell_type_intermediate = data_lungs_all$cell_type_intermediate, 
    immune_status = data_lungs_all$immune_status)
df_tobesummed_intermediate = subset(df_tobesummed_intermediate, immune_status == 
    "Non-immune")
df_tobesummed_intermediate$immune_status = NULL
df_summed_intermediate = df_tobesummed_intermediate %>% group_by(orig.ident, cell_type_intermediate) %>% 
    tally()
df_summed_intermediate = df_summed_intermediate %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

df_summed_intermediate$interval_death_symptoms_onset_days = 0
uniqueidents = unique(data_lungs_all$orig.ident)
for (i in 1:length(uniqueidents)) {
    onetime = unique(data_lungs_all$interval_death_symptoms_onset_days[data_lungs_all$orig.ident == 
        uniqueidents[i]])
    print(onetime)
    df_summed_intermediate$interval_death_symptoms_onset_days[df_summed_intermediate$orig.ident == 
        uniqueidents[i]] = onetime
}

df_summed_intermediate = subset(df_summed_intermediate, !(is.na(interval_death_symptoms_onset_days)))

scatterplotdf = subset(df_summed_intermediate, cell_type_intermediate == "AT1")
pearsoncorr = (cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = freq)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.00628, p-val = 0.747"))) + 
    xlab("Days from symptom onset to death") + ylab("AT1 frequency")
ggsave("AT1_freq_vs_symptom_death_interval_nonimmune_norm.pdf", width = 7, height = 7)

scatterplotdf = subset(df_summed_intermediate, cell_type_intermediate == "AT2")
pearsoncorr = (cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = freq)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.207, p-val = 0.0502"))) + 
    xlab("Days from symptom onset to death") + ylab("AT2 frequency")
ggsave("AT2_freq_vs_symptom_death_interval_nonimmune_norm.pdf", width = 7, height = 7)

# calculate combined frequencies of Pathological FB and Intermediate pathological
# FB, within the non-immune compartment
data_lungs_all$combined_pathological_fb = "no"
data_lungs_all$combined_pathological_fb[data_lungs_all$cell_type_fine %in% c("Pathological FB", 
    "Intermediate pathological FB")] = "yes"
df_tobesummed_fine = data.frame(orig.ident = data_lungs_all$orig.ident, cell_type_fine = data_lungs_all$cell_type_fine, 
    combined_pathological_fb = data_lungs_all$combined_pathological_fb, immune_status = data_lungs_all$immune_status)
df_tobesummed_fine = subset(df_tobesummed_fine, immune_status == "Non-immune")
df_tobesummed_fine$immune_status = NULL
df_summed_fine = df_tobesummed_fine %>% group_by(orig.ident, combined_pathological_fb) %>% 
    tally()
df_summed_fine = df_summed_fine %>% group_by(orig.ident) %>% mutate(freq = n/sum(n))

df_summed_fine$interval_death_symptoms_onset_days = 0
uniqueidents = unique(data_lungs_all$orig.ident)
for (i in 1:length(uniqueidents)) {
    onetime = unique(data_lungs_all$interval_death_symptoms_onset_days[data_lungs_all$orig.ident == 
        uniqueidents[i]])
    print(onetime)
    df_summed_fine$interval_death_symptoms_onset_days[df_summed_fine$orig.ident == 
        uniqueidents[i]] = onetime
}

df_summed_fine = subset(df_summed_fine, !(is.na(interval_death_symptoms_onset_days)))

# correlate combined Pathological FB and Intermediate pathological FB frequency
# with symptom-death interval
scatterplotdf = subset(df_summed_fine, combined_pathological_fb == "yes")
pearsoncorr = (cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$estimate)^2
pval = cor.test(scatterplotdf$freq, scatterplotdf$interval_death_symptoms_onset_days)$p.val
print("Pathological FB vs. interval")
print(pearsoncorr)
print(pval)
ggplot(scatterplotdf, aes(x = interval_death_symptoms_onset_days, y = freq)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.0449, p-val = 0.384"))) + 
    xlab("Days from symptom onset to death") + ylab("Pathological FB frequency")
ggsave("Pathological_FB_freq_vs_symptom_death_interval_nonimmune_norm.pdf", width = 7, 
    height = 7)
