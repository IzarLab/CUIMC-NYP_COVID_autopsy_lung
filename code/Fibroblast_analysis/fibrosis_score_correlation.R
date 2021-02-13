#!/usr/bin/env Rscript

### title: Print out scatterplots of 1. Fibrosis score in each patient vs. symptom
### to death interval, 2. Age vs. symptom to death interval, 3. Fibrosis score vs.
### age, for covid samples only 4. Fibrosis score vs. age, for control samples only
### author: Yiping Wang date: 02/08/2021

theme_set(theme_bw())

fibrosis_table = read.table("Fibrosis_score_sequencing_cohort.txt", header = T, sep = "\t", 
    quote = NULL)

fibrosis_table = fibrosis_table[!is.na(fibrosis_table$Sirius.Red) & !is.na(fibrosis_table$symptoms.death..days.) & 
    fibrosis_table$Group == "COVID-19", ]
pearsoncorr = (cor.test(fibrosis_table$Sirius.Red, fibrosis_table$symptoms.death..days.)$estimate)^2
pval = cor.test(fibrosis_table$Sirius.Red, fibrosis_table$symptoms.death..days.)$p.val
print(pearsoncorr)
print(pval)
ggplot(fibrosis_table, aes(x = symptoms.death..days., y = Sirius.Red)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.386, p-val = 0.0103"))) + 
    xlab("Days from symptom onset to death") + ylab("Fibrosis Score")
ggsave("Figure_4c.pdf", width = 7, height = 7)

fibrosis_table = read.table("Fibrosis_score_sequencing_cohort.txt", header = T, sep = "\t", 
    quote = NULL)

fibrosis_table = fibrosis_table[!is.na(fibrosis_table$Age) & !is.na(fibrosis_table$symptoms.death..days.) & 
    fibrosis_table$Group == "COVID-19", ]
pearsoncorr = (cor.test(fibrosis_table$Age, fibrosis_table$symptoms.death..days.)$estimate)^2
pval = cor.test(fibrosis_table$Age, fibrosis_table$symptoms.death..days.)$p.val
print(pearsoncorr)
print(pval)
ggplot(fibrosis_table, aes(x = symptoms.death..days., y = Age)) + geom_point() + 
    geom_smooth(method = "lm") + ggtitle(expression(paste(" ", R^2, "= 0.279, p-val = 0.0293"))) + 
    xlab("Days from symptom onset to death") + ylab("Age")
ggsave("Extended_Data_Fibroblasts_C.pdf", width = 7, height = 7)

fibrosis_table = read.table("Fibrosis_score_sequencing_cohort.txt", header = T, sep = "\t", 
    quote = NULL)

fibrosis_table = fibrosis_table[!is.na(fibrosis_table$Sirius.Red) & !is.na(fibrosis_table$Age) & 
    fibrosis_table$Group == "COVID-19", ]
pearsoncorr = (cor.test(fibrosis_table$Sirius.Red, fibrosis_table$Age)$estimate)^2
pval = cor.test(fibrosis_table$Sirius.Red, fibrosis_table$Age)$p.val
print(pearsoncorr)
print(pval)
ggplot(fibrosis_table, aes(x = Age, y = Sirius.Red)) + geom_point() + geom_smooth(method = "lm") + 
    ggtitle(expression(paste(" ", R^2, "= 0.177, p-val = 0.105"))) + xlab("Age") + 
    ylab("Fibrosis Score")
ggsave("Extended_Data_Fibroblasts_D.pdf", width = 7, height = 7)

fibrosis_table = read.table("Fibrosis_score_sequencing_cohort.txt", header = T, sep = "\t", 
    quote = NULL)

fibrosis_table = fibrosis_table[!is.na(fibrosis_table$Sirius.Red) & !is.na(fibrosis_table$Age) & 
    fibrosis_table$Group == "Control", ]
pearsoncorr = (cor.test(fibrosis_table$Sirius.Red, fibrosis_table$Age)$estimate)^2
pval = cor.test(fibrosis_table$Sirius.Red, fibrosis_table$Age)$p.val
print(pearsoncorr)
print(pval)
ggplot(fibrosis_table, aes(x = Age, y = Sirius.Red)) + geom_point() + geom_smooth(method = "lm") + 
    ggtitle(expression(paste(" ", R^2, "= 0.0122, p-val = 0.813"))) + xlab("Age") + 
    ylab("Fibrosis Score")
ggsave("Extended_Data_Fibroblasts_E.pdf", width = 7, height = 7)
