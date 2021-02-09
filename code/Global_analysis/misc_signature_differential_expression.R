#!/usr/bin/env Rscript

### title: Calculate differential expression of three important gene signatures
### across all intermediate cell types author: Yiping Wang date: 02/08/2021

# Libraries
library(textTinyR)

DefaultAssay(data_lungs_all) = "RNA"

celltypestotest = sort(unique(data_lungs_all$cell_type_intermediate))

misc_sigs <- read.csv("misc_signatures.csv", na.strings = "")
for (c in 1:ncol(misc_sigs)) {
    misc_sigs[, c] = toupper(misc_sigs[, c])
    data_lungs_all = AddModuleScore(object = data_lungs_all, features = list(na.omit(misc_sigs[, 
        c])), name = colnames(misc_sigs)[c], assay = "RNA", search = T)
}
misc_sigs_arr = c("Type.I.interferon.abbreviated1", "inflammasome.receptors1", "chemotaxis1")

# for each signature and celltype in genes and celltypestotest, calculate it's
# wilcoxon p-value for differential expression between covid and ctr samples, and
# store in dotplotdf if p-val<0.1, store above information as well as log fold
# change in dotplotdf2 if p-val<.01, store above information in dotplotdf3
dotplotdf = matrix(nrow = length(misc_sigs_arr), ncol = length(celltypestotest))
dotplotdf2 = data.frame(tissue = character(), sig = character(), pval = double(), 
    logfoldchange = double())
dotplotdf3 = data.frame(tissue = character(), sig = character())
rownames(dotplotdf) = misc_sigs_arr
colnames(dotplotdf) = celltypestotest
for (i in 1:length(celltypestotest)) {
    for (j in 1:length(misc_sigs_arr)) {
        misc_sigs_genes = na.omit(misc_sigs[, j])
        datacol = eval(parse(text = paste0("data_lungs_all$", misc_sigs_arr[j])))
        datacol1 = data_lungs_all$RNA@data[rownames(data_lungs_all$RNA) %in% misc_sigs_genes, 
            data_lungs_all$cell_type_intermediate == celltypestotest[i] & data_lungs_all$group == 
                "COVID-19"]
        datacol2 = data_lungs_all$RNA@data[rownames(data_lungs_all$RNA) %in% misc_sigs_genes, 
            data_lungs_all$cell_type_intermediate == celltypestotest[i] & data_lungs_all$group == 
                "Control"]
        if (!is.null(dim(datacol1)) && !is.null(dim(datacol2))) {
            arr1 = sparse_Sums(datacol1, rowSums = F)
            arr2 = sparse_Sums(datacol2, rowSums = F)
            if (length(arr1) > 2 && length(arr2) > 2 && (min(arr1) != max(arr1) || 
                min(arr2) != max(arr2)) && sum(arr1) != 0 && sum(arr2) != 0) {
                arr1save = arr1
                arr2save = arr2
                wilc_pval = wilcox.test(arr1, arr2)$p.value * length(celltypestotest) * 
                  length(misc_sigs_arr)
                print(celltypestotest[i])
                print(misc_sigs_arr[j])
                print(wilc_pval)
                dotplotdf[j, i] = wilc_pval
                if (wilc_pval < 0.1) {
                  tempdf = data.frame(tissue = str_replace_all(celltypestotest[i], 
                    "_", " "), sig = str_replace_all(str_replace_all(misc_sigs_arr[j], 
                    "\\.", " "), "1", ""), pval = wilc_pval, logfoldchange = log(mean(arr1)/mean(arr2), 
                    2))
                  dotplotdf2 = rbind(dotplotdf2, tempdf)
                  if (wilc_pval < 0.01) {
                    tempdf = data.frame(tissue = str_replace_all(celltypestotest[i], 
                      "_", " "), sig = str_replace_all(str_replace_all(misc_sigs_arr[j], 
                      "\\.", " "), "1$", ""))
                    dotplotdf3 = rbind(dotplotdf3, tempdf)
                  }
                }
            }
        }
    }
}
write.table(dotplotdf, "misc_sigs_revision_results.csv", row.names = T, col.names = T, 
    sep = ",", quote = F)

theme_set(theme_bw())

dotplotdf2$pval = dotplotdf2$pval

ggplot(dotplotdf2, aes(x = tissue, y = sig, color = logfoldchange, size = pval)) + 
    geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradientn(colours = c("blue", 
    "#EDFF00FF", "red"), breaks = c(min(dotplotdf2$logfoldchange), 0, max(dotplotdf2$logfoldchange))) + 
    scale_size(trans = "reverse", breaks = c(0.005, 0.025, 0.05, 0.075), labels = c("<=.005", 
        ".025", ".05", ">=.075")) + ylab("Pathway") + xlab("Cell Type") + labs(color = "Log 2 FC")
ggsave(paste0("Extended Data Figure Cytokines a.pdf"), width = 11, height = 5)

write.table(dotplotdf2, "misc_sigs_revision_sig_results.csv", row.names = F, col.names = T, 
    sep = ",", quote = F)
