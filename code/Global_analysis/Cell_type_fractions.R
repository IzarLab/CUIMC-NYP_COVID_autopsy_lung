#!/usr/bin/env Rscript

#### Tabels for cell type fractions and SARS-CoV-2 reads
library(Seurat)
library(dplyr)

# Load data
seu <- readRDS('data/lungs_all/data_lungs_all.rds')

# Total fractions by cell types
write.csv(prop.table(table(seu$cell_type_main)), 'data/lungs_all/table_freq_cell_type_main.csv', row.names = F)
write.csv(prop.table(table(seu$cell_type_intermediate)), 'data/lungs_all/table_freq_cell_type_intermediate.csv', row.names = F)
write.csv(prop.table(table(seu$cell_type_fine)), 'data/lungs_all/table_freq_cell_type_fine.csv', row.names = F)

# Main grouped by disease
write.csv(prop.table(table(subset(seu, group == 'COVID-19')$cell_type_main)), 'data/lungs_all/table_freq_cell_type_main_cov.csv', row.names = F)
write.csv(prop.table(table(subset(seu, group == 'Control')$cell_type_main)), 'data/lungs_all/table_freq_cell_type_main_ctr.csv', row.names = F)

# Intermediate grouped by disease
write.csv(prop.table(table(subset(seu, group == 'COVID-19')$cell_type_intermediate)), 'data/lungs_all/table_freq_cell_type_intermediate_cov.csv', row.names = F)
write.csv(prop.table(table(subset(seu, group == 'Control')$cell_type_intermediate)), 'data/lungs_all/table_freq_cell_type_intermediate_ctr.csv', row.names = F)

# Fine grouped by disease
write.csv(prop.table(table(subset(seu, group == 'COVID-19')$cell_type_fine)), 'data/lungs_all/table_freq_cell_type_fine_cov.csv', row.names = F)
write.csv(prop.table(table(subset(seu, group == 'Control')$cell_type_fine)), 'data/lungs_all/table_freq_cell_type_fine_ctr.csv', row.names = F)

# Count of SARS-CoV-2 reads in patient L19cov
L19 <- subset(seu, patient == 'L19cov')
hist(L19@assays$RNA@counts['SARS-COV-2', ], breaks = 200, ylim = c(0, 10))
dfL19 <- data.frame(SARS_COV_2 = L19@assays$RNA@counts['SARS-COV-2', ], barcode = colnames(L19@assays$RNA@counts), cell_type_fine = L19$cell_type_fine)
write.csv(table(dfL19$SARS_COV_2, dfL19$cell_type_fine), 'data/lungs_all/table_SARS-CoV-2_reads_in_L19cov.csv')
