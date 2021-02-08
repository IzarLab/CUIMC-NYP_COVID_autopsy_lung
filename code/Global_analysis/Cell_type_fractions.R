#!/usr/bin/env Rscript

#### Tabels for cell type fractions
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