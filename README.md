# CUIMC-NYP_COVID_autopsy_lung

Code for the manuscript:
## A single-cell lung atlas of lethal COVID-19 
Johannes C. Melms, Jana Biermann, Huachao Huang, Yiping Wang, Amit Dipak Amin, Igor Katsyv, Yinshan Fang, Adrienne M. Luoma, Meri Rogava, Sean W. Chen, Patricia Ho, Adam E. Kornberg, Arnold S. Han, Jay H. Lefkowitch, Charles Marboe, Stephen. M. Lagana, Armando Del Portillo, Emmanuel Zorn, Glen S. Markowitz, Anjali Saqi, Hanina Hibshoosh, Jianwen  Que, Benjamin Izar


### Overview of code files
File name	| Title |	Authors |	Input |	Output | Figures/tables generated
--- | --- | --- | --- | --- | ---
1_Initial_Seurat_anlysis |	Seurat analysis for CellBender output with sample name and expected doublet rate provided in arguments |	Jana Biermann |	CellBender output (filtered.h5) |	RDS file of individual Seurat object |	-
2_Integrate_all_lung_samples_Seurat |	Seurat integration for lung samples (7 COVID-19 and 1 control donors; input for diffusion component analysis) |	Jana Biermann |	RDS files of individual Seurat objects |	RDS file of integrated Seurat object (8 donors) |	-
3_Cell_type_identification |	Identification of cell types in lung samples (7 COVID-19 and 1 control donors) |	Yiping Wang, Jana Biermann |	RDS file of integrated Seurat object (8 donors) |	RDS file of integrated Seurat object with cell type assignment (8 donors) |	-
4_Myeloid_cell_analysis |	Analysis of myeloid cells in lung samples (7 COVID-19 and 1 control samples) |	Jana Biermann | RDS file of integrated Seurat object with cell type assignment (8 donors) |	- | Figure 2abd, Extended Data Figure 2, Extended Data Table 2 and 3
5_B_cell_analysis |	Analysis of B cells in lung samples (7 COVID-19 and 1 control samples) | Jana Biermann, Yiping Wang |	RDS file of integrated Seurat object with cell type assignment (8 donors) |	- | Figure2efh, Extended Data Figure 3acdf, Extended Data Table 4
6_T_cell_analysis |	Analysis of T cells in lung samples (7 COVID-19 and 1 control samples) | Jana Biermann | RDS file of integrated Seurat object with cell type assignment (8 donors) |	- | Figure 2ij
7_AT_cell_analysis | Analysis of AT1 and AT2 cells in lung samples (7 COVID-19 and 1 control samples) | Jana Biermann |	RDS file of integrated Seurat object with cell type assignment (8 donors) |	- | Figure 3f, Extended Data Figure 4a-k
8_XXX |	Seurat integration for lung samples (7 COVID-19 and 2 control donors), fractions, tuft cell analysis, and DGE |	Yiping Wang |	RDS files of individual Seurat objects | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | Figure 1cdefgh, Extended Data Figure 1, Figure 3cdehi, Extended Data Figure 5, Extended Data Table 5
