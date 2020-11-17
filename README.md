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
8_Load_software_libraries | Utility script, loads R libraries necessary for downstream analyses. | Yiping Wang |  |  | -
9_Integrate_all_lung_samples_Seurat_2 | Seurat integration for lung samples (7 COVID-19 and 2 control donors) | Yiping Wang | RDS files of individual Seurat objects | RDS file of integrated Seurat object (9 donors) | -
10_Overall_celltype_classification_driver_script | Driver script that uses 13_All_samples_UMAP to manually label cells under broad categories, then 11_Add_gene_signature_module_scores to assign cell types. Prints out UMAP of cells, grouped by covid vs. control, and overall cell identity. Code for Figures 1C and 1D | Yiping Wang | RDS file of integrated Seurat object (9 donors) | RDS file of integrated Seurat object with cell type assignment (9 donors) | Figure 1C and 1D
11_Add_gene_signature_module_scores | Add gene signatures module scores for non-immune cell types. Use 12_Overall_celltype_classification twice to determine degree of module score overlap and assign cell type. After first assignment, derives a new 30-gene signature for tuft cells from initial cell assignment. Uses this signature to reassign cell types. | Yiping Wang |  | - | -
12_Overall_celltype_classification | Determine degree of overlap between module scores for different cell types. Assign overall celltypes if overlap is not too significant, according to criteria described in Methods. | Yiping Wang |  | - | -
13_All_samples_UMAP | Manually label cells in UMAP plot under broad categories. | Yiping Wang | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | -
14_Cell_fraction_analysis | Analysis of cell fractions in covid and control samples. Code for Figure 1E-H, Extended Data Figure 1A-D | Yiping Wang | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | Figure 1E-H, Extended Data Figure 1A-D
15_IG_combination_heatmaps | Analysis of IG combinations and frequencies across all samples, and for each individual sample. Code for Figure 3F. | Yiping Wang |  | -
16_AT_heatmap | Analysis of differential expression between AT1 and AT2 cells. Code for Figure 3C, D and E, Extended Data Table 5. | Yiping Wang | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | Figure 3C, D and E, Extended Data Table 5
17_Differential_expression_dotplots | Analysis of differential expression genes and pathways between covid and control across all cell types. Code for Figures 3H and 3I | Yiping Wang | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | Figures 3H and 3I
18_Tuft_cell_violin_plots | Analysis of DCLK1 expression in tuft vs. non-tuft cells, and tuft-1 versus tuft-2 signature strength in tuft cells. Code for Extended Data Figures 6A-C | Yiping Wang | RDS file of integrated Seurat object with cell type assignment (9 donors) | - | Extended Data Figures 6A-C
