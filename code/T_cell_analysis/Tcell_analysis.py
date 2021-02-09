# Functions relating to analyses of T cell data

# Import necessary libraries / functions
import numpy as np
import pandas as pd
import seaborn as sns
import os,glob
import scipy.sparse as sp_sparse
import scipy.stats as sp_stats
from pdb import set_trace as bp
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from harmony import harmonize

# Names of T cell signatures
sig_names = ['ACT', 'EXH', 'MEM', 'TRM']

# Run Harmony batch correction
# Inputs:
#   adata : AnnData structure containing single cell expression and metadata
#       'X_pca' field in adata.obsm is PCA on expression matrix
#       'orig.ident' field in adata.obs DataFrame is patient ID used as batch key
# Output:
#   Batch-corrected PCA matrix
def run_Harmony(adata, batch_key='orig.ident'):
    return harmonize(adata.obsm['X_pca'],adata.obs,batch_key=batch_key)

# Determine perecent of TRM cells as a function of TRM score threshold
# Inputs:
#   TRM_scores : array of length (number of cells) with TRM program scores
#   control_inds : indices of control cells
#   covid_inds : indices of covid cells
#   start_sweep : value to start sweep of TRM threshold
#   end_sweep : value to end sweep of TRM threshold
#   sweep_step : step size for sweep of TRM threshold
# Outputs:
#   num_control_cells : number of control cells classified as TRM as a function of TRM threshold
#   num_covid_cells : number of covid cells classified as TRM as a function of TRM threshold
#   TRM_score_thresh_arr : TRM score thresholds tested
def vary_TRM_threshold(TRM_scores, control_inds, covid_inds, start_sweep=0.3, end_sweep=1.1, sweep_step=0.01):

    TRM_score_thresh_arr = np.arange(start_sweep, end_sweep, sweep_step)

    num_control_cells = np.zeros(TRM_score_thresh_arr.size)
    num_covid_cells = np.zeros(TRM_score_thresh_arr.size)

    # Loop through TRM scores thresholds + calculate number TRM cells at each threshold
    for TRM_i, TRM_score_thresh in enumerate(TRM_score_thresh_arr):
        num_control_cells[TRM_i] = np.intersect1d(control_inds, np.where(TRM_scores > TRM_score_thresh)[0]).size
        num_covid_cells[TRM_i] = np.intersect1d(covid_inds, np.where(TRM_scores > TRM_score_thresh)[0]).size

    return TRM_score_thresh_arr,num_control_cells, num_covid_cells

# Plot figure relating to number of TRM cells as a function of TRM score threshold
# Inputs:
#   TRM_score_thresh_arr : array of TRM thresholds tested in function vary_TRM_threshold()
#   num_control_cells : number of control cells classified as TRM as a function of TRM threshold
#   num_covid_cells : number of covid cells classified as TRM as a function of TRM threshold
#   control_inds : indices of control cells
#   covid_inds : indices of covid cells
# Output:
#   Saves pdf of figure in current working directory
def plot_TRM_cells_fig(TRM_score_thresh_arr, num_control_cells, num_covid_cells, control_inds, covid_inds):

    fig, ax = plt.subplots()
    ax.plot(TRM_score_thresh_arr, 100*(num_control_cells/control_inds.size), label='Control', color='#006E82')
    ax.plot(TRM_score_thresh_arr, 100*(num_covid_cells/covid_inds.size), label='COVID', color='#AA0A3C')
    ax.legend()
    ax.set_xlabel('TRM Score Threshold'), ax.set_ylabel('Percent TRM Cells')
    fig.savefig('./Fig_percent_TRM_cells.pdf', dpi=800, bbox_inches='tight')
    plt.close()

# Plot boxplots of T cell signature scores
# Inputs:
#   sig_pn : path name of directory to save signature scores
#   cell_type_arr : array of length (number of cells) with cell type assignments
#   control_inds : indices of control cells
#   covid_inds : indices of covid cells
# Output:
#   Saves figure of program scores of each T cell program in current working directoryß
#   Prints p-value of covid vs control for CD4 and CD8 T cells for each program
def plot_Tcell_boxplots(sig_pn, cell_type_arr, control_inds, covid_inds):

    # Find indices of CD4 and CD8 cells in covid and control
    CD4_inds = np.where(cell_type_arr=='CD4+ T cells')[0]
    CD8_inds = np.where(cell_type_arr=='CD8+ T cells')[0]
    CD4_covid = np.intersect1d(CD4_inds, covid_inds)
    CD4_control = np.intersect1d(CD4_inds, control_inds)
    CD8_covid = np.intersect1d(CD8_inds, covid_inds)
    CD8_control = np.intersect1d(CD8_inds, control_inds)

    color_dict = {'Control' : '#006E82', 'COVID' : '#AA0A3C'}
    for sig_name in sig_names:

        prog_scores_up = np.load(os.path.join(sig_pn, '%s_UP.Tcell_scores.npy' % sig_name))
        prog_scores_down = np.load(os.path.join(sig_pn, '%s_DOWN.Tcell_scores.npy' % sig_name))
        prog_scores_composite = prog_scores_up - prog_scores_down # Composite score is up-regulation score minus down-regulation score

        plot_df_comp = pd.DataFrame()
        plot_df_comp = plot_df_comp.append(pd.DataFrame({'Score' : prog_scores_composite, 'Type' : 'Composite', 'Cell Type' : cell_type_arr, 'Disease' : cov_label_arr}), ignore_index=True)

        fig, ax = plt.subplots()
        sns.boxplot(data=plot_df_comp, x='Cell Type', y='Score', hue='Disease', ax=ax, palette=color_dict)
        ax.set_xlabel(''), ax.set_title('%s Composite Score' % sig_name)
        ax.set_ylabel('Program Score')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
        fig.savefig('./%s_boxplot_celltypes.pdf' % sig_name, dpi=800, bbox_inches='tight')
        plt.close()

        # Print statistics by Welch's t-test
        print('%s CD4 p = %.3E' % (sig_name, sp_stats.ttest_ind(prog_scores_composite[CD4_control], prog_scores_composite[CD4_covid], equal_var=False).pvalue))
        print('%s CD8 p = %.3E' % (sig_name, sp_stats.ttest_ind(prog_scores_composite[CD8_control], prog_scores_composite[CD8_covid], equal_var=False).pvalue))

# Score T cell programs
# Inputs:
#   adata : AnnData structure containing single cell expression and metadata
#   sig_pn : path name of directory to save signature scores
#   sig_df : DataFrame of program scores with column headers as program names and entries as scores by cell
#   num_bins : number of expression bins for scoring
# Output:
#   Saves program scores for each program in directory sig_pn
def score_Tcell_progs(adata, sig_pn, sig_df, num_bins=50):

    # Loop through each program
    for sig_name in sig_names:

        curr_genes = sig_df[sig_name].to_numpy().astype(np.str)
        curr_genes = np.delete(curr_genes, np.where(curr_genes=='nan')[0]) # Ignore blank spaces in DataFrame

        sc.tl.score_genes(adata, curr_genes, ctrl_size=curr_genes.size, n_bins=num_bins, score_name=sig_name) # Score program using Scanpy implementation of Seurat gene set scoring

        np.save(os.path.join(sig_pn, '%s.Tcell_scores.npy' % sig_name), adata.obs[sig_name].to_numpy()) # Save program scores
