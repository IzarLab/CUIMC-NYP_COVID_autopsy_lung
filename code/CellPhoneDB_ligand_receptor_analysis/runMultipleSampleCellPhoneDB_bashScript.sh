#!/bin/bash
#$ -o output_cellPhoneDB_bash.txt -e error_cellPhoneDB_bash.txt

#giving the sample names to run
#All samples: L01cov L03cov L04cov L04covaddon L05cov L06cov L07cov L08cov L09cov L10cov L11cov L12cov L13cov L15cov L16cov L17cov L18cov L19cov L21cov L22cov

SAMPLE_NAME_ARRAY=(L01cov L03cov L04cov L04covaddon L05cov L06cov) #running a subset of samples
#doing the CellPhoneDB runs
echo "  "
echo "**************************"
echo "Running CellPhoneDB  $(date)"
for SAMPLE in ${SAMPLE_NAME_ARRAY[*]}; do
    echo "  "
    echo "**************************"
    echo "Patient:  ${SAMPLE}"
    cd /PATH/TO/DATA/FOLDER/Covid/${SAMPLE}
    cellphonedb method statistical_analysis data_cellphoneDB_metaData_hum_${SAMPLE}.txt data_cellphoneDB_countNorm_hum_${SAMPLE}.txt --output-path=out_cellTypeMain --counts-data=gene_name --threads=3
    #heatmap
    cellphonedb plot heatmap_plot --pvalues-path out_cellTypeMain/pvalues.txt --output-path out_cellTypeMain --count-name heatmap_cellTypeMain_${SAMPLE}.pdf data_cellphoneDB_metaData_hum_${SAMPLE}.txt
    #plotting
    cellphonedb plot dot_plot --means-path out_cellTypeMain/means.txt --pvalues-path out_cellTypeMain/pvalues.txt --output-path out_cellTypeMain --output-name dotplot_cellTypeMain_${SAMPLE}.pdf
done

echo "**************************"
echo "Finished runs CellPhoneDB  $(date)"

