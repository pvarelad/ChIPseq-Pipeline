#!/usr/bin/env nextflow


    process PLOTCORRELATION {
    label 'process_low'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path npz_file       

    output:
    path "correlation_plot.png", emit: plot
    path "correlation_values.tab", emit: data

    script:
    """
    # Generate correlation plot
    plotCorrelation \
        --corData ${npz_file} \
        --corMethod pearson \
        --skipZeros \
        --whatToPlot heatmap \
        --plotFile correlation_plot.png \
        --plotNumbers \
        --outFileCorMatrix correlation_values.tab
    """
}






