#!/usr/bin/env nextflow

process MULTIBWSUMMARY {
    label 'process_low'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path bigwig_files  
    val sample_labels

    output:
    path "results.npz", emit: npz

    script:
    def labels = sample_labels ? "--labels ${sample_labels.join(' ')}" : ""
    """
    # Compute summary of signals across the genome
    multiBigwigSummary bins \
        --bwfiles ${bigwig_files.join(' ')} \
        --binSize 10000 \
        --outFileName results.npz \
        ${labels}
    """

    stub:
    """
    touch bw_all.npz
    """
}