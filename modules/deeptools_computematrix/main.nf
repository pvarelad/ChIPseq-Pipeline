#!/usr/bin/env nextflow

process COMPUTEMATRIX {
    label 'process_high'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val (sample_id), path(bw)
    path(ucsc_genes)

    output:
    tuple val(sample_id), path("${sample_id}_matrix.gz")

    shell:
    """
    computeMatrix scale-regions \
        -S ${bw} \
        -R ${ucsc_genes} \
        -o ${sample_id}_matrix.gz \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --skipZeros \
        -p ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}_matrix.gz
    """
}