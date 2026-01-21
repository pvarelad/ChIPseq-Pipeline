#!/usr/bin/env nextflow

process PLOTPROFILE {
    label 'process_high'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(matrix)

    output:
    tuple val(sample_id), path("${sample_id}_signal_coverage.png")

    shell:
    """
    plotProfile -m ${matrix} -o ${sample_id}_signal_coverage.png --plotTitle "${sample_id} Signal Coverage" --perGroup
    """

    stub:
    """
    touch ${sample_id}_signal_coverage.png
    """
}