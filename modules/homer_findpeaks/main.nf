#!/usr/bin/env nextflow

process FINDPEAKS {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(ip_tagdir), path(input_tagdir)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.txt")

    shell:
    """
    findPeaks ${ip_tagdir} -i ${input_tagdir} -style factor -o ${sample_id}_peaks.txt
    """
    stub:
    """
    touch ${rep}_peaks.txt
    """
}


