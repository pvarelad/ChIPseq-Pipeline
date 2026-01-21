#!/usr/bin/env nextflow

process SAMTOOLS_IDX {
    label 'process_single'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sorted_bam}.bai"), emit: bai

    script:
    """
    samtools index ${sorted_bam}
    """

    stub:
    """
    touch ${sample_id}.stub.sorted.bam.bai
    """
}