#!/usr/bin/env nextflow

process BEDTOOLS_INTERSECT {
    label 'process_low'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input: 
    tuple val(sample_ids), path(bed_files)

    output: 
    path("reproducible_peaks.bed")

    shell:
    """
    bedtools intersect -a ${bed_files[0]} -b ${bed_files[1]} -f 0.50 -r > reproducible_peaks.bed
    """

    stub:
    """
    touch repr_peaks.bed
    """
}