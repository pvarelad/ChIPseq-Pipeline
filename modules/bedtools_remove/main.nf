#!/usr/bin/env nextflow

process BEDTOOLS_REMOVE {
    label 'process_low'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path(reproducible_peaks)
    path(blacklist)

    output:
    path("reproducible_peaks_filtered.bed")

    shell:
    """
    bedtools intersect -a ${reproducible_peaks} -b ${blacklist} -v > reproducible_peaks_filtered.bed
    """

    stub:
    """
    touch repr_peaks_filtered.bed
    """
}