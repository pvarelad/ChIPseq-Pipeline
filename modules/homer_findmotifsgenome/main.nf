#!/usr/bin/env nextflow

process FIND_MOTIFS_GENOME {
    label 'process_medium'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'

    input: 
    path(bed)
    path(fasta)

    output:
    path("motifs")

    shell:
    """
    mkdir -p custom_genome
    
    loadGenome.pl -name hg38 -fasta ${fasta} -org human -directory custom_genome
    
    findMotifsGenome.pl \
        ${bed} \
        ${fasta} \
        motifs \
        -size 200 \
        -mask \
        -preparsedDir custom_genome
    """

    stub:
    """
    mkdir motifs
    """
}


