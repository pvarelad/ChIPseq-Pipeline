#!/usr/bin/env nextflow

process TRIM {
    label 'process_medium'
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    path adapters_fa
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed
    tuple val(sample_id), path("${sample_id}_trim.log"), emit: log

    shell: 
    """
    trimmomatic SE -threads ${task.cpus} \
        ${reads} \
        ${sample_id}_trimmed.fastq.gz \
        ILLUMINACLIP:${adapters_fa}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${sample_id}_trim.log
    """


    stub:
    """
    touch ${sample_id}_stub_trim.log
    touch ${sample_id}_stub_trimmed.fastq.gz
    """
}
