#!/usr/bin/env nextflow

// Process BWA_MEM
process bwa_mem {
    label "bwa_mem"
    publishDir "${params.outdir}/bwa-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val index_file
    //val reference_genome

    output:
    path("${sample}_bwamem2.sam")

    script:
    """
    bwa-mem2 mem -t4 ${index_file} ${reads[0]} ${reads[1]} > "${sample}_bwamem2.sam"
    """
}

