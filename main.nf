#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*.fasta'
params.outdir = 'outputs/'
params.adapters = 'adapters.fa'
log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
"""

// Create read channel
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)

// Define fastqc process
process fastqc {
    label 'fastqc'
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    fastqc ${reads}
    """
}

// Process trimmomatic
process trimmomatic {
    label 'trimmomatic'
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    tuple val("${sample}"), path("${sample}*.trimmed.fq.gz")
    tuple val("${sample}"), path("${sample}*.discarded.fq.gz")

    script:
    """
    trimmomatic PE -phred33 ${reads[0]} ${reads[1]} ${sample}_1.trimmed.fq.gz ${sample}_1.discarded.fq.gz ${sample}_2.trimmed.fq.gz ${sample}_2.discarded.fq.gz
    ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
}
