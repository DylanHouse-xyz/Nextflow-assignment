#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq.gz'
params.outdir = './outputs/'
params.adapters = 'data/adapters.fa'
params.indexDir = 'data/'
params.indexFiles = '*data/LG12.{0123, amb, ann, bwt.2bit.64, pac}'
params.refgenome = 'data/*.fasta'
log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
IndexFiles       : ${params.indexFiles}
IndexDir         : ${params.indexDir}
Ref_Genome	 : ${params.refgenome}
"""

// Create read channel
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)
genome_ch = Channel.fromPath(params.refgenome)
index_ch = Channel.fromPath(params.indexFiles).collect()

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
    tuple val("${sample}"), path("${sample}*.trimmed.fq.gz"), emit: trimmed_fq
    tuple val("${sample}"), path("${sample}*.discarded.fq.gz"), emit: discarded_fq

    script:
    """
    trimmomatic PE -phred33 ${reads[0]} ${reads[1]} ${sample}_1.trimmed.fq.gz ${sample}_1.discarded.fq.gz ${sample}_2.trimmed.fq.gz ${sample}_2.discarded.fq.gz ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

// Process BWA_MEM
process bwa_mem {
    label "bwa_mem"
    publishDir "${params.outdir}/bwa-{sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path index_file
    //val reference_genome

    output:
    path("${sample}_bwamem2.sam")

    script:
    """
    bwa-mem2 mem -t4 ${index_file} ${reads[0]} ${reads[1]} > bwamem2.sam
    """
}

// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
    bwa_mem(trimmomatic.out.trimmed_fq, index_ch)
}
