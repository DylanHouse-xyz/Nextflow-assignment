#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = './data/fastas/*_{1,2}.fq.gz'
params.outdir = './outputs/'
params.adapters = './data/adapters.fa'
params.indexDir = './data/index'
params.indexFiles = './data/index*.{0123,amb,ann,bwt.2bit.64,pac}'
params.refgenome = './data/index/LG12.fasta'
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

include { fastqc } from './modules/fastqc'
include { trimmomatic } from './modules/trimmomatic'
include { bwa_mem } from './modules/bwa_mem'

// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
    bwa_mem(trimmomatic.out.trimmed_fq, file(params.refgenome).toAbsolutePath())
}
