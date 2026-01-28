#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include the new process
include { ASCAT } from './nextflow/ascat.nf'

workflow {
    // Prepare global reference files from params
    fasta      = params.fasta      ? file(params.fasta)      : []
    allele_res = params.allele_res ? file(params.allele_res) : []
    loci_res   = params.loci_res   ? file(params.loci_res)   : []
    gc_file    = params.gc_file    ? file(params.gc_file)    : []
    rt_file    = params.rt_file    ? file(params.rt_file)    : []
    bed_file   = params.bed_file   ? file(params.bed_file)   : []

    // Setup Input Channel
    ch_samples = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [ id: row.subjectID, gender: row.gender ]
            
            // This ensures we find 'sample.bai' instead of 'sample.bam.bai'
            def normal_bam = file(row.normalBAM)
            def normal_bai = file(row.normalBAM.replaceFirst(/\.bam$/, ".bai"))
            def tumor_bam  = file(row.tumorBAM)
            def tumor_bai  = file(row.tumorBAM.replaceFirst(/\.bam$/, ".bai"))

            return [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
        }

    // Run the ASCAT process
    ASCAT (
        ch_samples,
        allele_res,
        loci_res,
        bed_file,
        fasta,
        gc_file,
        rt_file
    )
}