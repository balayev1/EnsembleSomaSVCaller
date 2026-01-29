#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include the new process
include { ASCAT } from './nextflow/ascat.nf'
include { SEQUENZAUTILS_GCWIGGLE } from './nextflow/sequenza.nf'
include { SEQUENZAUTILS_BAM2SEQZ } from './nextflow/sequenza.nf'
include { SEQUENZAUTILS_BINNING }  from './nextflow/sequenza.nf'
include { SEQUENZA_RUN }           from './nextflow/sequenza.nf'

workflow {
    // Prepare global reference files from params
    fasta      = params.fasta      ? file(params.fasta)      : []
    allele_res = params.allele_res ? file(params.allele_res) : []
    loci_res   = params.loci_res   ? file(params.loci_res)   : []
    gc_file    = params.gc_file    ? file(params.gc_file)    : []
    rt_file    = params.rt_file    ? file(params.rt_file)    : []
    bed_file   = params.bed_file   ? file(params.bed_file)   : []

    // Gamma range for Sequenza
    ch_gammas  = Channel.fromList(params.gamma ?: [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000])

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
    /*
    ASCAT (
        ch_samples,
        allele_res,
        loci_res,
        bed_file,
        fasta,
        gc_file,
        rt_file
    )
    */

    // Run the Sequenza workflow
    // Generate GC Wiggle Reference
    // We use .first() so it only runs once and provides the file to all samples
    ch_wiggle_ref = SEQUENZAUTILS_GCWIGGLE( [ [id:'genome_ref'], fasta ] ).wig.map{ it[1] }.first()

    // Convert BAM to SEQZ
    ch_bam_pairs = ch_samples.map { meta, n_bam, n_bai, t_bam, t_bai -> [ meta, n_bam, t_bam ] }
    
    SEQUENZAUTILS_BAM2SEQZ (
        ch_bam_pairs,
        fasta,
        ch_wiggle_ref
    )

    // Binning SEQZ files
    SEQUENZAUTILS_BINNING ( SEQUENZAUTILS_BAM2SEQZ.out.seqz )

    // Run Sequenza (Sample x Gamma)
    // .combine creates the cartesian product of samples and gammas
    ch_run_input = SEQUENZAUTILS_BINNING.out.binned_seqz.combine(ch_gammas)

    SEQUENZA_RUN ( ch_run_input )
}