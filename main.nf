#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include the new process
include { ASCAT } from './modules/ascat.nf'
include { SEQUENZAUTILS_GCWIGGLE } from './modules/sequenza.nf'
include { SEQUENZA_RUN } from './modules/sequenza.nf'
include { FACETS } from './modules/facets.nf'

workflow {
    // Prepare global reference files from params
    fasta = params.fasta ? [file(params.fasta), file("${params.fasta}.fai")] : []
    allele_res = params.allele_res ? file(params.allele_res) : []
    loci_res   = params.loci_res   ? file(params.loci_res)   : []
    gc_file    = params.gc_file    ? file(params.gc_file)    : []
    rt_file    = params.rt_file    ? file(params.rt_file)    : []
    bed_file   = params.bed_file   ? file(params.bed_file)   : []

    // Setup Input Channel
    ch_samples = Channel.fromPath(params.input)
        .splitCsv(sep: '\t')
        .map { row -> 
            def meta = [ id: row[0], gender: row[5] ]
            
            // This ensures we find 'sample.bai' instead of 'sample.bam.bai'
            def normal_bam = file(row[1])
            def normal_bai = file(row[2])
            def tumor_bam  = file(row[3])
            def tumor_bai  = file(row[4])

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

    // Run the SEQUENZA process
    // Generate GC Wiggle Reference
    // We use .first() so it only runs once and provides the file to all samples
    ch_wiggle_ref = SEQUENZAUTILS_GCWIGGLE( [ [id:'genome_ref'], fasta ] ).wig.map{ it[1] }.first()

    // Run SEQUENZA
    SEQUENZA_RUN ( 
        ch_samples, 
        fasta, 
        ch_wiggle_ref 
        )

    // Run the FACETS process
    // Set FACETS-specific parameters
    // VCF file with common polymoprhisms
    if (!params.facets_snp_vcf) {
        error "FACETS requires a SNP VCF file! Please specify 'params.facets_snp_vcf' in your config."
    }
    snp_vcf = file(params.facets_snp_vcf, checkIfExists: true)
    idx_tbi = file("${params.facets_snp_vcf}.tbi")
    idx_csi = file("${params.facets_snp_vcf}.csi")
    snp_vcf_index = idx_tbi.exists() ? idx_tbi : (idx_csi.exists() ? idx_csi : null)
    if (!snp_vcf_index) {
        error "Could not find an index file for ${snp_vcf}. \nExpected: ${snp_vcf}.tbi OR ${snp_vcf}.csi"
    }

    // Optional target and annotation BED files
    targets_bed = params.facets_targets_bed ? file(params.facets_targets_bed) : []
    annotation_bed = params.facets_annotation_bed ? file(params.facets_annotation_bed) : []
    
    // Run FACETS
    FACETS (
        ch_samples,
        snp_vcf,
        snp_vcf_index,
        targets_bed,
        annotation_bed
    )


}