#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    VALIDATE INPUTS
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.allele_res, params.loci_res, params.gc_file, params.rt_file, params.facets_snp_vcf ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
    CHANNEL SETUP
*/
// Universal parameters
refgen                  = Channel.fromPath([params.fasta,params.fasta_fai], checkIfExists: true).collect()

// ASCAT parameters
allele_res              = Channel.fromPath(params.allele_res, checkIfExists: true)
loci_res                = Channel.fromPath(params.loci_res, checkIfExists: true)
gc_file                 = Channel.fromPath(params.gc_file, checkIfExists: true)
rt_file                 = Channel.fromPath(params.rt_file, checkIfExists: true)
ascat_bed_file          = params.bed_file                       ? Channel.fromPath(params.bed_file, checkIfExists: true)
                                                                : Channel.empty()

// FACETS parameters
facets_snp_vcf          = Channel.fromPath([params.facets_snp_vcf, params.facets_snp_vcf + '.tbi'], checkIfExists: true).collect()
facets_targets_bed      = params.facets_targets_bed             ? Channel.fromPath(params.facets_targets_bed, checkIfExists: true) 
                                                                : Channel.empty()
facets_annotation_bed   = params.facets_annotation_bed          ? Channel.fromPath(params.facets_annotation_bed, checkIfExists: true)
                                                                : Channel.empty()

// Hetpileups parameters
hapmap_sites       = params.hapmap_sites       ? Channel.fromPath(params.hapmap_sites).collect()      : Channel.empty()
filter_hets         = params.filter_hets       ?: Channel.empty()
max_depth           = params.max_depth         ?: Channel.empty()

// fragCounter parameters
gcmapdir_frag      = params.gcmapdir_frag      ? Channel.fromPath(params.gcmapdir_frag).collect() : Channel.empty()
windowsize_frag    = params.windowsize_frag    ?: Channel.empty()
minmapq_frag       = params.minmapq_frag       ?: Channel.empty()
midpoint_frag      = params.midpoint_frag      ?: Channel.empty()
paired_frag        = params.paired_frag        ?: Channel.empty()
exome_frag         = params.exome_frag         ?: Channel.empty()

// dryclean parameters
pon_dryclean           = params.pon_dryclean           ? Channel.fromPath(params.pon_dryclean).collect() : Channel.empty()
center_dryclean        = params.center_dryclean        ?: Channel.empty()
cbs_dryclean           = params.cbs_dryclean           ?: Channel.empty()
cnsignif_dryclean      = params.cnsignif_dryclean      ?: Channel.empty()
wholeGenome_dryclean   = params.wholeGenome_dryclean   ?: Channel.empty()
blacklist_dryclean     = params.blacklist_dryclean     ?: Channel.empty()
blacklist_path_dryclean= params.blacklist_path_dryclean? Channel.fromPath(params.blacklist_path_dryclean).collect() : Channel.empty()
germline_filter_dryclean= params.germline_filter_dryclean ?: Channel.empty()
germline_file_dryclean = params.germline_file_dryclean ? Channel.fromPath(params.germline_file_dryclean).collect() : Channel.empty()
field_dryclean         = params.field_dryclean         ?: Channel.empty()
build_dryclean         = params.build_dryclean         ?: Channel.empty()

/*
    IMPORT SUBWORKFLOWS
*/
include { INPUT_PREP } from './subworkflows/input_prep.nf'
include { ZERO_SHOT_CNV_CALL } from './subworkflows/zeroshot_cnv_calling.nf'
include { BAM_HETPILEUPS } from './subworkflows/bam_hetpileups/main.nf'

/*
    MAIN WORKFLOW
*/
workflow SOMASV_CALLER {
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_samples = INPUT_PREP(ch_input)

    //println "The samples: "
    ch_samples.view()

    //
    // SUBWORKFLOW: Run zero-shot CNV calling with ASCAT, SEQUENZA, and FACETS
    //  
    ZERO_SHOT_CNV_CALL (
        ch_samples, 
        refgen, 
        allele_res, 
        loci_res, 
        gc_file, 
        rt_file, 
        ascat_bed_file,
        facets_snp_vcf,
        facets_targets_bed,
        facets_annotation_bed
    )
    ch_versions    = ch_versions.mix(ZERO_SHOT_CNV_CALL.out.versions)

    //
    // SUBWORKFLOW: Generate GC- and mappability-corrected denoised genomic coverage files for JAbBA
    //
    COVERAGE_JABBA (
        ch_samples, 
        midpoint_frag, 
        windowsize_frag, 
        gcmapdir_frag, 
        minmapq_frag, 
        paired_frag, 
        exome_frag, 
        pon_dryclean, 
        center_dryclean, 
        cbs_dryclean, 
        cnsignif_dryclean, 
        wholeGenome_dryclean, 
        blacklist_dryclean, 
        blacklist_path_dryclean, 
        germline_filter_dryclean, 
        germline_file_dryclean, 
        field_dryclean, 
        build_dryclean
    )
    ch_versions    = ch_versions.mix(COVERAGE_JABBA.out.versions)

    //
    // SUBWORKFLOW: Generate heterozygous pileups file for JabBa
    //
    BAM_HETPILEUPS (
        ch_samples, 
        filter_hets, 
        max_depth,
        hapmap_sites
    )
    ch_versions = ch_versions.mix(BAM_HETPILEUPS.out.versions)
    sites_from_het_pileups_wgs = Channel.empty().mix(BAM_HETPILEUPS.out.het_pileups_wgs)
}