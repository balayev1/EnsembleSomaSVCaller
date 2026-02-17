#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    VALIDATE INPUTS
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, 
    params.fasta, 
    params.fasta_fai, 
    params.ascat_alleles, 
    params.ascat_loci, 
    params.ascat_gc, 
    params.ascat_rt, 
    params.dbsnp, 
    params.dbsnp_tbi, 
    params.gcmapdir_frag,
    params.hapmap_sites,
    params.pon_dryclean
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
    CHANNEL SETUP
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
fasta               = params.fasta              ? Channel.fromPath(params.fasta).first()                : Channel.empty()
fasta_fai           = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()          : Channel.empty()
target_regs         = params.target_regs        ? Channel.fromPath(params.target_regs).collect()         : Channel.empty()
dbsnp               = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()              : Channel.empty()
dbsnp_tbi               = params.dbsnp_tbi              ? Channel.fromPath(params.dbsnp_tbi).collect()              : Channel.empty()

// Ascat
ascat_alleles     = params.ascat_alleles          ? Channel.fromPath(params.ascat_alleles).collect()          : Channel.empty()
ascat_loci       = params.ascat_loci            ? Channel.fromPath(params.ascat_loci).collect()            : Channel.empty()
ascat_gc         = params.ascat_gc              ? Channel.fromPath(params.ascat_gc).collect()              : Channel.empty()
ascat_rt         = params.ascat_rt              ? Channel.fromPath(params.ascat_rt).collect()              : Channel.empty()

// Facets
facets_annotation_bed   = params.facets_annotation_bed          ? Channel.fromPath(params.facets_annotation_bed).collect() : Channel.empty()

// FragCounter
gcmapdir_frag      = params.gcmapdir_frag      ? Channel.fromPath(params.gcmapdir_frag).collect()     : Channel.empty()   // This is the GC/Mappability directory for fragCounter. (Must contain gc* & map* .rds files)

// HetPileups
hapmap_sites       = params.hapmap_sites       ? Channel.fromPath(params.hapmap_sites).collect()      : Channel.empty()

// Dryclean
pon_dryclean      = params.pon_dryclean      ? Channel.fromPath(params.pon_dryclean).collect()     : Channel.empty()   // This is the path to the PON for Dryclean.
blacklist_path_dryclean      = params.blacklist_path_dryclean      ? Channel.fromPath(params.blacklist_path_dryclean).collect()     : Channel.empty()   // This is the path to the blacklist for Dryclean (optional).
germline_file_dryclean      = params.germline_file_dryclean      ? Channel.fromPath(params.germline_file_dryclean).collect()     : Channel.empty()   // This is the path to the germline mask for dryclean (optional).

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
//Ascat
ascat_genome       = params.ascat_genome       ?: Channel.empty()

//Facets
facets_genome      = params.facets_genome      ?: Channel.empty()

// Hetpileups
filter_hets         = params.filter_hets       ?: Channel.empty()
max_depth_hets           = params.max_depth_hets         ?: Channel.empty()

// fragCounter
windowsize_frag    = params.windowsize_frag    ?: Channel.empty()
minmapq_frag       = params.minmapq_frag       ?: Channel.empty()
midpoint_frag      = params.midpoint_frag      ?: Channel.empty()
paired_frag        = params.paired_frag        ?: Channel.empty()
exome_frag         = params.exome_frag         ?: Channel.empty()

// dryclean
center_dryclean        = params.center_dryclean        ?: Channel.empty()
cbs_dryclean           = params.cbs_dryclean           ?: Channel.empty()
cnsignif_dryclean      = params.cnsignif_dryclean      ?: Channel.empty()
wholeGenome_dryclean   = params.wholeGenome_dryclean   ?: Channel.empty()
blacklist_dryclean     = params.blacklist_dryclean     ?: Channel.empty()
germline_filter_dryclean= params.germline_filter_dryclean ?: Channel.empty()
field_dryclean         = params.field_dryclean         ?: Channel.empty()
build_dryclean         = params.build_dryclean         ?: Channel.empty()

/*
    IMPORT SUBWORKFLOWS
*/
include { INPUT_PREP } from '../subworkflows/input_prep.nf'
include { ZERO_SHOT_CNV_CALL } from '../subworkflows/zeroshot_cnv_calling.nf'
include { BAM_HETPILEUPS } from '../subworkflows/bam_hetpileups/main.nf'
include { COVERAGE_JABBA } from '../subworkflows/coverage_jabba.nf'

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
    ZERO_SHOT_CNV_CALL(
        ch_samples, 
        fasta, 
        ascat_alleles, 
        ascat_loci, 
        ascat_gc, 
        ascat_rt, 
        target_regs,
        dbsnp,
        dbsnp_tbi,
        facets_annotation_bed
    )
    ch_versions    = ch_versions.mix(ZERO_SHOT_CNV_CALL.out.versions)

    //
    // SUBWORKFLOW: Generate GC- and mappability-corrected denoised genomic coverage files for JAbBA
    //
    COVERAGE_JABBA(
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
    BAM_HETPILEUPS(
        ch_samples, 
        filter_hets, 
        max_depth_hets,
        hapmap_sites
    )
    ch_versions = ch_versions.mix(BAM_HETPILEUPS.out.versions)
    sites_from_het_pileups_wgs = Channel.empty().mix(BAM_HETPILEUPS.out.het_pileups_wgs)
}