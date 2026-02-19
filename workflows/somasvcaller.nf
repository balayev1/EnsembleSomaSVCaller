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
target_regs         = params.target_regs        ? Channel.fromPath(params.target_regs).collect()         : Channel.value([])
dbsnp               = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()              : Channel.empty()
dbsnp_tbi           = params.dbsnp_tbi          ? Channel.fromPath(params.dbsnp_tbi).collect()              : Channel.empty()

// Ascat
ascat_alleles     = params.ascat_alleles        ? Channel.fromPath(params.ascat_alleles).collect()          : Channel.empty()
ascat_loci       = params.ascat_loci            ? Channel.fromPath(params.ascat_loci).collect()            : Channel.empty()
ascat_gc         = params.ascat_gc              ? Channel.fromPath(params.ascat_gc).collect()              : Channel.empty()
ascat_rt         = params.ascat_rt              ? Channel.fromPath(params.ascat_rt).collect()              : Channel.empty()

// Facets
facets_annotation_bed   = params.facets_annotation_bed          ? Channel.fromPath(params.facets_annotation_bed).collect() : Channel.value([])

// FragCounter
gcmapdir_frag      = params.gcmapdir_frag      ? Channel.fromPath(params.gcmapdir_frag).collect()     : Channel.empty()   // This is the GC/Mappability directory for fragCounter. (Must contain gc* & map* .rds files)

// Delly
delly_blacklist    = params.delly_blacklist    ? Channel.fromPath(params.delly_blacklist).collect()     : Channel.value([])

// HetPileups
hapmap_sites       = params.hapmap_sites       ? Channel.fromPath(params.hapmap_sites).collect()      : Channel.empty()

// Dryclean
pon_dryclean      = params.pon_dryclean      ? Channel.fromPath(params.pon_dryclean).collect()     : Channel.empty()   // This is the path to the PON for Dryclean.
blacklist_path_dryclean      = params.blacklist_path_dryclean      ? Channel.fromPath(params.blacklist_path_dryclean).collect()     : Channel.empty()   // This is the path to the blacklist for Dryclean (optional).
germline_file_dryclean      = params.germline_file_dryclean      ? Channel.fromPath(params.germline_file_dryclean).collect()     : Channel.empty()   // This is the path to the germline mask for dryclean (optional).

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
//Ascat
ascat_genome       = params.ascat_genome       ?: Channel.empty()

//Sequenzautils
windowsize_sequtils = params.windowsize_sequtils ?: Channel.empty()
hom_sequtils        = params.hom_sequtils        ?: Channel.empty()
het_sequtils        = params.het_sequtils        ?: Channel.empty()
het_f_sequtils      = params.het_f_sequtils      ?: Channel.empty()
qlimit_sequtils     = params.qlimit_sequtils     ?: Channel.empty()
qformat_sequtils    = params.qformat_sequtils    ?: Channel.empty()
rd_thr_sequtils     = params.rd_thr_sequtils     ?: Channel.empty()
seqz_bin_size_sequtils = params.seqz_bin_size_sequtils ?: Channel.empty()

// Sequenza
sequenza_genome      = params.sequenza_genome      ?: Channel.empty()
sequenza_rd_thr_tumor   = params.sequenza_rd_thr_tumor   ?: Channel.empty()
sequenza_rd_thr_normal  = params.sequenza_rd_thr_normal  ?: Channel.empty()
sequenza_purity_range   = params.sequenza_purity_range   ?: Channel.empty()
sequenza_ploidy_range   = params.sequenza_ploidy_range   ?: Channel.empty()

//Facets
facets_genome      = params.facets_genome      ?: Channel.empty()

// Delly
delly_mode                  = params.delly_mode                 ?: Channel.empty()
delly_min_svqual            = params.delly_min_svqual           ?: Channel.empty()
delly_altaf                 = params.delly_altaf               ?: Channel.empty()
delly_min_svsize            = params.delly_min_svsize           ?: Channel.empty()
delly_max_svsize            = params.delly_max_svsize           ?: Channel.empty()
delly_ratio_geno            = params.delly_ratio_geno           ?: Channel.empty()
delly_pass                  = params.delly_pass                 ?: Channel.empty()
delly_tags                  = params.delly_tags                 ?: Channel.empty()
delly_cov                   = params.delly_cov                  ?: Channel.empty()
delly_ctrl_contamination    = params.delly_ctrl_contamination   ?: Channel.empty()
delly_geno_qual             = params.delly_geno_qual            ?: Channel.empty()
delly_rddel                 = params.delly_rddel                ?: Channel.empty()
delly_rddup                 = params.delly_rddup                ?: Channel.empty()

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
    // ch_samples.view()

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

    //
    // SUBWORKFLOW: Estimate genomic breakpoints for JabBa
    //
    BREAKPOINT_ESTIMATOR (
        ch_samples,
        fasta,
        fasta_fai,
        delly_blacklist,
        delly_mode,
        delly_min_svqual,
        delly_altaf,
        delly_min_svsize,
        delly_max_svsize,
        delly_ratio_geno,
        delly_pass,
        delly_tags,
        delly_cov,
        delly_ctrl_contamination,
        delly_geno_qual,
        delly_rddel,
        delly_rddup
    )
    ch_versions = ch_versions.mix(BREAKPOINT_ESTIMATOR.out.versions)
}