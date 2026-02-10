#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    VALIDATE INPUTS
*/

// Validate input parameters
WorkflowSomasvcaller.initialise(params, log)

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

/*
    IMPORT SUBWORKFLOWS
*/
include { INPUT_PREP } from './subworkflows/input_prep'
include { ZERO_SHOT_CNV_CALL } from './subworkflows/zeroshot_cnv_calling.nf'

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
}