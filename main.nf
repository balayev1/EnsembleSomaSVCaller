#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EnsembleSomaSVCaller
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/balayev1/EnsembleSomaSVCaller

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai        = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.target_regs      = WorkflowMain.getGenomeAttribute(params, 'target_regs')
params.dbsnp            = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi        = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.hapmap_sites     = WorkflowMain.getGenomeAttribute(params, 'hapmap_sites')
params.ascat_alleles    = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome     = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci       = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_gc         = WorkflowMain.getGenomeAttribute(params, 'ascat_gc')
params.ascat_rt         = WorkflowMain.getGenomeAttribute(params, 'ascat_rt')
params.sequenza_genome  = WorkflowMain.getGenomeAttribute(params, 'sequenza_genome')
params.facets_genome    = WorkflowMain.getGenomeAttribute(params, 'facets_genome')
params.delly_blacklist  = WorkflowMain.getGenomeAttribute(params, 'delly_blacklist')
params.gcmapdir_frag    = WorkflowMain.getGenomeAttribute(params, 'gcmapdir_frag')
params.build_dryclean   = WorkflowMain.getGenomeAttribute(params, 'build_dryclean')
params.pon_dryclean     = WorkflowMain.getGenomeAttribute(params, 'pon_dryclean')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SOMASV_CALLER } from './workflows/somasvcaller.nf'

//
// WORKFLOW: Run main ENSEMBLESOMASVCALLER analysis pipeline
//
workflow NFCORE_SOMASV {
    SOMASV_CALLER ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    NFCORE_SOMASV ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/