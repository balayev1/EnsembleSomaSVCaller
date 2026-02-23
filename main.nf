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

params.fasta                = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai            = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.target_regs          = WorkflowMain.getGenomeAttribute(params, 'target_regs')
params.dbsnp                = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.hapmap_sites         = WorkflowMain.getGenomeAttribute(params, 'hapmap_sites')
params.ascat_alleles        = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome         = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci           = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_gc             = WorkflowMain.getGenomeAttribute(params, 'ascat_gc')
params.ascat_rt             = WorkflowMain.getGenomeAttribute(params, 'ascat_rt')
params.sequenza_genome      = WorkflowMain.getGenomeAttribute(params, 'sequenza_genome')
params.facets_genome        = WorkflowMain.getGenomeAttribute(params, 'facets_genome')
params.delly_blacklist      = WorkflowMain.getGenomeAttribute(params, 'delly_blacklist')
params.gridss_blacklist     = WorkflowMain.getGenomeAttribute(params, 'gridss_blacklist')
params.gridss_pon           = WorkflowMain.getGenomeAttribute(params, 'gridss_pon')
params.indel_mask           = WorkflowMain.getGenomeAttribute(params, 'indel_mask')
params.germ_sv_db           = WorkflowMain.getGenomeAttribute(params, 'germ_sv_db')
params.simple_seq_db        = WorkflowMain.getGenomeAttribute(params, 'simple_seq_db')
params.brass_genome_cache   = WorkflowMain.getGenomeAttribute(params, 'brass_genome_cache')
params.brass_depth          = WorkflowMain.getGenomeAttribute(params, 'brass_depth')
params.brass_viral          = WorkflowMain.getGenomeAttribute(params, 'brass_viral')
params.brass_microbe        = WorkflowMain.getGenomeAttribute(params, 'brass_microbe')
params.brass_gcbins         = WorkflowMain.getGenomeAttribute(params, 'brass_gcbins')
params.brass_cytoband       = WorkflowMain.getGenomeAttribute(params, 'brass_cytoband')
params.brass_centtel        = WorkflowMain.getGenomeAttribute(params, 'brass_centtel')
params.brass_np             = WorkflowMain.getGenomeAttribute(params, 'brass_np')
params.brass_protocol       = WorkflowMain.getGenomeAttribute(params, 'brass_protocol')
params.brass_sp             = WorkflowMain.getGenomeAttribute(params, 'brass_sp')
params.brass_genome         = WorkflowMain.getGenomeAttribute(params, 'brass_genome')
params.gcmapdir_frag        = WorkflowMain.getGenomeAttribute(params, 'gcmapdir_frag')
params.build_dryclean       = WorkflowMain.getGenomeAttribute(params, 'build_dryclean')
params.pon_dryclean         = WorkflowMain.getGenomeAttribute(params, 'pon_dryclean')

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