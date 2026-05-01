#!/usr/bin/env nextflow

import WorkflowMain

nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

def initGenomeParams() {
    [
        'fasta',
        'fasta_fai',
        'target_regs',
        'dbsnp',
        'dbsnp_tbi',
        'chrom_sizes',
        'bwa_index',
        'hapmap_sites',
        'ascat_alleles',
        'ascat_genome',
        'ascat_loci',
        'ascat_gc',
        'ascat_rt',
        'sequenza_genome',
        'facets_genome',
        'delly_blacklist',
        'gridss_blacklist',
        'gridss_pon',
        'indel_mask',
        'germ_sv_db',
        'simple_seq_db',
        'simple_seq_db_slop',
        'brass_cache_dir',
        'brass_depth',
        'brass_viral',
        'brass_microbe',
        'brass_gcbins',
        'brass_cytoband',
        'brass_centtel',
        'brass_np',
        'brass_np_tbi',
        'brass_protocol',
        'brass_sp',
        'brass_genome',
        'gcmapdir_frag',
        'build_dryclean',
        'pon_dryclean'
    ].each { attribute ->
        params[attribute] = WorkflowMain.getGenomeAttribute(params, attribute)
    }

    params.blacklist_coverage_jabba = WorkflowMain.getGenomeAttribute(params, 'blacklist_coverage_jabba') ?: params.blacklist_coverage_jabba
    params.whitelist_genes_jabba = WorkflowMain.getGenomeAttribute(params, 'whitelist_genes_jabba') ?: params.whitelist_genes_jabba
}

initGenomeParams()

include { SOMASV_CALLER } from './workflows/somasvcaller.nf'

workflow NFCORE_SOMASV {
    SOMASV_CALLER()
}

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)
    NFCORE_SOMASV()
}
