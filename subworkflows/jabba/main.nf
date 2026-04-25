#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// JaBbA
//

include { JABBA } from '../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_COV } from '../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_JUNCTION } from '../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_HETS } from '../../modules/local/jabba/main.nf'
include { RETIER_WHITELIST_JUNCTIONS } from '../../modules/local/jabba/main.nf'

workflow JABBA_RUN {
    take:
    junctions
    cov_rds
    het_pileups_wgs
    cbs_seg_rds
    cbs_nseg_rds
    purity_ploidy_consensus

    main:
    versions = Channel.empty()

    junctions_for_jabba = Channel.empty()
    if (params.is_retier_whitelist_junctions) {
        if (!params.whitelist_genes_jabba || params.whitelist_genes_jabba in ['NULL', 'NA']) {
            error "params.whitelist_genes_jabba must be provided when params.is_retier_whitelist_junctions is true"
        }
        RETIER_WHITELIST_JUNCTIONS(
            junctions,
            params.tfield_jabba,
            file(params.whitelist_genes_jabba, checkIfExists: true)
        )
        junctions_for_jabba = RETIER_WHITELIST_JUNCTIONS.out.retiered_junctions
    } else {
        junctions_for_jabba = junctions
    }

    COERCE_SEQNAMES_JUNCTION(junctions_for_jabba)
    coerced_junctions = COERCE_SEQNAMES_JUNCTION.out.file
        .map { meta, junction_file -> [meta.id, meta, junction_file] }

    COERCE_SEQNAMES_COV(cov_rds)
    coerced_cov_rds = COERCE_SEQNAMES_COV.out.file
        .map { meta, cov_file -> [meta.id, cov_file] }

    COERCE_SEQNAMES_HETS(het_pileups_wgs)
    coerced_hets = COERCE_SEQNAMES_HETS.out.file
        .map { meta, hets_file -> [meta.id, hets_file] }

    cbs_seg_by_id = cbs_seg_rds.map { meta, seg_file -> [meta.id, seg_file] }
    cbs_nseg_by_id = cbs_nseg_rds.map { meta, nseg_file -> [meta.id, nseg_file] }

    purity_ploidy_values = purity_ploidy_consensus.map { meta, consensus_file ->
        def lines = java.nio.file.Files.readAllLines(consensus_file).findAll { it?.trim() }
        if (lines.size() < 2) {
            error "Invalid purity/ploidy consensus file for ${meta.id}: ${consensus_file}"
        }

        def header = lines[0].split('\t', -1)
        def values = lines[1].split('\t', -1)
        def row = [:]
        header.eachWithIndex { key, idx ->
            row[key] = idx < values.size() ? values[idx] : ''
        }

        def purity = (row.status == 'consensus' && row.purity && row.purity != 'NA') ? row.purity : (params.purity_jabba ?: 'NA')
        def ploidy = (row.status == 'consensus' && row.ploidy && row.ploidy != 'NA') ? row.ploidy : (params.ploidy_jabba ?: 'NA')

        [meta.id, purity, ploidy]
    }

    jabba_inputs = coerced_junctions
        .join(coerced_cov_rds)
        .join(coerced_hets)
        .join(cbs_seg_by_id)
        .join(cbs_nseg_by_id)
        .join(purity_ploidy_values)
        .map { sample_id, meta, junction_file, cov_file, hets_file, seg_file, nseg_file, purity, ploidy ->
            [meta, junction_file, cov_file, 'NULL', hets_file, purity, ploidy, seg_file, nseg_file]
        }

    JABBA(
        jabba_inputs,
        params.blacklist_junctions_jabba ?: 'NULL',
        params.geno_jabba,
        params.indel_jabba,
        params.tfield_jabba,
        params.iter_jabba,
        params.rescue_window_jabba,
        params.rescue_all_jabba,
        params.nudgebalanced_jabba,
        params.edgenudge_jabba,
        params.strict_jabba,
        params.allin_jabba,
        params.field_jabba,
        params.maxna_jabba,
        params.blacklist_coverage_jabba ?: 'NULL',
        params.pp_method_jabba,
        params.cnsignif_jabba,
        params.slack_jabba,
        params.linear_jabba,
        params.tilim_jabba,
        params.epgap_jabba,
        params.fix_thres_jabba,
        params.lp_jabba,
        params.ism_jabba,
        params.filter_loose_jabba,
        params.gurobi_jabba,
        params.gurobi_license_path,
        params.verbose_jabba
    )

    versions = versions.mix(JABBA.out.versions)

    emit:
    jabba_rds     = JABBA.out.jabba_rds
    jabba_gg      = JABBA.out.jabba_gg
    jabba_vcf     = JABBA.out.jabba_vcf
    jabba_raw_rds = JABBA.out.jabba_raw_rds
    opti          = JABBA.out.opti
    jabba_seg     = JABBA.out.jabba_seg
    karyograph    = JABBA.out.karyograph
    versions
}
