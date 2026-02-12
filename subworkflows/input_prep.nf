#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Input Preparation Subworkflow
//


workflow INPUT_PREP {
    take:
    csv_file // file: /path/to/samplesheet.csv

    main:
    ch_samples = Channel.fromPath(csv_file)
        .splitCsv( header:true, sep:',' )
        .map { row ->
            def meta = [
                id: row.sample,
                gender: row.sex,
                sv_type: row.sv ?: []
            ]

            // Validate If Files Exist
            if (!file(row.control).exists()) error "Control BAM not found: ${row.control}"
            if (!file(row.control_index).exists()) error "Control Index not found: ${row.control_index}"
            if (!file(row.tumor).exists()) error "Tumor BAM not found: ${row.tumor}"
            if (!file(row.tumor_index).exists()) error "Tumor Index not found: ${row.tumor_index}"

            // Construct [ val(meta), [control], [control_index],[ tumor], [tumor_index]]
            return [
                meta,
                file(row.control),       // Input Normal
                file(row.control_index), // Index Normal
                file(row.tumor),         // Input Tumor
                file(row.tumor_index)    // Index Tumor
            ]
        }

    emit:
    ch_samples // channel: [ val(meta), [control], [control_index],[ tumor], [tumor_index]]
}