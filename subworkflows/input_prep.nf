#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Input Preparation Subworkflow
//

def resolveInputPath(value, baseDir) {
    def raw = value?.toString()?.trim()
    if (!raw) {
        return null
    }

    def candidate = java.nio.file.Paths.get(raw)
    def resolved = candidate.isAbsolute() ? candidate : baseDir.resolve(candidate).normalize()
    file(resolved.toString())
}


workflow INPUT_PREP {
    take:
    csv_file // file: /path/to/samplesheet.csv

    main:
    def csv_base_dir = csv_file.parent
    ch_samples = Channel.fromPath(csv_file)
        .splitCsv( header:true, sep:',' )
        .map { row ->
            def meta = [
                id: row.sample,
                gender: row.sex,
                sv_type: row.sv ?: []
            ]

            def control = resolveInputPath(row.control, csv_base_dir)
            def control_index = resolveInputPath(row.control_index, csv_base_dir)
            def tumor = resolveInputPath(row.tumor, csv_base_dir)
            def tumor_index = resolveInputPath(row.tumor_index, csv_base_dir)

            // Validate If Files Exist
            if (!control.exists()) error "Control BAM not found: ${row.control}"
            if (!control_index.exists()) error "Control Index not found: ${row.control_index}"
            if (!tumor.exists()) error "Tumor BAM not found: ${row.tumor}"
            if (!tumor_index.exists()) error "Tumor Index not found: ${row.tumor_index}"

            // Construct [ val(meta), [normal],[normal.bai],[ tumor], [tumor.bai]]
            return [
                meta,
                control,       // Input Normal
                control_index, // Index Normal
                tumor,         // Input Tumor
                tumor_index    // Index Tumor
            ]
        }

    emit:
    ch_samples // channel: [ val(meta), [normal],[normal.bai],[ tumor], [tumor.bai]]
}
