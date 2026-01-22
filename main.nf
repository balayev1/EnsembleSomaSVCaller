#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { 
    GENERATE_PER_SAMPLE_SNPS; 
    MERGE_AND_GC_CORRECT_SNPS; 
    RUN_ASCAT 
} from './ascat.nf'

workflow {
    // 1. Setup Input Channel
    ch_samples = Channel.fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> tuple(row.subjectID, file("${row.tumorBAM}*"), file("${row.normalBAM}*"), row.gender) }
    
    ref_bundle = [ file(params.fasta), file("${params.fasta}.fai") ]

    // 2. Parallel SNP generation
    individual_hets = GENERATE_PER_SAMPLE_SNPS(ch_samples, ref_bundle)

    // 3. Merge and GC-correct SNPs
    // .collect() gathers all individual_hets files into a single list
    custom_panel = MERGE_AND_GC_CORRECT_SNPS(individual_hets.het_file.collect(), ref_bundle)

    // 4. Final Step: Run ASCAT for each sample
    RUN_ASCAT(
        ch_samples, 
        ref_bundle, 
        custom_panel.gc_file)
}