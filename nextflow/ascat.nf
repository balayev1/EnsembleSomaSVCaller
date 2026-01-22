#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Identify SNPs from BAM files using ascatSnpPanelGenerator.pl per sample.
 */
 process GENERATE_PER_SAMPLE_SNPS {
    tag "$subjectID"

    input:
    tuple val(subjectID), path(tumor), path(normal), val(gender)
    path fasta

    output:
    path "results/ascat_snps/${subjectID}/${subjectID}-hets.tsv.0", emit: het_file

    script:
    """
    mkdir -p results/ascat_snps/${subjectID}
    ascatSnpPanelGenerator.pl -ref ${fasta[0]} -hf ${normal[0]} > results/ascat_snps/${subjectID}/${subjectID}-hets.tsv.0
    """
}

/*
 * Merge SNPs from all sample files using ascatSnpPanelMerge.pl and GC-correct using ascatSnpPanelGcCorrections.pl
 */
process MERGE_AND_GC_CORRECT_SNPS {
    label 'process_high'

    input:
    path all_hets
    path fasta

    output:
    path "SnpPositions.tsv", emit: master_list
    path "SnpGcCorrections.tsv", emit: gc_file

    script:
    """
    # 1. Merge all .tsv.0 files into one master snp list
    ascatSnpPanelMerge.pl ${fasta[0]} $all_hets > SnpPositions.tsv

    # 2. Perform GC Correction for the master list
    ascatSnpPanelGcCorrections.pl ${fasta[0]} SnpPositions.tsv > SnpGcCorrections.tsv
    """
}

/*
 * Run ascat.pl for each sample using the GC-corrected SNPs
 */
process RUN_ASCAT {
    tag "$subjectID"

    input:
    tuple val(subjectID), path(tumor), path(normal), val(gender)
    path fasta
    path snp_gc

    output:
    tuple val(subjectID), path("results/ascat/${subjectID}/*.samplestatistics.txt"), path("results/ascat/${subjectID}/*.copynumber.txt.gz"), emit: ascat_results
    path "results/ascat/${subjectID}/*", emit: all_files

    script:
    """
    mkdir -p results/ascat/${subjectID}
    ascat.pl \\
        -o results/ascat/${subjectID} \\
        -t ${tumor[0]} \\
        -n ${normal[0]} \\
        -r ${fasta[0]} \\
        -sg $snp_gc \\
        -pr WGS \\
        -g $gender \\
        -c ${task.cpus} \\
        -rs ${params.species} \\
        -ra ${params.assembly}
    """
}