
//
// GRIDSS SV CALLING
//
//
//

include { GRIDSS_GRIDSS   } from '../../modules/nf-core/gridss/gridss/main.nf'
include { GRIDSS_SOMATIC_FILTER  } from '../../modules/nf-core/gridss/somaticfilter/main.nf'

workflow GRIDSS_SV_CALLING {
    take:
    input                                              // channel: [mandatory] [ meta, normalbam, normalbai, tumorbam, tumorbai ]
    bwa_index                                          // channel: [mandatory] bwa index path
    fasta                                              // channel: [mandatory] [ path(fasta) ]
    fasta_fai                                          // channel: [mandatory] [ path(fasta_fai) ]
    blacklist_gridss                                   // channel: [optional] gridss blacklist bed file based on genome

    main:
    versions               = Channel.empty()
    vcf                    = Channel.empty()
    vcf_index              = Channel.empty()
    assembly_bam           = Channel.empty()

    GRIDSS_GRIDSS(input, fasta, fasta_fai, bwa_index, blacklist_gridss)

    vcf                    = GRIDSS_GRIDSS.out.filtered_vcf
    vcf_index              = GRIDSS_GRIDSS.out.filtered_vcf_index


    versions = versions.mix(GRIDSS_GRIDSS.out.versions)

    emit:
    vcf     // channel: [mandatory] [ meta, pass filltered vcf, pass filtered vcf tbi ]
    vcf_index


    versions

}



workflow GRIDSS_SOMATIC_FILTER_STEP {
    take:
    vcf
    pondir_gridss

    main:
    versions                = Channel.empty()
    somatic_all             = Channel.empty()
    somatic_high_confidence = Channel.empty()

    GRIDSS_SOMATIC_FILTER(vcf, pondir_gridss)

    somatic_high_confidence = GRIDSS_SOMATIC_FILTER.out.somatic_high_vcf
    somatic_all             = GRIDSS_SOMATIC_FILTER.out.somatic_all_vcf

    versions                = GRIDSS_SOMATIC_FILTER.out.versions

    emit:
    somatic_high_confidence // channel: [mandatory] [ meta, high confidence somatic vcf ]
    somatic_all
    versions
}