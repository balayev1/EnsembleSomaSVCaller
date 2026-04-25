//
// CBS
//
//

include { CBS } from '../../modules/local/cbs/main.nf'


// Define the main workflow process
workflow COV_CBS {
    // Define the input parameters for the main workflow
    take:
    cov         // channel: [mandatory] [ meta, tumor_cov, normal_cov ]
    cnsignif
    field
    name

    main:
    cbs_cov_rds             = Channel.empty()
    cbs_seg_rds             = Channel.empty()
    cbs_nseg_rds            = Channel.empty()
    versions                = Channel.empty()

    CBS(cov, cnsignif, field, name)

    cbs_cov_rds              = CBS.out.cbs_cov_rds
    cbs_seg_rds              = CBS.out.cbs_seg_rds
    cbs_nseg_rds             = CBS.out.cbs_nseg_rds

    versions = versions.mix(CBS.out.versions)

    emit:
    cbs_cov_rds
    cbs_seg_rds
    cbs_nseg_rds

    versions
}