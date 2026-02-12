//
// BAM HETPILEUPS
//

include { HETPILEUPS } from '../../../modules/nf-core/hetpileups/main.nf'

workflow BAM_HETPILEUPS {
    // defining inputs
    take:
    input                                             // required: Format should be [meta, tumor bam, normal bam] : BAM only!!!
	filter                                                               
	max_depth
	hapmap_sites

    //Creating empty channels for output
    main:
    versions            = Channel.empty()
    het_pileups_wgs     = Channel.empty()

    //Remapping the input based on Hetpileups module
    bam_hets = input.map { meta, control, control_index, tumor, tumor_index ->
        [meta, tumor, tumor_index, control, control_index]
    }

    HETPILEUPS(bam_hets, filter, max_depth, hapmap_sites)

    // initializing outputs from Hetpileups
    versions          = HETPILEUPS.out.versions
    het_pileups_wgs   = HETPILEUPS.out.het_pileups_wgs

    //
    emit:
    het_pileups_wgs

    versions
}