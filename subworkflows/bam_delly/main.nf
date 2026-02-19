//
// BAM DELLY
//

include { DELLY_CALL } from '../../modules/nf-core/delly/call/main.nf'
include { DELLY_FILTER } from '../../modules/local/delly/filter/main.nf'
include { CREATE_DELLY_SAMPLES_TSV } from '../../modules/local/delly/create_sample_tsv/main.nf'

workflow BAM_DELLY {
    // defining inputs
    take:
    input                                               // channel: [ val(meta), [bams], [bais], exclude_bed ]
    fasta                                               // channel: [ val(meta2), path(fasta) ]
    fai                                                 // channel: [ val(meta3), path(fai) ]
    mode                                                // channel: val(mode)
    min_qual                                            // channel: val(min_qual)
    altaf                                               // channel: val(altaf)
    min_sv_size                                         // channel: val(min_sv_size)
    max_sv_size                                         // channel: val(max_sv_size)
    ratio_geno                                          // channel: val(ratio_geno)
    pass                                                // channel: val(pass)
    tags                                                // channel: val(tags)
    cov                                                 // channel: val(cov)
    contamination_in_control                            // channel: val(contamination_in_control)
    geno_qual                                           // channel: val(geno_qual)
    rd_del                                              // channel: val(rd_del)
    rd_dup                                              // channel: val(rd_dup)

    //Creating empty channels for output
    main:
    versions                = Channel.empty()
    delly_call_bcf          = Channel.empty()
    delly_call_bcf_csi      = Channel.empty()
    delly_filter_bcf        = Channel.empty()
    delly_samples_tsv       = Channel.empty()

    delly_call_input = input.map { meta, bams, bais, exclude ->
        def tumor_id   = bams[1].getSimpleName() 
        def control_id = bams[0].getSimpleName()
        
        def new_meta = meta.clone()
        new_meta.tumor_id   = tumor_id
        new_meta.control_id = control_id

        [ new_meta, bams, bais, [], [], exclude ]
    }

    DELLY_CALL(
        delly_call_input, 
        fasta, 
        fai
    )
    versions            =   versions.mix(DELLY_CALL.out.versions)
    delly_call_bcf      =   DELLY_CALL.out.bcf
    delly_call_bcf_csi  =   DELLY_CALL.out.csi

    meta_for_tsv = DELLY_CALL.out.bcf.map { meta, bcf -> [meta] }
    CREATE_DELLY_SAMPLES_TSV(meta_for_tsv)
    delly_samples_tsv = CREATE_DELLY_SAMPLES_TSV.out.tsv

    filter_input = delly_call_bcf.join(delly_call_bcf_csi).join(delly_samples_tsv)
    DELLY_FILTER(
        filter_input, 
        mode,
        min_qual,
        altaf,
        min_sv_size,
        max_sv_size,
        ratio_geno,
        pass,
        tags,
        cov,
        contamination_in_control,
        geno_qual,
        rd_del,
        rd_dup
    )
    delly_filter_bcf        = DELLY_FILTER.out.bcf
    versions                = versions.mix(DELLY_FILTER.out.versions)

    //
    emit:
    
    delly_call_bcf
    delly_call_bcf_csi
    delly_filter_bcf

    versions  
}