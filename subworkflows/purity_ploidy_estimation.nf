//
// Calling CNVs using a zero-shot approach with ASCAT, SEQUENZA, and FACETS.
// External ACEseq outputs are validated through the manifest passed at launch time,
// then merged into a per-sample purity/ploidy consensus.
//

include { ASCAT }                  from '../modules/nf-core/ascat/main.nf'
include { FACETS }                 from '../modules/local/facets/main.nf'
include { PURITY_PLOIDY_MERGE }    from '../modules/local/purity_ploidy_merge/main.nf'
include { SEQUENZAUTILS_GCWIGGLE } from '../modules/local/sequenza/main.nf'
include { SEQUENZA_PREP }          from '../modules/local/sequenza/main.nf'
include { SEQUENZA_MERGE }         from '../modules/local/sequenza/main.nf'
include { SEQUENZA_RUN }           from '../modules/local/sequenza/main.nf'

def readTsvRows(path) {
    def lines = java.nio.file.Files.readAllLines(path).findAll { it?.trim() }
    if (lines.size() < 2) {
        return []
    }

    def header = lines[0].split('\t', -1)
    lines.tail().collect { line ->
        def values = line.split('\t', -1)
        def row = [:]
        header.eachWithIndex { key, idx ->
            row[key] = idx < values.size() ? values[idx] : ''
        }
        row
    }
}

def resolveManifestPath(value, baseDir) {
    def raw = value?.toString()?.trim()
    if (!raw) {
        return null
    }

    def candidate = java.nio.file.Paths.get(raw)
    def resolved = candidate.isAbsolute() ? candidate : baseDir.resolve(candidate).normalize()
    file(resolved.toString(), checkIfExists: true)
}

workflow PURITY_PLOIDY_ESTIMATION {
    take:
        ch_samples
        fasta
        fai
        ascat_alleles
        ascat_loci
        ascat_gc
        ascat_rt
        target_regs
        dbsnp
        dbsnp_tbi
        facets_annotation
        aceseq_manifest_tsv

    main:
        versions = Channel.empty()

        ASCAT(
            ch_samples,
            ascat_alleles,
            ascat_loci,
            target_regs,
            fasta,
            ascat_gc,
            ascat_rt
        )
        ascat_purityploidy = ASCAT.out.purityploidy
        versions = versions.mix(ASCAT.out.versions)

        ch_gc_wiggle_input = Channel.value([id: 'reference_genome']).combine(fasta)
        SEQUENZAUTILS_GCWIGGLE(
            ch_gc_wiggle_input,
            params.windowsize_sequtils
        )
        ch_wiggle_file = SEQUENZAUTILS_GCWIGGLE.out.wig.map { meta, wig -> wig }.first()
        versions = versions.mix(SEQUENZAUTILS_GCWIGGLE.out.versions)

        ch_chromosomes = fai
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .filter { it ==~ /^(chr)?([0-9]+|[XY])$/ }

        ch_sequenza_prep_input = ch_samples.combine(ch_chromosomes)
        SEQUENZA_PREP(
            ch_sequenza_prep_input,
            fasta,
            ch_wiggle_file,
            params.hom_sequtils,
            params.het_sequtils,
            params.het_f_sequtils,
            params.qlimit_sequtils,
            params.qformat_sequtils,
            params.rd_thr_sequtils,
            params.seqz_bin_size_sequtils
        )
        versions = versions.mix(SEQUENZA_PREP.out.versions)

        ch_to_merge = SEQUENZA_PREP.out.seqz_file.groupTuple(by: 0)
        SEQUENZA_MERGE(ch_to_merge)
        versions = versions.mix(SEQUENZA_MERGE.out.versions)

        SEQUENZA_RUN(
            SEQUENZA_MERGE.out.merged_seqz,
            params.sequenza_genome,
            params.sequenza_rd_thr_tumor,
            params.sequenza_rd_thr_normal,
            params.sequenza_purity_range,
            params.sequenza_ploidy_range
        )
        versions = versions.mix(SEQUENZA_RUN.out.versions)

        FACETS(
            ch_samples,
            dbsnp,
            dbsnp_tbi,
            target_regs,
            facets_annotation,
            params.facets_mapq,
            params.facets_baq,
            params.facets_mindepth,
            params.facets_maxdepth,
            params.facets_cval_pre,
            params.facets_cval_proc,
            params.facets_nbhd_snp,
            params.facets_genome,
            params.facets_count_orphans,
            params.facets_unmatched,
            params.facets_no_cov_plot
        )
        versions = versions.mix(FACETS.out.versions)

        if (params.skip_aceseq) {
            purity_ploidy_merge_input = ASCAT.out.purityploidy
                .map { meta, ascat_file -> [meta.id, meta, ascat_file] }
                .join(FACETS.out.vcf.map { meta, facets_vcf -> [meta.id, facets_vcf] })
                .join(SEQUENZA_RUN.out.purity_ploidy_est.map { meta, sequenza_file -> [meta.id, sequenza_file] })
                .map { sample_id, meta, ascat_file, facets_vcf, sequenza_file ->
                    [meta, ascat_file, facets_vcf, sequenza_file, '']
                }
        } else {
            aceseq_ploidy_purity = aceseq_manifest_tsv
                .flatMap { manifest_path ->
                    def manifest_base_dir = manifest_path.parent
                    readTsvRows(manifest_path).collect { row ->
                        [
                            row.sample_id.toString(),
                            resolveManifestPath(row.ploidy_purity_2d.toString(), manifest_base_dir)
                        ]
                    }
                }

            purity_ploidy_merge_input = ASCAT.out.purityploidy
                .map { meta, ascat_file -> [meta.id, meta, ascat_file] }
                .join(FACETS.out.vcf.map { meta, facets_vcf -> [meta.id, facets_vcf] })
                .join(SEQUENZA_RUN.out.purity_ploidy_est.map { meta, sequenza_file -> [meta.id, sequenza_file] })
                .join(aceseq_ploidy_purity)
                .map { sample_id, meta, ascat_file, facets_vcf, sequenza_file, aceseq_file ->
                    [meta, ascat_file, facets_vcf, sequenza_file, aceseq_file]
                }
        }

        PURITY_PLOIDY_MERGE(purity_ploidy_merge_input)
        purity_ploidy_candidates = PURITY_PLOIDY_MERGE.out.candidates
        purity_ploidy_consensus  = PURITY_PLOIDY_MERGE.out.consensus
        versions                 = versions.mix(PURITY_PLOIDY_MERGE.out.versions)

    emit:
        ascat_purityploidy
        purity_ploidy_candidates
        purity_ploidy_consensus
        versions
}
