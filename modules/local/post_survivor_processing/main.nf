process POST_SURVIVOR_PROCESSING {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bcftools=1.20" : null)
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: vcf.name.replaceFirst(/\.vcf(\.(gz|bgz))?$/, '') + '.postprocessed'
    def output_vcf = "${prefix}.vcf"
    def normalizedSex = ((meta.gender ?: meta.sex ?: '') as String).trim().toLowerCase()
    def isFemale = ['female', 'f', 'xx'].contains(normalizedSex)

    if ("${vcf}" == "${output_vcf}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    awk -v is_female=${isFemale ? 1 : 0} '
    BEGIN {
        FS = OFS = "\\t"
        tier_header = "##INFO=<ID=TIER,Number=1,Type=Integer,Description=\\"Support-derived tier: 1 if INFO\\/SUPP>1, 2 if INFO\\/SUPP=1\\">"
        mate_header = "##INFO=<ID=MATEID,Number=.,Type=String,Description=\\"ID of mate breakends\\">"
        event_header = "##INFO=<ID=EVENT,Number=1,Type=String,Description=\\"ID of event associated to breakend\\">"
        gt_header = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">"
    }
    function add_field(info, field) {
        return info == "" ? field : info ";" field
    }
    function clamp_pos(pos) {
        return pos < 1 ? 1 : int(pos)
    }
    function ct_to_strands(ct_value) {
        if (ct_value == "5to5") return "++"
        if (ct_value == "5to3") return "+-"
        if (ct_value == "3to5") return "-+"
        if (ct_value == "3to3") return "--"
        return ""
    }
    function bnd_alt(local_strand, mate_strand, mate_chrom, mate_pos, ref, token) {
        token = mate_strand == "+" ? "[" mate_chrom ":" mate_pos "[" : "]" mate_chrom ":" mate_pos "]"
        return local_strand == "-" ? ref token : token ref
    }
    function emit_bnd_record(chrom, pos, id, mateid, eventid, mate_chrom, mate_pos, local_strand, mate_strand, qual, filter, info_base, gt, tier, info, alt) {
        info = info_base
        info = add_field(info, "SVTYPE=BND")
        info = add_field(info, "MATEID=" mateid)
        info = add_field(info, "EVENT=" eventid)
        info = add_field(info, "TIER=" tier)
        alt = bnd_alt(local_strand, mate_strand, mate_chrom, mate_pos, "N")
        print chrom, clamp_pos(pos), id, "N", alt, qual, filter, info, "GT", gt
    }
    function emit_bnd_pair(id_prefix, eventid, pair_tag, chrom1, pos1, strand1, chrom2, pos2, strand2, qual, filter, info_base, gt, tier, id1, id2) {
        id1 = id_prefix "_" pair_tag "_1"
        id2 = id_prefix "_" pair_tag "_2"
        emit_bnd_record(chrom1, pos1, id1, id2, eventid, chrom2, pos2, strand1, strand2, qual, filter, info_base, gt, tier)
        emit_bnd_record(chrom2, pos2, id2, id1, eventid, chrom1, pos1, strand2, strand1, qual, filter, info_base, gt, tier)
    }
    /^##INFO=<ID=TIER,/ {
        has_tier_header = 1
        print
        next
    }
    /^##INFO=<ID=MATEID,/ {
        has_mate_header = 1
        print
        next
    }
    /^##INFO=<ID=EVENT,/ {
        has_event_header = 1
        print
        next
    }
    /^##INFO=<ID=CT,/ {
        next
    }
    /^##FORMAT=<ID=GT,/ {
        has_gt_header = 1
        print
        next
    }
    /^##FORMAT=/ {
        next
    }
    /^##contig=<ID=chrY([,>])|^##contig=<ID=Y([,>])/ {
        if (is_female) {
            next
        }
    }
    /^#/ {
        if (\$0 ~ /^#CHROM[[:space:]]/) {
            if (!has_tier_header) {
                print tier_header
            }
            if (!has_mate_header) {
                print mate_header
            }
            if (!has_event_header) {
                print event_header
            }
            if (!has_gt_header) {
                print gt_header
            }
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, "FORMAT", "sample"
        } else {
            print
        }
        next
    }
    {
        base_id = (\$3 == "" || \$3 == ".") ? "event" NR : \$3
        chrom1 = \$1
        start = \$2 + 0
        raw_alt = \$5
        qual = \$6
        filter = \$7

        svtype = ""
        chr2 = ""
        end = -1
        strands = ""
        ct = ""
        supp = ""
        gt_value = "./."
        fallback_gt = "./."
        info_base = ""

        info_fields = split(\$8, info_tokens, ";")
        for (i = 1; i <= info_fields; i++) {
            field = info_tokens[i]
            if (field == "" || field == ".") {
                continue
            }

            split(field, kv, "=")
            key = kv[1]
            value = substr(field, length(key) + 2)

            if (key == "SVTYPE") {
                svtype = value
                continue
            }
            if (key == "CHR2") {
                chr2 = value
                continue
            }
            if (key == "END") {
                end = value + 0
                continue
            }
            if (key == "STRANDS") {
                strands = value
                continue
            }
            if (key == "CT") {
                ct = value
                continue
            }
            if (key == "SUPP") {
                supp = value
            }
            if (key == "MATEID" || key == "EVENT" || key == "TIER") {
                continue
            }

            info_base = add_field(info_base, field)
        }

        if (strands == "" && ct != "") {
            strands = ct_to_strands(ct)
        }

        if (is_female) {
            chrom_is_y = (chrom1 == "chrY" || chrom1 == "Y")
            chr2_is_y = (chr2 == "chrY" || chr2 == "Y")
            alt_has_y = (raw_alt ~ /(\\[|\\])chrY:/ || raw_alt ~ /(\\[|\\])Y:/)

            if (chrom_is_y || chr2_is_y || alt_has_y) {
                next
            }
        }

        if (supp == "" || supp !~ /^[0-9]+\$/) {
            printf("ERROR: Missing or invalid SUPP in INFO for %s:%s\\n", chrom1, start) > "/dev/stderr"
            exit 1
        }

        gt_idx = 0
        if (NF >= 10) {
            format_fields = split(\$9, format_tokens, ":")
            for (i = 1; i <= format_fields; i++) {
                if (format_tokens[i] == "GT") {
                    gt_idx = i
                    break
                }
            }
        }

        if (gt_idx > 0) {
            for (j = NF; j >= 10; j--) {
                sample_fields = split(\$j, sample_tokens, ":")
                if (gt_idx <= sample_fields) {
                    candidate_gt = sample_tokens[gt_idx]
                    if (fallback_gt == "./." && candidate_gt != "" && candidate_gt != ".") {
                        fallback_gt = candidate_gt
                    }
                    if (candidate_gt != "" && candidate_gt != "." && candidate_gt != "./." && candidate_gt != ".|.") {
                        gt_value = candidate_gt
                        break
                    }
                }
            }
            if (gt_value == "./." && fallback_gt != "./.") {
                gt_value = fallback_gt
            }
        }

        tier = supp > 1 ? 1 : 2

        if (svtype == "DEL") {
            if (end < 0) {
                printf("ERROR: Missing END for DEL at %s:%s\\n", chrom1, start) > "/dev/stderr"
                exit 1
            }
            emit_bnd_pair(base_id, base_id, "A", chrom1, start - 1, "-", chrom1, end + 1, "+", qual, filter, info_base, gt_value, tier)
        } else if (svtype == "DUP") {
            if (end < 0) {
                printf("ERROR: Missing END for DUP at %s:%s\\n", chrom1, start) > "/dev/stderr"
                exit 1
            }
            emit_bnd_pair(base_id, base_id, "A", chrom1, end, "-", chrom1, start, "+", qual, filter, info_base, gt_value, tier)
        } else if (svtype == "INV") {
            if (end < 0) {
                printf("ERROR: Missing END for INV at %s:%s\\n", chrom1, start) > "/dev/stderr"
                exit 1
            }
            emit_bnd_pair(base_id, base_id, "A", chrom1, start - 1, "-", chrom1, end, "-", qual, filter, info_base, gt_value, tier)
            emit_bnd_pair(base_id, base_id, "B", chrom1, start, "+", chrom1, end + 1, "+", qual, filter, info_base, gt_value, tier)
        } else if (svtype == "TRA") {
            if (end < 0 || chr2 == "") {
                printf("ERROR: Missing CHR2/END for TRA at %s:%s\\n", chrom1, start) > "/dev/stderr"
                exit 1
            }
            if (strands !~ /^[+-][+-]\$/) {
                printf("ERROR: Missing or invalid STRANDS/CT for TRA at %s:%s\\n", chrom1, start) > "/dev/stderr"
                exit 1
            }
            emit_bnd_pair(base_id, base_id, "A", chrom1, start, substr(strands, 1, 1), chr2, end, substr(strands, 2, 1), qual, filter, info_base, gt_value, tier)
        } else if (svtype == "BND") {
            info_base = add_field(info_base, "SVTYPE=BND")
            info_base = add_field(info_base, "TIER=" tier)
            print chrom1, clamp_pos(start), base_id, "N", raw_alt, qual, filter, info_base, "GT", gt_value
        } else {
            printf("WARN: Skipping unsupported SVTYPE=%s for %s\\n", svtype, base_id) > "/dev/stderr"
        }
    }' ${vcf} > ${output_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | sed -n '1p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "output"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
