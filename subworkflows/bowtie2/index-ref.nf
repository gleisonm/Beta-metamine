include {   BOWTIE2_BUILD as BOWTIE2_REF_IDX    } from "../../modules/Index_build.nf"
include {   GUNZIP as GUNZIP_GENOME_REF         } from '../../modules/gunzip'
include {   GUNZIP as GUNZIP_GFF                } from '../../modules/gunzip'

workflow GENOME_REF {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //

    if (params.genome_ref.endsWith('.gz')) {
        GUNZIP_GENOME_REF (
            [ [:], params.genome_ref ]
        )
        ch_genome_ref   = GUNZIP_GENOME_REF.out.gunzip.map { it[1] }
        ch_versions     = ch_versions.mix(GUNZIP_GENOME_REF.out.versions)
    }

    //
    // Uncompress GFF annotation file
    //

    ch_gff_ref = Channel.empty()
    if (params.gff_ref) {
        if (params.gff_ref.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], params.gff_ref ]
            )
            ch_gff_ref  = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        }
    }

    ch_bowtie2_index = Channel.empty()

    ch_genome_ref = Channel.value(file(params.genome_ref))
    
    BOWTIE2_REF_IDX (
        ch_genome_ref.map { [ [:], it ] },
        )

    ch_bt2_index = BOWTIE2_REF_IDX.out.index
    ch_versions  = ch_versions.mix(BOWTIE2_REF_IDX.out.versions)

    emit:
    genome_ref           = ch_genome_ref           // path: genome.genome_ref
    gff_ref              = ch_gff_ref              // path: genome.gff_ref
    bowtie2_index        = ch_bowtie2_index        // path: bowtie2/index/

    versions             = ch_versions             // channel: [ versions.yml ]
}

