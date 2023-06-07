include {   BOWTIE2_BUILD as BOWTIE2_HOST_IDX   } from "../../modules/Index_build.nf"
include {   GUNZIP as GUNZIP_HOST               } from '../../modules/gunzip'
include {   GUNZIP as GUNZIP_GFF                } from '../../modules/gunzip'

workflow GENOME_HOST {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //

    if (params.genome_host.endsWith('.gz')) {
        GUNZIP_HOST (
            [ [:], params.genome_host ]
        )
        ch_genome_host    = GUNZIP_HOST.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_HOST.out.versions)
    }

    //
    // Uncompress GFF annotation file
    //

    ch_gff_host = Channel.empty()
    if (params.gff_host) {
        if (params.gff_host.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], params.gff_host ]
            )
            ch_gff_host = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        }
    }

    ch_genome_host = Channel.value(file(params.genome_host))
    
    BOWTIE2_HOST_IDX (
        ch_genome_host.map { [ [:], it ] },
        )

    ch_bt2_index = BOWTIE2_HOST_IDX.out.index
    ch_versions  = ch_versions.mix(BOWTIE2_HOST_IDX.out.versions)

    emit:
    genome_host          = ch_genome_host          // path: genome.fasta
    gff_host             = ch_gff_host             // path: genome.gff_host
    bowtie2_index        = ch_bt2_index            // path: bowtie2/index/

    versions             = ch_versions             // channel: [ versions.yml ]
}

