// Validate input parameters
//WorkflowMetamine.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.genome_host, params.genome_ref ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


include {    HOST_ALIGN_BOWTIE2        } from '../subworkflows/bowtie2/align_sample2host'
include {    REF_ALIGN_BOWTIE2         } from '../subworkflows/bowtie2/align_sample2ref'
include {    INPUT_CHECK                                  } from '../modules/input_check'
include {    GENOME_HOST                                  } from '../subworkflows/bowtie2/index-host'
include {    GENOME_REF                                   } from '../subworkflows/bowtie2/index-ref'

workflow METAMINE {

    ch_versions = Channel.empty()

    //ch_bam                      = Channel.empty()
    //ch_bai                      = Channel.empty()
    //ch_bowtie2_multiqc          = Channel.empty()
    //ch_bowtie2_flagstat_multiqc = Channel.empty()

    INPUT_CHECK (ch_input)
    ch_input_check = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    GENOME_HOST ()
    ch_versions = ch_versions.mix(GENOME_HOST.out.versions)

    GENOME_REF ()
    ch_versions = ch_versions.mix(GENOME_REF.out.versions)

    HOST_ALIGN_BOWTIE2 (
            ch_input_check,
            GENOME_HOST.out.bowtie2_index,
            params.save_unaligned,
            params.sort_bam,
            GENOME_HOST.out.genome_host
        )
    ch_versions = ch_versions.mix(HOST_ALIGN_BOWTIE2.out.versions)
    
    ch_meta = HOST_ALIGN_BOWTIE2.out.fastq
    
    REF_ALIGN_BOWTIE2 (
            ch_meta,
            GENOME_REF.out.bowtie2_index,
            params.save_unaligned,
            params.sort_bam,
            GENOME_REF.out.genome_ref
        )
        
    ch_teste = REF_ALIGN_BOWTIE2.out.fastq
    ch_teste.dump(tag: 'meta')    
    ch_versions  = ch_versions.mix(REF_ALIGN_BOWTIE2.out.versions)
}
    




























/*include {BOWTIE2_BUILD as BOWTIE2_HOST_IDX } from "../modules/Index_build.nf"
include { GUNZIP as GUNZIP_FASTA           } from '../modules/gunzip'
include { GUNZIP as GUNZIP_GFF             } from '../modules/gunzip'

workflow METAMINE {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //

    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], params.fasta ]
        )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }

    //
    // Uncompress GFF annotation file
    //

    ch_gff = Channel.empty()
    if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], params.gff ]
            )
            ch_gff      = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        }
    }

    ch_bowtie2_index = Channel.empty()

    ch_fasta = Channel.value(file(params.fasta))
    
    BOWTIE2_HOST_IDX (
        ch_fasta.map { [ [:], it ] },
        )

    ch_bt2_index = BOWTIE2_HOST_IDX.out.index
    ch_versions  = ch_versions.mix(BOWTIE2_HOST_IDX.out.versions)

    emit:
    fasta                = ch_fasta                // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    bowtie2_index        = ch_bowtie2_index        // path: bowtie2/index/

    versions             = ch_versions             // channel: [ versions.yml ]
}*/

