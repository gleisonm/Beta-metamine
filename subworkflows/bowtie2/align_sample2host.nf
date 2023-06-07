//
// Align sample to host with Bowtie2
//
include { BAM_SORT_STATS_SAMTOOLS } from '../../modules/bam_sort_stats_samtools'
include { BOWTIE2_ALIGN           } from '../../modules/bowtie2_align'


workflow ALIGN_2_HOST_BOWTIE2 {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_index          // channel: /path/to/bowtie2/index/
    save_unaligned    // val
    sort_bam          // val
    ch_fasta          // channel: /path/to/reference.fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    BOWTIE2_ALIGN ( ch_reads, ch_index, save_unaligned, sort_bam)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

  
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    //BOWTIE2_ALIGN.out.aligned.dump(tag: 'out')
    //BAM_SORT_STATS_SAMTOOLS ( BOWTIE2_ALIGN.out.aligned, ch_fasta )
    //ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam             = BOWTIE2_ALIGN.out.bam          // channel: [ val(meta), bam   ]
    logs            = BOWTIE2_ALIGN.out.log          // channel: [ val(meta), log   ]
    fastq           = BOWTIE2_ALIGN.out.fastq        // channel: [ val(meta), fastq ]

    //bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    //bai              = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    //csi              = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    //stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    //flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    //idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                      // channel: [ versions.yml ]
}