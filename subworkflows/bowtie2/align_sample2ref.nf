//
// Align sample to ref with Bowtie2
//
include { BAM_SORT_STATS_SAMTOOLS } from '../../modules/bam_sort_stats_samtools'
include { BOWTIE2_ALIGN as ALIGN_2_REF_BOWTIE2   } from '../../modules/bowtie2_align'


workflow REF_ALIGN_BOWTIE2 {
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
    ALIGN_2_REF_BOWTIE2 ( ch_reads, ch_index, save_unaligned, sort_bam)
    ch_versions = ch_versions.mix(ALIGN_2_REF_BOWTIE2.out.versions.first())
    
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    //BOWTIE2_ALIGN.out.aligned.dump(tag: 'out')
    //BAM_SORT_STATS_SAMTOOLS ( BOWTIE2_ALIGN.out.aligned, ch_fasta )
    //ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam             = ALIGN_2_REF_BOWTIE2.out.bam          // channel: [ val(meta), bam   ]
    logs            = ALIGN_2_REF_BOWTIE2.out.log          // channel: [ val(meta), log   ]
    fastq           = ALIGN_2_REF_BOWTIE2.out.fastq        // channel: [ val(meta), fastq ]

    //bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    //bai              = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    //csi              = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    //stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    //flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    //idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                      // channel: [ versions.yml ]
}