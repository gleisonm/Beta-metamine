#!/usr/bin/env nextflow

log.info """\
    
    ╔═══════════════════════════════════════╗
    ║ ▂▃▄▅▆▇█▓▒░ M E T A M I N E ░▒▓█▇▆▅▄▃▂ ║
    ║ N E X T F L O W   -   P I P E L I N E ║
    ╚══════════════════════════════════════ ☣
    genome       : ${params.fasta}
    reads        : 
    outdir       : ${params.outdir}
    """
    .stripIndent()

include { METAMINE } from './workflow/metamine.nf'

    workflow {
        METAMINE()
    }