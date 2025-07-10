#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Remove_Common_Flora_rRNA_reads
} from './modules/local/remove_rRNA_from_reads/main.nf'

workflow Remove_rRNA_Reads_Workflow {
    take:
    input_ch  // tuple(sample_id, sr_file, lr_file, cpus, mem)

    main:
    input_ch.view { "🔬 Input to Remove_Common_Flora_rRNA_reads: $it" }
    Remove_Common_Flora_rRNA_reads(input_ch)

    emit:
    rRNA_host_remove_short = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_short
    rRNA_host_remove_long  = Remove_Common_Flora_rRNA_reads.out.rRNA_host_remove_long
    rRNA_all_files_ch      = Remove_Common_Flora_rRNA_reads.out.rRNA_all_files_ch
}
