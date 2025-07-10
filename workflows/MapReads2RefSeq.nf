#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Map_Reads_2_RefSeq
} from './modules/local/map2refseq/main.nf'

workflow MapReads_2_RefSeq {
    take:
    trimmed_short_ch
    trimmed_long_ch

    main:
    // Convert short reads to 5-tuple
    //def short_ch = trimmed_short_ch.map { sample_id, fq1, fq2 -> 
    //    tuple(sample_id, fq1, fq2, null, 'short') 
    //}
    // * TEST* Apply file-size-based resource scaling
    def short_ch = trimmed_short_ch.map { sample_id, fq1, fq2 -> 
                   def total_size = [fq1, fq2].findAll().sum { it.size() }
                   def size_gb = total_size / 1e9

                   def cpus = size_gb < 1 ? 2 : 
                              size_gb < 5 ? 4 : 
                              size_gb < 10 ? 8 : 16

                   def mem = size_gb < 1 ? '8GB' : 
                              size_gb < 5 ? '16GB' : 
                              size_gb < 10 ? '32GB' : '64GB'   
                   tuple(sample_id, fq1, fq2, null, 'short', cpus, mem) 
    }                      
    // Convert long reads to 5-tuple
    def long_ch = trimmed_long_ch.map { sample_id, long_read -> 
        def size_gb = long_read.size() / 1e9
        def cpus = size_gb < 1 ? 2 : 
                   size_gb < 5 ? 4 : 
                   size_gb < 10 ? 8 : 16

        def mem =  size_gb < 1 ? '8GB' : 
                   size_gb < 5 ? '16GB' : 
                   size_gb < 10 ? '32GB' : '64GB'   
        tuple(sample_id, null, null, long_read, 'long', cpus, mem) 
    }

    // Merge both types
    def merged_ch = short_ch.mix(long_ch)

    // Send to mapping process
    Map_Reads_2_RefSeq(merged_ch)

    emit:
    host_removed_fq_ch_sr = Map_Reads_2_RefSeq.out.host_removed_fq_ch_sr
    host_removed_fq_ch_lr = Map_Reads_2_RefSeq.out.host_removed_fq_ch_lr
}
