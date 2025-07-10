#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    Quality_Control;
    PoreChop
} from './modules/local/qc/main.nf'

workflow QC_Workflow {
    take:
    input_samples_ch

    main:
    Quality_Control(input_samples_ch)
    PoreChop(input_samples_ch)
    trimmed_short_ch = Quality_Control.out.trimmed_short_ch
    trimmed_long_ch  = PoreChop.out.porechop_ch
    interleaved_trimmed_short_ch  = Quality_Control.out.interleaved_trimmed_short_ch
    //Read_Distribution(input_samples_ch,interleaved_trimmed_short_ch,trimmed_long_ch)
    emit:
    trimmed_short_ch = trimmed_short_ch
    trimmed_long_ch  = trimmed_long_ch
    interleaved_trimmed_short_ch = interleaved_trimmed_short_ch
}
