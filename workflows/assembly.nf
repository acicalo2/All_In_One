#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    Spades;
    MetaSPAdes_Assembly;
    Plasmid_Spades;
    Dragonflye as Dragonflye_LR;
    Dragonflye as Dragonflye_Hybrid;
    Dragonflye_Medaka as Dragonflye_Medaka_LR;
    Dragonflye_Medaka as Dragonflye_Medaka_Hybrid;
    Unicycler_Assembly;
} from './modules/local/assembler/main_working.nf'

workflow Assembly_Workflow {

    take:
    assembly_input_ch

    main:

    assembly_input_ch.view { "✈️ Assembly input (received): $it" }

    short_input_ch  = assembly_input_ch.filter { it[4] in ['short', 'hybrid'] }
    long_input_ch   = assembly_input_ch.filter { it[4] == 'long' }
    hybrid_input_ch = assembly_input_ch.filter { it[4] == 'hybrid' }

    if (params.spades) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2) }
            .set { spades_input_ch }

        Spades(spades_input_ch)
    }

    if (params.plasmidspades) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2) }
            .set { plasmid_spades_input_ch }

        Plasmid_Spades(plasmid_spades_input_ch)
    }

    if (params.metaspades) {
        short_input_ch
            .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2) }
            .set { metaspades_input_ch }

        MetaSPAdes_Assembly(metaspades_input_ch)
    }


    def dragonflye_ch = Channel.empty()

    if (params.dragonflye) {
        if (long_input_ch) {
            long_input_ch
                .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, null, null, lr) }
                .set { dragonflye_input_ch_long }

            Dragonflye_LR(dragonflye_input_ch_long)
            dragonflye_ch = Dragonflye_LR.out.dragonflye_assembly_ch
        }

        if (hybrid_input_ch) {
            hybrid_input_ch
                .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2, lr) }
                .set { hybrid_dragonflye_input_ch }

            Dragonflye_Hybrid(hybrid_dragonflye_input_ch)
            dragonflye_ch = dragonflye_ch.mix(Dragonflye_Hybrid.out.dragonflye_assembly_ch)
        }
    }

    def dragonflye_medaka_ch = Channel.empty()  

    if (params.medaka) {
        if (long_input_ch) {
            long_input_ch
                .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, null, null, lr) }
                .set { dragonflye_medaka_input_long }   

            Dragonflye_Medaka_LR(dragonflye_medaka_input_long)
            dragonflye_medaka_ch = Dragonflye_Medaka_LR.out.dragonflye_assembly_medaka_ch
        }   

        if (hybrid_input_ch) {
            hybrid_input_ch
                .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2, lr) }
                .set { dragonflye_medaka_input_hybrid } 

            Dragonflye_Medaka_Hybrid(dragonflye_medaka_input_hybrid)
            dragonflye_medaka_ch = dragonflye_medaka_ch.mix(Dragonflye_Medaka_Hybrid.out.dragonflye_assembly_medaka_ch)
        }
    }

    def unicycler_ch = Channel.empty()

    if (params.unicycler) {
        hybrid_input_ch
            .map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2, lr) }
            .set { unicycler_input_ch }

        Unicycler_Assembly(unicycler_input_ch)
        unicycler_ch = Unicycler_Assembly.out.unicycler_assembly_ch
    }

    emit:
    spades_assembly_ch        = params.spades        ? Spades.out.spades_assembly_ch : Channel.empty()
    metaspades_assembly_ch    = params.metaspades    ? MetaSPAdes_Assembly.out.metaspades_assembly_ch : Channel.empty()
    dragonflye_assembly_ch    = dragonflye_ch
    unicycler_assembly_ch     = unicycler_ch
    plasmidspades_assembly_ch = params.plasmidspades ? Plasmid_Spades.out.plasmidspades_assembly_ch : Channel.empty()
    dragonflye_medaka_assembly_ch = params.medaka ? dragonflye_medaka_ch : Channel.empty()
}

