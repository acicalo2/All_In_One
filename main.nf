#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ─── Imports ──────────────────────────────────────────────────
include { QC_Workflow }                        from './workflows/qc.nf'
include { MapReads_2_RefSeq }                  from './workflows/MapReads2RefSeq.nf'
include { Assembly_Workflow }                  from './workflows/assembly.nf'
include { Split_Interleave }                   from './workflows/modules/local/split_interleave/main.nf'
include { BlastX_Reads_Workflow }              from './workflows/BlastX_Reads.nf'
include { Split_BamFiles_Workflow }            from './workflows/dorado_polish.nf'
include { Dorado_Aligner }                     from './workflows/modules/local/dorado_polish/main.nf'
include { Dorado_Polisher }                    from './workflows/modules/local/dorado_polish/main.nf'
include { DAA2INFO_contigs_daa_file }          from './workflows/modules/local/blastx/main.nf'
include { Meganize_ShortReads_BlastX }         from './workflows/modules/local/blastx/main.nf'
include { Meganize_LongReads_BlastX }          from './workflows/modules/local/blastx/main.nf'
include { Hecatomb} 		                   from './workflows/modules/local/hecatomb/main.nf'
include { Quast } 	                           from './workflows/modules/local/quast/main.nf'
include { BlastX_Dragonflye_contigs }          from './workflows/modules/local/blastx/main.nf'
include { BlastX_metaspades_contigs }          from './workflows/modules/local/blastx/main.nf'
include { BlastX_spades_contigs }              from './workflows/modules/local/blastx/main.nf'
include { BlastX_unicycler_contigs }           from './workflows/modules/local/blastx/main.nf'
include { Meganize_BlastX_Contigs_Dragonflye } from './workflows/modules/local/megan/main.nf'
include { Meganize_BlastX_Contigs_Metaspades } from './workflows/modules/local/megan/main.nf'
include { Meganize_BlastX_Contigs_Unicycler }  from './workflows/modules/local/megan/main.nf'
include { Parse_BlastX_Dragonflye }            from './workflows/modules/local/megan/main.nf'
include { Parse_BlastX_metaspades }            from './workflows/modules/local/megan/main.nf'
include { Parse_BlastX_unicycler }             from './workflows/modules/local/megan/main.nf'
include { MMSEQ_contigs_against_NT_Dragonflye }from './workflows/modules/local/mmseq/main.nf'
include { MMSEQ_contigs_against_NT_metaspades }from './workflows/modules/local/mmseq/main.nf'
include { MMSEQ_contigs_against_NT_unicycler } from './workflows/modules/local/mmseq/main.nf'
include { Parse_MMSEQ_dragonflye_contigs }     from './workflows/modules/local/mmseq/main.nf'
include { Parse_MMSEQ_metaspades_contigs }     from './workflows/modules/local/mmseq/main.nf'
include { Parse_MMSEQ_unicycler_contigs }      from './workflows/modules/local/mmseq/main.nf'
include { Remove_rRNA_Reads_Workflow }         from './workflows/Remove_rRNA_wf.nf'
include { CheckM_Dragonflye }                  from './workflows/modules/local/checkm/main.nf'
include { CheckM_Dragonflye_Medaka }           from './workflows/modules/local/checkm/main.nf'
include { CheckM_Unicycler }                   from './workflows/modules/local/checkm/main.nf'

// ─── Main Workflow ─────────────────────────────────────
workflow {

    def raw_samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def fq1 = row.fastq_1?.contains('No_Read') ? null : file(row.fastq_1, checkIfExists: false)
            def fq2 = row.fastq_2?.contains('No_Read') ? null : file(row.fastq_2, checkIfExists: false)
            def lr  = row.long_read?.contains('No_Read') ? null : file(row.long_read, checkIfExists: false)

            def mode = fq1 && fq2 && lr ? 'hybrid'
                     : fq1 && fq2       ? 'short'
                     : lr               ? 'long'
                     : 'none'

            // Estimate size and set resources for QC
            def read_files = [fq1, fq2, lr].findAll()
            def total_size_gb = read_files.sum { it?.size() ?: 0 } / 1e9

            // Dynamically assign cpus and memory
            def cpus = total_size_gb < 1 ? 2 :
                       total_size_gb < 5 ? 4 :
                       total_size_gb < 10 ? 6 : 8
            def mem = total_size_gb < 1 ? '8 GB' :
                       total_size_gb < 5 ? '16 GB' :
                       total_size_gb < 10 ? '24 GB' : '32 GB'
                       
            tuple(row.sample_id, fq1, fq2, lr, mode, cpus, mem)
        }
        .filter { it[4] != 'none' }

    raw_samples_ch.view { "Raw sample data: $it" }

    // QC Workflow
    QC_Workflow(raw_samples_ch)


    // Ensure raw_samples_ch is consumed first (now used by QC_Workflow)

    def pipeline_mode_ch = raw_samples_ch
        .map { it[4] }
        .distinct()
        .collect()
        .map { modes ->
            if (modes.contains('hybrid')) return 'hybrid'
            else if (modes == ['short']) return 'short'
            else if (modes == ['long']) return 'long'
            else throw new IllegalStateException("Unsupported pipeline modes: ${modes}")
        }
    pipeline_mode_ch.view { "Detected pipeline mode: $it" }
    
    def short_ch
    def long_ch
    def short_ch_final
    def long_ch_final
    def sr_ch
    def lr_ch

    if (params.map2reads) {
        MapReads_2_RefSeq(
            QC_Workflow.out.trimmed_short_ch,
            QC_Workflow.out.trimmed_long_ch
        )   

        MapReads_2_RefSeq.out.host_removed_fq_ch_sr.view { " host_removed_fq_ch_sr: $it" }
        MapReads_2_RefSeq.out.host_removed_fq_ch_lr.view { " host_removed_fq_ch_lr: $it" } 

        if (params.remove_rRNA_reads) {
            def rrna_input_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .map { id, sr -> tuple(id, [sr, null]) }
                .mix(
                    MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                        .map { id, lr -> tuple(id, [null, lr]) }
                )
                .groupTuple()
                .map { id, values ->
                    def sr_file = values.collect { it[0] }.find { it != null }
                    def lr_file = values.collect { it[1] }.find { it != null }

                    def files = [sr_file, lr_file].findAll()
                    def total_size_gb = files.sum {it?.size() ?: 0 } / 1e9

                    def cpus = total_size_gb < 1 ? 8 :
                               total_size_gb < 5 ? 16 :
                               total_size_gb < 10 ? 32 :
                               total_size_gb < 20 ? 48 :
                               total_size_gb < 40 ? 64 : 80
                    
                    def mem = total_size_gb < 1 ? '32 GB' :
                               total_size_gb < 5 ? '64 GB' :
                               total_size_gb < 10 ? '128 GB' :
                               total_size_gb < 20 ? '256 GB' :
                               total_size_gb < 40 ? '384 GB' : '512 GB'

                    tuple(id, sr_file, lr_file, cpus, mem)
                }
            rrna_input_ch.view { "🧪RNA input tuple: $it" }    

            Remove_rRNA_Reads_Workflow(rrna_input_ch)   

            short_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            long_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_long
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }  

            long_ch_final = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_long
                .map { sample_id, lr -> tuple(sample_id, lr) }  

            def split_input_ch = Remove_rRNA_Reads_Workflow.out.rRNA_host_remove_short
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            Split_Interleave(split_input_ch)    

            short_ch_final = Split_Interleave.out.split_interleave_ch
                .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }  

        } else {
            short_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            long_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }  

            long_ch_final = MapReads_2_RefSeq.out.host_removed_fq_ch_lr
                .map { sample_id, lr -> tuple(sample_id, lr) }  

            def split_input_ch = MapReads_2_RefSeq.out.host_removed_fq_ch_sr
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') } 

            Split_Interleave(split_input_ch)    

            short_ch_final = Split_Interleave.out.split_interleave_ch
                .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }
        }   

    } else {
            short_ch = QC_Workflow.out.interleaved_trimmed_short_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') }     

            long_ch = QC_Workflow.out.trimmed_long_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'long') }      

            long_ch_final = QC_Workflow.out.trimmed_long_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq) }      

            def split_input_ch = QC_Workflow.out.interleaved_trimmed_short_ch
                .filter { sample_id, fq -> fq != null }
                .map { sample_id, fq -> tuple(sample_id, fq, 'short') }     

            Split_Interleave(split_input_ch)        

            short_ch_final = Split_Interleave.out.split_interleave_ch
                .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }
    }
    // ─── Enforce tuple structure ──────────────────────────────────────────────   

    // Ensure short_ch_final is a 3-tuple (sample_id, R1, R2)
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }    

    // Ensure long_ch_final is a 2-tuple (sample_id, LR)
    long_ch_final = long_ch_final.map { it -> tuple(it[0], it[1]) } 

    // Optional: debug before building hybrid
    short_ch_final.view { "🧪short_ch_final (tuple): $it" }
    long_ch_final.view  { "🧪long_ch_final (tuple): $it" }  

    // ─── Build conditional assembly input ─────────────────────────────────────   

    // Short-only input channel: (sample_id, R1, R2, null, 'short')
    def short_only_ch = short_ch_final.map { id, r1, r2 ->
        tuple(id, r1, r2, null, 'short')
    }   

    // Long-only input channel: (sample_id, null, null, LR, 'long')
    def long_only_ch = long_ch_final.map { id, lr ->
        tuple(id, null, null, lr, 'long')
    }   

    // First, ensure sample_id is a plain String in both channels
    short_ch_final = short_ch_final.map { it -> tuple(it[0], it[1], it[2]) }
    long_ch_final  = long_ch_final.map { it -> tuple(it[0], it[1]) }

    short_ch_final.view { " short tuple: $it (${it.getClass()})" }
    long_ch_final.view  { " long tuple:  $it (${it.getClass()})" }
    def hybrid_ch = short_ch_final
        .map { id, r1, r2 -> tuple(id.toString(), ['short', r1, r2]) }
        .mix(
            long_ch_final.map { id, lr -> tuple(id.toString(), ['long', lr]) }
        )
        .groupTuple()
        .filter { id, values ->
            values.find { it[0] == 'short' } && values.find { it[0] == 'long' }
        }
        .map { id, values ->
            def short_vals = values.find { it[0] == 'short' }
            def long_vals  = values.find { it[0] == 'long' }    

            def r1 = short_vals[1]
            def r2 = short_vals[2]
            def lr = long_vals[1]   

            tuple(id, r1, r2, lr, 'hybrid')
        }
        
    hybrid_ch.view { it -> println " hybrid_ch: ${it} (${it.getClass()})" }

    def assembly_input_ch   

    if (params.hybrid) {
        assembly_input_ch = hybrid_ch
    } else if (params.shortreads) {
        assembly_input_ch = short_only_ch
    } else if (params.longreads) {
        assembly_input_ch = long_only_ch
    } else {
        error "You must set one of: --shortreads, --longreads, or --hybrid"
    }
    
    assembly_input_ch.view { " assembly_input_ch: $it" }

    // ─── Call the Assembly Workflow ────────────────────────────────────────  

    Assembly_Workflow(
        assembly_input_ch
    )
    short_ch.view { "short_ch: $it" }
    long_ch.view  { "long_ch: $it" }
    if (params.medaka) {
        Assembly_Workflow.out.dragonflye_medaka_assembly_ch.view { " Dragonflye_Medaka output: $it" }
    }

    if ( params.dorado_polish ) {
        log.info " Running Split_BamFiles_Workflow"
        log.info " Using BAM file: ${params.bam_file}"

        def dorado_split_input_ch = Channel.fromPath(params.bam_file, checkIfExists: true)
        dorado_split_input_ch.view { " Input to Split_BamFiles_Workflow: $it" }

        Split_BamFiles_Workflow(dorado_split_input_ch)

        // 🛠 Extract sample_id from BAM filename
        def dorado_split_tuples_ch = Split_BamFiles_Workflow.out.dorado_split_output_ch
            .flatten()
            .map { bam_file ->
                def sample_id = bam_file.getBaseName().replaceFirst(/\.bam$/, '')
                tuple(sample_id, bam_file)
            }

        dorado_split_tuples_ch.view { " BAM as tuples: $it" }
        Assembly_Workflow.out.dragonflye_assembly_ch.view { " Dragonflye output: $it" }

        //  Join on sample_id and feed to Dorado_Aligner
        def dorado_aligner_input_ch = dorado_split_tuples_ch
            .join(Assembly_Workflow.out.dragonflye_assembly_ch)
            .map { sample_id, bam_file, contigs_fasta -> tuple(sample_id, bam_file, contigs_fasta) }

        dorado_aligner_input_ch.view { "🧪Joined input to Dorado_Aligner: $it" }

        Dorado_Aligner(dorado_aligner_input_ch)
        Dorado_Aligner.out.dorado_aligner_ch.view { " Dorado aligner output: $it" }
        // Step 1: Get original BAM input (before alignment)
        def dorado_input_bam_ch = dorado_aligner_input_ch
            .map { sample_id, bam_file, contigs_fasta -> tuple(sample_id, bam_file) }       

        // Step 2: Wait for aligner to produce .bam.bai
        def dorado_bai_ch = Dorado_Aligner.out.dorado_aligner_ch
            .map { sample_id, bai_file -> tuple(sample_id, bai_file) }      

        // Step 3: Join BAM + .bai on sample_id
        def dorado_bam_with_index_ch = dorado_input_bam_ch
            .join(dorado_bai_ch)
            .map { sample_id, bam_file, bai_file -> tuple(sample_id, bam_file, bai_file) }      

        // Step 4: Join with contigs.fa from Dragonflye
        def dorado_polisher_input_ch = dorado_bam_with_index_ch
            .join(Assembly_Workflow.out.dragonflye_assembly_ch)
            .map { sample_id, bam_file, bai_file, contigs_fasta -> tuple(sample_id, bam_file, bai_file, contigs_fasta) }        

        dorado_polisher_input_ch.view { "Input to Dorado_Polisher: $it" }     

        Dorado_Polisher(dorado_polisher_input_ch)
        Dorado_Polisher.out.dorado_polisher_ch.view { " Dorado_Polisher output: $it" }
    }


    // CheckM
    if ( params.dragonflye) {
        CheckM_Dragonflye(Assembly_Workflow.out.dragonflye_assembly_ch)
    }
    if (params.medaka) {
        CheckM_Dragonflye_Medaka(Assembly_Workflow.out.dragonflye_medaka_assembly_ch)
    }

    if (params.unicycler) {
        CheckM_Unicycler(Assembly_Workflow.out.unicycler_assembly_ch)
    }


    if ( params.nanopore_assembly ) {
        
	// Remove 'mode' field from assembly_input_ch before passing to Hecatomb
	def hecatomb_input_ch = assembly_input_ch.map { id, r1, r2, lr, mode -> 
	    tuple(id, r1, r2, lr)
	}	

	Hecatomb(hecatomb_input_ch)
        
    // Input to Stage_7a_quast: tuple(sample_id, merged_fasta_path)
	def quast_input_ch = assembly_input_ch
        	.map { sample_id, fq1, fq2, lr, mode -> tuple(sample_id, fq1, fq2, lr) } // remove mode
        	.join(Hecatomb.out.hecatomb_merged_fasta)
        	.map { sample_id, fq1, fq2, lr, merged_fasta ->
            	tuple(sample_id, fq1, fq2, lr, merged_fasta)
      	} 

	Quast(quast_input_ch)
    
    }

    Assembly_Workflow.out.dragonflye_assembly_ch.view { " Dragonflye output: $it" }
    Assembly_Workflow.out.spades_assembly_ch.view { " SPAdes output: $it" }
    Assembly_Workflow.out.metaspades_assembly_ch.view { " MetaSPAdes output: $it" }
    Assembly_Workflow.out.unicycler_assembly_ch.view { " Unicycler output: $it" }
    Assembly_Workflow.out.plasmidspades_assembly_ch.view { " PlasmidSPAdes output: $it" }
    
    if ( params.run_blastx ) {
        if ( params.dragonflye ){
            BlastX_Dragonflye_contigs(Assembly_Workflow.out.dragonflye_assembly_ch)
            MMSEQ_contigs_against_NT_Dragonflye(Assembly_Workflow.out.dragonflye_assembly_ch)
            Parse_MMSEQ_dragonflye_contigs(MMSEQ_contigs_against_NT_Dragonflye.out.mmseqs_dragonflye_contigs_ch)
            BlastX_Dragonflye_contigs.out.dragonflye_blastx_contigs_ch.view { " Dragonflye BlastX Contigs Output: $it" }
            Meganize_BlastX_Contigs_Dragonflye(BlastX_Dragonflye_contigs.out.dragonflye_blastx_contigs_ch)
            Parse_BlastX_Dragonflye(Meganize_BlastX_Contigs_Dragonflye.out.dragonflye_meganize_blastx_contigs_ch)
        }
        if ( params.spades ){
            BlastX_spades_contigs(Assembly_Workflow.out.spades_assembly_ch)
        }
        if ( params.metaspades ){
            BlastX_metaspades_contigs(Assembly_Workflow.out.metaspades_assembly_ch)
            MMSEQ_contigs_against_NT_metaspades(Assembly_Workflow.out.metaspades_assembly_ch)
            Parse_MMSEQ_metaspades_contigs(MMSEQ_contigs_against_NT_metaspades.out.mmseqs_metaspades_contigs_ch)
            BlastX_metaspades_contigs.out.metaspades_blastx_contigs_ch.view { " Metaspades BlastX Contigs Output: $it" }
            Meganize_BlastX_Contigs_Metaspades(BlastX_metaspades_contigs.out.metaspades_blastx_contigs_ch)
            Parse_BlastX_metaspades(Meganize_BlastX_Contigs_Metaspades.out.metaspades_meganize_blastx_contigs_ch)
        }
        if ( params.unicycler ){
            BlastX_unicycler_contigs(Assembly_Workflow.out.unicycler_assembly_ch)
            MMSEQ_contigs_against_NT_unicycler(Assembly_Workflow.out.unicycler_assembly_ch)
            Parse_MMSEQ_unicycler_contigs(MMSEQ_contigs_against_NT_unicycler.out.mmseqs_unicycler_contigs_ch)
            BlastX_unicycler_contigs.out.unicycler_blastx_contigs_ch.view { " Unicycler BlastX Contigs Output: $it" }
            Meganize_BlastX_Contigs_Unicycler(BlastX_unicycler_contigs.out.unicycler_blastx_contigs_ch)
            Parse_BlastX_unicycler(Meganize_BlastX_Contigs_Unicycler.out.unicycler_meganize_blastx_contigs_ch)
        } 
        //if ( params.plasmidspades ){
        //    BlastX_Contigs_Workflow(Assembly_Workflow.out.plasmidspades_assembly_ch)
        //} 
        BlastX_Reads_Workflow(short_ch, long_ch)
	// Assume you already have these defined:
	BlastX_Reads_Workflow.out.short_reads_blastx_channel.view { " SHORT: $it" }
	BlastX_Reads_Workflow.out.long_reads_blastx_channel.view  { " LONG: $it" }
	//short_reads_ch_blastx = BlastX_Reads_Workflow.out.short_reads_blastx_channel
	//long_reads_ch_blastx  = BlastX_Reads_Workflow.out.long_reads_blastx_channel	

	if (params.shortreads || params.hybrid) {
	BlastX_Reads_Workflow.out.short_reads_blastx_channel
	    .ifEmpty(Channel.empty())
	    .view { " Short DAA: $it" }
	    | Meganize_ShortReads_BlastX	
        }
        if (params.longreads || params.hybrid) {
	BlastX_Reads_Workflow.out.long_reads_blastx_channel
	    .ifEmpty(Channel.empty())
	    .view { " Long DAA: $it" }
	    | Meganize_LongReads_BlastX
	}

    }
}

