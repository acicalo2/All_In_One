#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// QC MODULE: Mode-aware processes

process Interleave {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/interleave_fastq/", mode: 'copy'
    label 'optimized_qc_workflow'

    input:
    tuple val(sample_id), val(fastq_1), val(fastq_2), val(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output:
    tuple val(sample_id), file("${sample_id}_IR.fastq.gz"), emit: interleave_ch

    when:
    mode in ['short', 'hybrid']

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate bbmap
    reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t >> Interleave.log 2>&1
    """
}

process Fastqc {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/fastqc/${sample_id}/${point}/", mode: 'copy'
    label 'optimized_qc_workflow'

    input:
    tuple val(sample_id), file(short_reads), file(long_read), val(point), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output:
    tuple val(sample_id), file("*.zip"), emit: fastqc_ch

    script:
    def short_reads_str = short_reads ? (short_reads instanceof List ? short_reads.join(' ') : short_reads) : ""
    def long_read_str = long_read ?: ""

    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"

    conda activate md

    if [ "${mode}" == "short" ] || [ "${mode}" == "hybrid" ]; then
        [ -n "${short_reads_str}" ] && fastqc --outdir . ${short_reads_str}
    fi

    if [ "${mode}" == "long" ] || [ "${mode}" == "hybrid" ]; then
        [ -n "${long_read_str}" ] && fastqc --outdir . ${long_read_str} 
    fi
    """
}


process Quality_Control {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    //label 'optimized_qc_workflow'

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 


    output:
    tuple val(sample_id), file("${sample_id}_fastp_R1.fastq.gz"), file("${sample_id}_fastp_R2.fastq.gz"), optional: true, emit: trimmed_short_ch
    tuple val(sample_id), file("${sample_id}_fastp_merged.fastq.gz"), optional: true, emit: interleaved_trimmed_short_ch
    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/

    if [ "${mode}" == "short" ] || [ "${mode}" == "hybrid" ]; then
        conda activate md
        fastp --verbose -xyp --dedup --thread=${task.cpus} --report_title "${sample_id}_fastp_report" -q 20 -e 20 --cut_front --cut_tail -w 16 \
              -j ${sample_id}_fastp_ILLUMINA.json \
              -h ${sample_id}_fastp_ILLUMINA.html \
              --in1 ${fastq_1} --in2 ${fastq_2} \
              --out1 ${sample_id}_fastp_R1.fastq.gz \
              --out2 ${sample_id}_fastp_R2.fastq.gz >> fastp.log 2>&1

        conda activate bbmap
        reformat.sh in1=${sample_id}_fastp_R1.fastq.gz in2=${sample_id}_fastp_R2.fastq.gz out=${sample_id}_fastp_merged.fastq.gz ow=t
    fi

    echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished
    """
}


process PoreChop {
    tag {sample_id} 
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    //label 'optimized_qc_workflow'
    //conda './md.yml'

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output: 
    tuple val(sample_id), file("${sample_id}.LR.trimmed.fastq.gz"), emit: porechop_ch
    
    when:
    params.longreads || params.hybrid
    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate dragonflye
    porechop -i ${long_read} \
             -o ${sample_id}.LR.trimmed.fastq.gz  \
             --format fastq.gz \
             -t ${task.cpus} \
             --require_two_barcodes \
             --discard_unassigned --discard_middle >> PoreChop.log 2>&1
    """

}


process Multiqc_QC_Stats {
    publishDir "${params.outdir}/${params.project_id}/multiqc", mode: 'copy'
    label 'optimized_qc_workflow'

    input:
    path pretrim_fastqc_ch
    path posttrim_fastqc_ch

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate multiqc

    mkdir -p ${params.outdir}/${params.project_id}/multiqc/pretrim
    multiqc ${pretrim_fastqc_ch.join(' ')} --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/pretrim

    mkdir -p ${params.outdir}/${params.project_id}/multiqc/post_trim
    multiqc ${posttrim_fastqc_ch.join(' ')} --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/post_trim

    conda activate nextflow
    mkdir -p ${params.outdir}/${params.project_id}/qc_stats
    Rscript ${params.scripts}/create_qc_stats.R \
        -i ${params.outdir}/${params.project_id}/multiqc/pretrim/multiqc_data/multiqc_general_stats.csv \
        -p ${params.outdir}/${params.project_id}/multiqc/post_trim/multiqc_data/multiqc_general_stats.csv \
        -o ${params.outdir}/${params.project_id}/qc_stats/

    echo "Post Trim Multiqc Complete" > post_trim_multiqc.finished
    """
}

//process Read_Distribution {
//    publishDir "${params.outdir}/${params.project_id}/read_distribution", mode: 'copy'
//    label 'optimized_qc_workflow'

//    input:
//    tuple val(sample_id), file(short_reads), file(long_read), val(point), val(mode), val(cpus), val(mem)

//    cpus {cpus} // setting slurm allocation dynamically
//    memory {mem} // setting slurm allocation dynamically
//    output:

//    script:
//    """

//    """
//}
