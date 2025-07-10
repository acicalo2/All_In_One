#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Map 2 Reference
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Map_Reads_2_RefSeq {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/", mode: 'copy'
    //label 'optimized_map2reads' // can erase if new tuple works

    input:
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read), val(mode), val(cpus), val(mem)

    cpus {cpus}  // setting slurm allocation dynamically
    memory {mem} // setting slurm allocation dynamically

    output:
    tuple val(sample_id), file("${sample_id}_host_removed_sr.fastq.gz"), optional: true, emit: host_removed_fq_ch_sr
    tuple val(sample_id), file("${sample_id}_host_removed_LR.fastq.gz"), optional: true, emit: host_removed_fq_ch_lr

    when:
    params.map2reads

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/

    if [[ "${mode}" == "long" || "${mode}" == "hybrid" ]]; then
        minimap2 -ax map-ont -t $task.cpus ${params.host_db} \\
            ${long_read} > ${sample_id}_host_removed_LR.sam
        samtools fastq -f 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_LR.fastq
        pigz *.fastq
        echo "Host Removal LR complete" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/LR_hostRem.finished
    fi

    if [[ "${mode}" == "short" || "${mode}" == "hybrid" ]]; then
        reformat.sh in1=${fastq_1} in2=${fastq_2} out=${sample_id}_IR.fastq.gz ow=t
        minimap2 -ax sr -t $task.cpus ${params.host_db} \\
            ${sample_id}_IR.fastq.gz > ${sample_id}_host_removed_sr.sam
        samtools fastq -f 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_removed_sr.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_sr.fastq
        pigz *.fastq
        echo "Host Removal SR complete" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/SR_hostRem.finished
    fi
    """
}
