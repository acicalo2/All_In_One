#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Remove_Common_Flora_rRNA_reads {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/rRNA_removed/", mode: 'copy'
    label 'optimized_remove_rRNA_Reads_Workflow'
    
    input:
    tuple val(sample_id), file(sr_file), file(lr_file), val(cpus), val(mem)
    cpus {cpus} // setting slurm allocation dynamically 
    memory {mem} // setting slurm allocation dynamically 

    output:
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz"), optional: true, emit: rRNA_host_remove_short
    tuple val(sample_id), file("${sample_id}_host_contaminant_rRNA_removed_LR.fastq"), optional: true, emit: rRNA_host_remove_long
    tuple val(sample_id), file("*"), optional: true, emit: rRNA_all_files_ch

    when:
    params.remove_rRNA_reads

    script:
    def java_mem = mem.tokenize()[0] + 'g'
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    if [[ -s "$lr_file" ]]; then
        minimap2 -ax map-ont -t ${task.cpus} ${params.host_db_silva} ${lr_file} > ${sample_id}_host_contaminant_rRNA_removed_LR.sam >> minimap2_removeRNA.log 2>&1
        samtools fastq -f 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam > ${sample_id}_host_contaminant_rRNA_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam > ${sample_id}_rRNA_LR.fastq
    fi

    if [[ -s "$sr_file" ]]; then
        bbmap.sh ${params.bbmap_args} \\
            in=${sr_file} \\
            path=${params.bbmap_ref_silva} tossbrokenreads printunmappedcount=t \\
            covstats=${sample_id}.rRNA_covstats.txt \\
            outm=${sample_id}_rRNA.fastq.gz usejni=t \\
            -Xmx${java_mem} \\
            outu=${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz overwrite=true >> bbmap_removeRNA.log 2>&1

        seqtk seq -1 ${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz | gzip -1 > ${sample_id}_host_contaminant_rRNA_removed_R1.fastq.gz
        seqtk seq -2 ${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz | gzip -1 > ${sample_id}_host_contaminant_rRNA_removed_R2.fastq.gz
    fi
    """
}
