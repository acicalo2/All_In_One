#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Dorado Polish Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Dorado_Split_BamFiles {

    tag "bamfile"
    publishDir "${params.outdir}/${params.project_id}/bam_files/", mode: 'copy'
    label 'normal'

    input:
    path(bamfile)

    output:
    path("*.bam")

    script:
    """
    ${params.scripts}/split_rename_bams.sh $bamfile
    """
}

process Dorado_Aligner {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dorado_aligner/", mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), path(bam_file), path(contigs_fasta)

    output:
    tuple val(sample_id), file("${sample_id}*"), emit: dorado_aligner_ch

    script:
    """
    ${params.dorado_software}/dorado aligner ${contigs_fasta} ${bam_file} -o .
    """
}

/*
process Dorado_Aligner {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dorado_aligner/", mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), path(bam_file), path(contigs_fasta)

    output:
    tuple val(sample_id),file("contigs.fasta"), file("${sample_id}.bam"),file("${sample_id}.bam.bai") emit: dorado_aligner_ch

    script:
    """
    ${params.dorado_software}/dorado aligner ${contigs_fasta} ${bam_file} -o .
    """
}
*/
process Dorado_Polisher {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dorado_polish/", mode: 'copy'
    label 'normal'
    clusterOptions = [ '--partition="normal"','--gpus=2' ]
    errorStrategy 'ignore'
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: dorado_polisher_ch

    script:
    """
    ${params.dorado_software}/dorado polish \\
        -x 'cuda:all' \\
        --bacteria \\
        -m ${params.model} \\
        -o . \\
        ${bam_file} \\
        ${contigs_fasta}
    """
}
