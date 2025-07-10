#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   CheckM Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process CheckM_Dragonflye {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/checkm/assemblies/dragonflye/", mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: checkm_dragonflye
    when:
    params.dragonflye

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate checkm-genome
    checkm lineage_wf -t 10 -x fa --tab_table -f checkm_summary.tsv . . --tmpdir /export/tmp
    checkm qa -t 10 --tab_table -f checkm_full_stats.tsv -o 2 lineage.ms .
    """
}

process CheckM_Dragonflye_Medaka {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/checkm/assemblies/dragonflye_medaka/", mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: checkm_dragonflye_medaka
    when:
    params.dragonflye || params.medaka

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate checkm-genome
    checkm lineage_wf -t 10 -x fa --tab_table -f checkm_summary.tsv . . --tmpdir /export/tmp
    checkm qa -t 10 --tab_table -f checkm_full_stats.tsv -o 2 lineage.ms .
    """
}

process CheckM_Unicycler {
    tag { sample_id }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/checkm/assemblies/unicycler/", mode: 'copy'
    label 'normal'
    errorStrategy 'ignore'
    cpus { 32 }
    memory { '64 GB'}
    time '36h'

    input:
    tuple val(sample_id), file(assembly_fasta)

    output:
    tuple val(sample_id), file("*"), emit: checkm_unicycler
    when:
    params.unicycler

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate checkm-genome
    checkm lineage_wf -t 10 -x fasta --tab_table -f checkm_summary.tsv . . --tmpdir /export/tmp
    checkm qa -t 10 --tab_table -f checkm_full_stats.tsv -o 2 lineage.ms .
    """
}
