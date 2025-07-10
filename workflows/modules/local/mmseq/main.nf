#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   MMSEQ
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process MMSEQ_contigs_against_NT_Dragonflye {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_MMSEQ_contigs_against_NT_metaspades'
    //cpus { 128 }
    //memory { '900 GB'}
    //time '48h'
    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id),file("${sample_id}_mmseq_dragonflye_contig_against_NT.out") ,emit: mmseqs_dragonflye_contigs_ch
    
    when:
    params.dragonflye

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    mmseqs easy-search ${contigs} \
      ${params.db_BN} \
      ${sample_id}_mmseq_dragonflye_contig_against_NT.out \
      ${params.mmseq_tmp_dir} \
      --max-accept 25 \
      --threads ${task.cpus} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen >> MMSEQ_contigs_against_NT_Dragonflye.log 2>&1
    """
}

process MMSEQ_contigs_against_NT_metaspades {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_MMSEQ_contigs_against_NT_metaspades'
    //cpus { 128 }
    //memory { '900 GB'}
    //time '48h'
    conda './env/vs.yml'  // ← Use a conda environment here instead of eval manually

    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id), file("${sample_id}_mmseq_metaspades_contig_against_NT.out"), emit: mmseqs_metaspades_contigs_ch

    when:
    params.metaspades

    script:
    """
    mmseqs easy-search ${contigs} \
        ${params.db_BN} \
        ${sample_id}_mmseq_metaspades_contig_against_NT.out \
        ${params.mmseq_tmp_dir} \
        --max-accept 25 \
        --threads ${task.cpus} \
        -s 7.0 -e 1.0E-8 --search-type 3 \
        --split-memory-limit 120G \
        --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen >> MMSEQ_contigs_against_NT_metaspades.log 2>&1
    """
}

process MMSEQ_contigs_against_NT_unicycler {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_MMSEQ_contigs_against_NT_metaspades'
    //cpus { 128 }
    //memory { '900 GB'}
    time '48h'
    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id),file("${sample_id}_mmseq_unicycler_contig_against_NT.out") ,emit: mmseqs_unicycler_contigs_ch
    
    when:
    params.unicycler

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    mmseqs easy-search ${contigs} \
      ${params.db_BN} \
      ${sample_id}_mmseq_unicycler_contig_against_NT.out \
      ${params.mmseq_tmp_dir} \
      --max-accept 25 \
      --threads ${task.cpus} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen >> MMSEQ_contigs_against_NT_unicycler.log 2>&1
    """
}


process Parse_MMSEQ_dragonflye_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id),file("*") ,emit: mmseqs_dragonflye_contigs_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${sample_id}_mmseq_dragonflye_contig_against_NT.out -t mmseqs -r species >> Parse_MMSEQ_dragonflye_contigs.log 2>&1
    """
}

process Parse_MMSEQ_metaspades_contigs {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id),file("*") ,emit: mmseqs_metaspades_contigs_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${sample_id}_mmseq_metaspades_contig_against_NT.out -t mmseqs -r species >> Parse_MMSEQ_metaspades_contigs.log 2>&1
    """
}

process Parse_MMSEQ_unicycler_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(contigs)

    output:
    tuple val(sample_id),file("*") ,emit: mmseqs_unicycler_contigs_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${sample_id}_mmseq_unicycler_contig_against_NT.out -t mmseqs -r species >> Parse_MMSEQ_unicycler_contigs.log 2>&1
    """
}
