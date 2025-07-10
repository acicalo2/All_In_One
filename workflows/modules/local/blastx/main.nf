#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   BlastX
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process BlastX_Dragonflye_contigs {
    tag { sample_id }
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_contigs'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(contigs_fasta)

    output:
    tuple val(sample_id), file("${sample_id}_dragonflye_blastx.daa"), emit: dragonflye_blastx_contigs_ch
    when: 
    params.dragonflye
    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md

    diamond blastx ${params.diamond_args} \\
        --threads ${task.cpus} \\
        --db ${params.diamond_dbdir} \\
        --query ${contigs_fasta} \\
        --range-culling -F 15 \\
        --evalue 1e-5 \\
        --outfmt 100 \\
        --out ${sample_id}_dragonflye_blastx.daa >> BlastX_Dragonflye_contigs.log 2>&1

    echo "blastx complete: ${contigs_fasta}" > ${sample_id}_dragonflye_blastx.finished
    """
}

process BlastX_metaspades_contigs {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_contigs'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(contigs_fasta)

    output:
    tuple val(sample_id), file("${sample_id}_metaspades_blastx.daa"), emit: metaspades_blastx_contigs_ch
    
    when: 
    params.metaspades

    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md

    diamond blastx ${params.diamond_args} \\
        --threads ${task.cpus} \\
        --db ${params.diamond_dbdir} \\
        --query ${contigs_fasta} \\
        --range-culling -F 15 \\
        --evalue 1e-5 \\
        --outfmt 100 \\
        --out ${sample_id}_metaspades_blastx.daa >> BlastX_metaspades_contigs.log 2>&1

    echo "blastx complete: ${contigs_fasta}" > ${sample_id}_metaspades_blastx.finished
    """
}

process BlastX_spades_contigs {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_contigs'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(contigs_fasta)

    output:
    tuple val(sample_id), file("${sample_id}_spades_blastx.daa"), emit: spades_blastx_contigs_ch
    
    when: 
    params.spades

    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md

    diamond blastx ${params.diamond_args} \\
        --threads ${task.cpus} \\
        --db ${params.diamond_dbdir} \\
        --query ${contigs_fasta} \\
        --range-culling -F 15 \\
        --evalue 1e-5 \\
        --outfmt 100 \\
        --out ${sample_id}_spades_blastx.daa >> BlastX_spades_contigs.log 2>&1

    echo "blastx complete: ${contigs_fasta}" > ${sample_id}_spades_blastx.finished
    """
}

process BlastX_unicycler_contigs {
    tag { sample_id }
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_contigs'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(contigs_fasta)

    output:
    tuple val(sample_id), file("${sample_id}_unicycler_blastx.daa"), emit: unicycler_blastx_contigs_ch
    
    when: 
    params.spades

    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md

    diamond blastx ${params.diamond_args} \\
        --threads ${task.cpus} \\
        --db ${params.diamond_dbdir} \\
        --query ${contigs_fasta} \\
        --range-culling -F 15 \\
        --evalue 1e-5 \\
        --outfmt 100 \\
        --out ${sample_id}_unicycler_blastx.daa >> BlastX_unicycler_contigs.log 2>&1

    echo "blastx complete: ${contigs_fasta}" > ${sample_id}_unicycler_blastx.finished
    """
}

process BlastX_reads {
    tag { "${sample_id}_${mode}" }
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'optimized_blastx_reads'
    conda './env/diamond.yml'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(fastq_file), val(mode)

    output:
    tuple val(sample_id), file("*"), emit: blastx_reads_ch
    tuple val(sample_id), file("${sample_id}_short_reads_blastx.daa"), optional: true, emit: short_reads_blastx_channel
    tuple val(sample_id), file("${sample_id}_long_reads_blastx.daa"), optional: true, emit: long_reads_blastx_channel
    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    mkdir -p logs

    if [[ "$mode" == "long" && -s "$fastq_file" ]]; then
        echo "Running diamond on long read: $fastq_file" | tee logs/blastx_long.log
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "$fastq_file" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_long_reads_blastx.daa" >> BlastX_longreads.log 2>&1
    elif [[ "$mode" == "short" && -s "$fastq_file" ]]; then 
        echo "Running diamond on short reads: $fastq_file" | tee logs/blastx_short.log
        diamond blastx ${params.diamond_args} \\
            --threads ${task.cpus} \\
            --db ${params.diamond_dbdir} \\
            --query "$fastq_file" \\
            --evalue 1e-5 \\
            --outfmt 100 \\
            --out "${sample_id}_short_reads_blastx.daa" >> BlastX_shortreads.log 2>&1
    else
        echo "SKIP: No valid reads provided for $mode mode" > "${sample_id}_blastx.skipped.txt"
    fi

    echo "OK" > "${sample_id}_blastx_reads.finished"
    """
}

// DAA2INFO contigs .daa file

process DAA2INFO_contigs_daa_file {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'

    input:

    tuple val(sample_id), file(blastx_contigs)

    output:

    tuple val(sample_id), file("${sample_id}_c2c.txt"), emit: c2c_txt_file


    script:

    """
    if [ "${params.metaspades}" == "true" ]; then
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_metaspades_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u

    elif [ "${params.hybrid}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.dragonflye_isolate}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.unicycler}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    else 
      "skip this process"
  
    fi
    """
}

/* Meganize Blastx Reads */

process Meganize_ShortReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_reads", mode: 'copy'
    label 'optimized_blastx_reads'
    conda './env/md.yml'
    cpus { 64 }
    memory { '260 GB'}
    time '36h'

    input:
    tuple val(sample_id), file(blastx_short)

    output:
    tuple val(sample_id), file("${sample_id}_short_reads_blastx_daa_summary_count.tsv"), emit: meganize_short_reads_ch
    tuple val(sample_id), file("*"), emit: all_shortreads_meganized_files
    tuple val(sample_id), file(blastx_short), emit: meganized_short_daa_ch
    when:
    params.shortreads || params.hybrid

    script:
    """
    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${blastx_short} \\
        --only Taxonomy \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads false \\
        --verbose >> meganize_short_reads.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_short} -c2c Taxonomy | awk '{print \$1}') \\
          <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_short} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \\
          > ${sample_id}_short_reads_blastx_daa_summary_count.tsv
    """
}

process Meganize_LongReads_BlastX {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_reads", mode: 'copy'
    label 'optimized_blastx_reads'
    conda './env/md.yml'
    cpus { 64 }
    memory { '260 GB'}
    time '36h'

    input:
    tuple val(sample_id), file(blastx_long)

    output:
    tuple val(sample_id), file("${sample_id}_long_reads_blastx_daa_summary_count.tsv"), emit: meganize_long_reads_ch
    tuple val(sample_id), file("*"), emit: all_longreads_meganized_files
    tuple val(sample_id), file(blastx_long), emit: long_reads_daa_meganized_ch
    when:
    params.longreads || params.hybrid

    script:
    """
    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${blastx_long} \\
        --only Taxonomy \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --minSupportPercent 0 \\
        --topPercent 0.5 \\
        --lcaAlgorithm weighted \\
        --longReads false \\
        --verbose >> meganize_longreads.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_long} -c2c Taxonomy | awk '{print \$1}') \\
          <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${blastx_long} -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') \\
          > ${sample_id}_long_reads_blastx_daa_summary_count.tsv
    """
}

