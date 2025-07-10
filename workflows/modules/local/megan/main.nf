#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Megan 
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

process Meganize_BlastX_Contigs_Dragonflye {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_contigs/", mode: 'copy'
    label 'optimized_Meganize_Contigs_BlastX'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(daa_file)
    
    output:
    tuple val(sample_id), file("${sample_id}_dragonflye_contigs_blastx_diamondview.tsv"), emit: dragonflye_meganize_blastx_contigs_ch
    tuple val(sample_id), file("*"), emit: other_dragonflye_meganize_blastx_contigs_ch
    when:
    params.dragonflye

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${daa_file} \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --topPercent 0.5 \\
        --minSupportPercent 0 \\
        --lcaAlgorithm longReads \\
        --longReads true \\
        --verbose >> Meganize_BlastX_Contigs_Dragonflye.log 2>&1 

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \\
      -i ${daa_file} \\
      -c2c Taxonomy | awk '{print \$1}') \\
      <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${daa_file} \\
      -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_dragonflye_contigs_blastx_daa_summary_count.tsv

    diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
                 --daa ${daa_file} > ${sample_id}_dragonflye_contigs_blastx_diamondview.tsv
    """
}


process Meganize_BlastX_Contigs_Metaspades {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_contigs", mode: 'copy'
    label 'optimized_Meganize_Contigs_BlastX'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(daa_file)
    
    output:
    tuple val(sample_id), file("${sample_id}_metaspades_contigs_blastx_diamondview.tsv"), emit: metaspades_meganize_blastx_contigs_ch
    tuple val(sample_id), file("*"), emit: other_metaspades_meganize_blastx_contigs_ch
    tuple val(sample_id), file(daa_file), emit: meganized_blastx_contig_metaspades_file_ch
    when:
    params.metaspades

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${daa_file} \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --topPercent 0.5 \\
        --minSupportPercent 0 \\
        --lcaAlgorithm longReads \\
        --longReads true \\
        --verbose >> Meganize_BlastX_Contigs_Metaspades.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \\
      -i ${daa_file} \\
      -c2c Taxonomy | awk '{print \$1}') \\
      <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${daa_file} \\
      -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_metaspades_contigs_blastx_daa_summary_count.tsv

    diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
                 --daa ${daa_file} > ${sample_id}_metaspades_contigs_blastx_diamondview.tsv
    """
}


process Meganize_BlastX_Contigs_Unicycler {
    tag { sample_id }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/meganized_contigs", mode: 'copy'
    label 'optimized_Meganize_Contigs_BlastX'
    conda './env/blastx.yml'

    input:
    tuple val(sample_id), path(daa_file)
    
    output:
    tuple val(sample_id), file("${sample_id}_unicycler_contigs_blastx_diamondview.tsv"), emit: unicycler_meganize_blastx_contigs_ch
    tuple val(sample_id), file("*"), emit: other_unicycler_meganize_blastx_contigs_ch
    when:
    params.unicycler

    script:
    """
    eval "\$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate md

    ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
    ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

    ${params.meganpath}/daa-meganizer \\
        --in ${daa_file} \\
        --mapDB ${params.megandb}/megan-map.db \\
        --threads ${task.cpus} \\
        --topPercent 0.5 \\
        --minSupportPercent 0 \\
        --lcaAlgorithm longReads \\
        --longReads true \\
        --verbose >> Meganize_BlastX_Contigs_Unicycler.log 2>&1

    paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \\
      -i ${daa_file} \\
      -c2c Taxonomy | awk '{print \$1}') \\
      <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${daa_file} \\
      -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_unicycler_contigs_blastx_daa_summary_count.tsv

    diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
                 --daa ${daa_file} > ${sample_id}_unicycler_contigs_blastx_diamondview.tsv
    """
}


process Parse_BlastX_Dragonflye {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(diamondview_file)

    output:
    tuple val(sample_id), file("*"), emit: parse_blastx_dragonflye_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${diamondview_file} -t blastx -r species >> Parse_BlastX_Dragonflye.log 2>&1
    """
}

process Parse_BlastX_metaspades {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(diamondview_file)

    output:
    tuple val(sample_id), file("*"), emit: parse_blastx_metaspades_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${diamondview_file} -t blastx -r species >> Parse_BlastX_metaspades.log 2>&1
    """
}


process Parse_BlastX_unicycler {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'lowmem'

    input:
    tuple val(sample_id), file(diamondview_file)

    output:
    tuple val(sample_id), file("*"), emit: parse_blastx_unicycler_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${params.scripts}/VS_MD_diamond_parser_linFilt_working.py -i ${diamondview_file} -t blastx -r species >> Parse_BlastX_unicycler.log 2>&1
    """
}
