nextflow run  /export/nextflow/bdrd/all_in_one_pipeline/md_vs_workflow.nf --samplesheet_csv=/export/nextflow/bdrd/all_in_one_pipeline/samplesheet.csv \
                          --outdir=/export/nextflow/bdrd/all_in_one_pipeline/results_test2/ \
                          --fastqc_threads=100 \
                          --project_id="TEST_DawnBats" \
                          -profile cluster_normal 
