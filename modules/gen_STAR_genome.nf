process genGenome {
    label 'process_high'
    publishDir "${projectDir}/modules/STAR/"

    script:
    """
    sh arriba/download_references.sh $params.arriba_genome_version

    STAR --runThreadN $task.cpus \
        --runMode genomeGenerate \
        --genomeDir $projectDir/STAR/STAR_index_GRCh38_GENCODE38 \
        --genomeFastaFiles $projectDir/arriba/GRCh38.fa\
        --sjdbGTFfile $projectDir/arriba/GENCODE38.gtf\
        --sjdbOverhang $params.R2_length - 1
    """
}
