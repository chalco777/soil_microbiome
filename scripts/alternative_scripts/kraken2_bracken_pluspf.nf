#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "/media/crowfoot2/DATOS/alen_belen_241109/Downsampling_1.6/fastq_concat" 
params.kraken2DB = "/media/crowfoot2/DATOS/alen_belen_241109/databases/kraken/aria2c_PlusPF" 

process Kraken2 {
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1'
    publishDir '/media/crowfoot2/DATOS/alen_belen_241109/Downsampling_1.6/krakenplusPF/kraken2'
    cpus 15
    maxForks 1

    input:
    tuple val(sample_id), path(query_file)

    output:
    val sample_id
    path "${sample_id}/${sample_id}_kraken2_report.tsv"
    path "${sample_id}/${sample_id}_kraken2_classification.tsv"

    script:
    """
    mkdir -p ${sample_id}

    kraken2 --db $params.kraken2DB \\
    --threads ${task.cpus} \\
    --report "${sample_id}/${sample_id}_kraken2_report.tsv" \\
    --output "${sample_id}/${sample_id}_kraken2_classification.tsv" \\
    $query_file

    """
}

process Bracken {
    container 'quay.io/biocontainers/bracken:2.9--py311h2a4ad6c_1'
    publishDir '/media/crowfoot2/DATOS/alen_belen_241109/Downsampling_1.6/krakenplusPF/bracken', mode: 'copy' 
    cpus 15
    maxForks 1

    input:
    val sample_id  
    path query_file
    path classification
    path kraken2DB
    
    output:
    val sample_id
    path "${sample_id}_bracken"

    script:
    """
    bracken \\
    -d $kraken2DB \\
    -i $query_file \\
    -o ${sample_id}_bracken \\
    -r 300 \\
    -l S
    """
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/barcode*/*.fastq", type:"file",size:1)
    //channel_fastq.view()
    Kraken2(channel_fastq)
    Bracken(Kraken2.out, params.kraken2DB)
}

