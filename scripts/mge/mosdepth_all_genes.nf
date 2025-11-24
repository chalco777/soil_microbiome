#!/usr/bin/env nextflow
//no funciono porq mosdepth esta mal
//arreglar los publish dir noma
params.fastqDir="/media/nova/datos/alen-belen-2/alen-belen/results_reads_filtered/samtools_extract"
params.genome="/media/nova/datos/alen-belen-2/alen-belen/results/bakta"

params.GENOMAD_BAKTA_PATH = "/media/nova/datos/alen-belen-2/alen-belen/results/bakta"
params.MINIMAP2_PATH = "/media/nova/datos/alen-belen-2/alen-belen/results/minimap2_bakta"
params.OUTPUT_PATH = "/media/nova/datos/alen-belen-2/alen-belen/results/mosdepth_all_genes"
params.LOG_PATH = "/media/nova/datos/alen-belen-2/alen-belen/scripts/mge/mosdepth_all_genes_logs"

process SortIndexBAM {
    conda '/home/nova/anaconda3/envs/minimap2'
    cpus 20
    maxForks 2
    input:
    tuple val(barcode), path(reads)

    output:
    path("${barcode}_sorted.bam")
    path("${barcode}_sorted.bam.bai")
    
    publishDir "${params.MINIMAP2_PATH}/${barcode}", mode: 'copy'

    script:
    """
    echo "Aligning and Sorting BAM for ${barcode}"
    minimap2 -ax map-ont -t ${task.cpus} ${params.genome}/${barcode}/consensus.fna ${reads}/${barcode}_filtered.fastq.gz | \
    samtools sort -@ ${task.cpus} -o ${barcode}_sorted.bam -
    samtools index ${barcode}_sorted.bam
    """
}

process GenerateSampleList {

    output:
    path 'samples_list.txt', emit: samples_list_ch

    script:
    """
    GENOMAD_BAKTA_PATH="${params.GENOMAD_BAKTA_PATH}"
    MINIMAP2_PATH="${params.MINIMAP2_PATH}"
    mkdir -p ${params.LOG_PATH}
    for barcode_dir in "\${GENOMAD_BAKTA_PATH}"/*; do
        barcode=\$(basename "\$barcode_dir")
        # Only process 'barcode01', QUITAR EL IF PARA CORRER LOS DEMAS
        gff3_file="\${GENOMAD_BAKTA_PATH}/\${barcode}/consensus.gff3"
        bam_file="\${MINIMAP2_PATH}/\${barcode}/\${barcode}_sorted.bam"
        bai_file="\${MINIMAP2_PATH}/\${barcode}/\${barcode}_sorted.bam.bai"
        if [[ -f "\$gff3_file" && -f "\$bam_file" ]]; then
            echo "\$barcode \$gff3_file \$bam_file \$bai_file"
        else
            echo "Error: Missing files for \$barcode" >> "${params.LOG_PATH}/error.log"
        fi
    done > samples_list.txt
    """
}

process GenerateBED {
    input:
    tuple val(barcode), path(gff3_file), path(bam_file), path(bai_file)

    output:
    tuple val(barcode), path("${barcode}_all_genes_bakta.bed"), path(gff3_file), path(bam_file), path(bai_file)
    
    script:
    """
    echo "Generating BED file for ${barcode}"
    sed '/##FASTA/,\$d' ${gff3_file} | grep -P -v '\tregion\t' | grep -v '^#' | awk '{split(\$9, a, ";"); split(a[1], b, "="); print \$1"\t"\$4"\t"\$5"\t"b[2]}' > ${barcode}_all_genes_bakta.bed
    """
    //sed elimina (d es de delete)desde fasta hasta el final del archivo ($)
}

process RunMosdepth {
    conda '/home/nova/anaconda3/envs/mosdepth'
    cpus 30
    maxForks 2
    input:
    tuple val(barcode), path(bed_file), path(gff3_file), path(bam_file), path(bai_file)

    output:
    path("${barcode}_mosdepth.regions.bed.gz")
    tuple val(barcode), path(bed_file), path(gff3_file), path(bam_file), path(bai_file)

    publishDir "${params.OUTPUT_PATH}/${barcode}/", pattern: "*_mosdepth.regions.bed.gz", mode: 'copy'

    script:
    """
    echo "Running mosdepth for ${barcode}"
    mosdepth -F 2304 -t ${task.cpus} --by ${bed_file} ${barcode}_mosdepth ${bam_file}
    """
}

process CombineDepthGFF3 {
    input:
    path(mosdepth)
    tuple val(barcode), path(bed_file), path(gff3_file), path(bam_file), path(bai_file)
    output:
    path("${barcode}_all_genes_bakta_mosdepth.gff3")
    path("${barcode}_mosdepth.regions.bed")
    path(bed_file)
    
    publishDir "${params.GENOMAD_BAKTA_PATH}/${barcode}", pattern: "*.gff3", mode:'copy'
    publishDir "${params.OUTPUT_PATH}/${barcode}", pattern: "*regions*", mode: 'copy'
    publishDir "${params.GENOMAD_BAKTA_PATH}/${barcode}", pattern: "*bakta.bed", mode: 'copy'

    script:
    regions_bed = "${barcode}_mosdepth.regions.bed"
    output_gff3 = "${barcode}_all_genes_bakta_mosdepth.gff3"
    """
    echo "Combining depth into GFF3 for ${barcode}"
    zcat ${mosdepth}> ${regions_bed}
    awk 'NR==FNR {
            depth[\$4] = \$5;
            next
        }
        {
            split(\$9, a, ";");
            split(a[1], id, "=");
            if (id[2] in depth)
                print \$0 ";depth=" depth[id[2]];
            else
                print \$0
        }' ${regions_bed} ${gff3_file} > ${output_gff3}
    """
}
workflow {
    //channel_reads=channel.fromFilePairs("$params.fastqDir/*", type:'dir',size:1)
    //SortIndexBAM(channel_reads)
    //def samples_list_ch = GenerateSampleList().out
    //samples_list_ch.view()
    GenerateSampleList()
    GenerateSampleList.out.samples_list_ch.flatMap { file -> file.readLines() }.map { line ->
            def (barcode, gff3_file, bam_file, bai_file) = line.tokenize()
            tuple(barcode, file(gff3_file), file(bam_file), file(bai_file))
        }.set { samples_ch }
    samples_ch.view()
    //Continue with the pipeline
    GenerateBED(samples_ch)
    RunMosdepth(GenerateBED.out)
    //OJO: ACÁ PRIMERO CORRER RUN MOSDEPTH Y LUEGO COMBINE DEPTH EN OTRO,
    //NO ENCONTRÉ LA FORMA DE QUE COMBINE DEPTH CORRA DEPUES DE QUE TODO RUNMOSDEPTH TERMINE, POR LOS PROCESOS SON INDEPENDIENTES SI ES QUE NO SE CONECTAN POR OUTPUT

    CombineDepthGFF3(RunMosdepth.out)

}