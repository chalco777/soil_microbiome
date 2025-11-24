#!/bin/bash

# Directorios base
BED_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/eggnog"
ALIGNMENT_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/minimap_result"
GENOME_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/bakta_2"
OUTPUT_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/count_bedtools_3"

# Crear el directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Iterar sobre cada directorio de barcode en BED_DIR
for barcode_dir in "$BED_DIR"/barcode*; do
    barcode=$(basename "$barcode_dir")  # e.g., barcode03

    # Buscar el archivo GFF decorado en el directorio del barcode
    gff_file="$barcode_dir/resultados_emapper_${barcode}.emapper.decorated.gff"
    if [[ ! -f "$gff_file" ]]; then
        echo "Archivo GFF no encontrado para $barcode en $barcode_dir. Saltando..."
        continue
    fi

    # Ordenar el archivo GFF y guardar en OUTPUT_DIR
    sorted_gff_file="$OUTPUT_DIR/${barcode}_decorated.sorted.gff"
        grep -v "^#" "$gff_file" | grep -v -P '\tregion\t' > "$sorted_gff_file"
    echo "Archivo BED ordenado: $sorted_bed_file"

    # Convertir el GFF a formato BED, usando el ID como nombre
    bed_file="$OUTPUT_DIR/${barcode}_decorated.bed"
    awk 'BEGIN {FS=OFS="\t"} {
        # Extraer el ID del campo de atributos (columna 9)
        match($9, /ID=([^;]+)/, arr);
        id = arr[1];
        # Imprimir en formato BED
        print $1, $4-1, $5, id
    }' "$sorted_gff_file" > "$bed_file"

    sorted_bed_file="$OUTPUT_DIR/${barcode}_sorted.bed"
    sort -k1.8,1n -k2,2n "$bed_file" > "$sorted_bed_file"

    echo "Archivo BED ordenado: $sorted_bed_file"

    # Buscar el archivo BAM correspondiente
    bam_dir="$ALIGNMENT_DIR/$barcode"
    bam_file="$bam_dir/${barcode}_sorted.bam"

    # Excluir lecturas no mapeadas y alineamientos suplementarios
    filtered_bam_file="$bam_dir/${barcode}_filtered.bam"
    samtools view -F 2308 -b "$bam_file" > "$filtered_bam_file"
    #samtools sort -o ${bam_file%.bam}_sorted.bam "$bam_file"
    #mv ${bam_file%.bam}_sorted.bam "$bam_file"
    # Verificar si el archivo BAM existe antes de ejecutar bedtools
    if [[ -f "$filtered_bam_file" ]]; then
        # Ruta al archivo de genoma (consensus.fna)
        genome_dir="$GENOME_DIR/$barcode"
        genome_fasta="$genome_dir/consensus.fna"

        # Verificar si el archivo de genoma existe
        if [[ -f "$genome_fasta" ]]; then
            # Crear índice .fai si no existe
            fai_file="$genome_fasta.fai"
            if [[ ! -f "$fai_file" ]]; then
                samtools faidx "$genome_fasta"
                echo "Archivo FAI creado: $fai_file"
            fi

            # Crear el genome file para bedtools
            genome_file="$OUTPUT_DIR/${barcode}_genome.txt"
            awk '{print $1"\t"$2}' "$fai_file" > "$genome_file"
            echo "Genome file creado: $genome_file"

            # Ejecutar bedtools coverage y guardar el resultado en la carpeta de salida
            coverage_output="$OUTPUT_DIR/${barcode}_coverage.txt"
            bedtools coverage -a "$sorted_bed_file" -b "$filtered_bam_file" -g "$genome_file" -f 0.5 -sorted > "$coverage_output"

            # Añadir las columnas de cobertura al archivo GFF original basado en el ID
            gff_with_coverage="$OUTPUT_DIR/${barcode}_decorated_with_coverage.gff"

            # Crear un archivo temporal con los datos de cobertura y el ID como clave
            coverage_data="$OUTPUT_DIR/${barcode}_coverage_data.txt"
            awk 'BEGIN {FS=OFS="\t"} {print $4, $0}' "$coverage_output" > "$coverage_data"

            # Procesar el archivo GFF y añadir los datos de cobertura
            awk -v cov_file="$coverage_data" 'BEGIN {
                FS=OFS="\t";
                # Cargar los datos de cobertura en un array asociativo
                while ((getline < cov_file) > 0) {
                    id = $1;
                    coverage[id] = $0;
                }
                close(cov_file);
            }
            {
                # Extraer el ID del campo de atributos
                match($9, /ID=([^;]+)/, arr);
                id = arr[1];

                # Si el ID existe en los datos de cobertura, añadir los atributos
                if (id in coverage) {
                    split(coverage[id], cov_fields, FS);
                    conteo = cov_fields[6];
                    coverage_bases = cov_fields[7];
                    longitud_feature = cov_fields[8];
                    coverage_fraction = cov_fields[9];

                    # Añadir los nuevos atributos al campo 9
                    if ($9 == ".") $9 = "";
                    $9 = $9 ";CONTEO=" conteo ";coverage_bases=" coverage_bases ";longitud_feature=" longitud_feature ";coverage_fraction=" coverage_fraction;
                }
                print
            }' "$sorted_gff_file" > "$gff_with_coverage"
            echo "Archivo GFF con cobertura añadida: $gff_with_coverage"

            # Eliminar archivo temporal
            rm "$coverage_data"

        else
            echo "Archivo de genoma no encontrado: $genome_fasta"   
        fi
    else
        echo "Archivo BAM no encontrado: $bam_file"
    fi
done