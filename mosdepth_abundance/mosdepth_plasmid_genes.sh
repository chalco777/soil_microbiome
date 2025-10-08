#!/bin/bash
# NO FUNCIONAeval "$(conda shell.bash hook)"  # This initializes conda for the current bash shell
# Activate mob_suit conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mosdepth

# Define rutas absolutas
GENOMAD_BAKTA_PATH="/home/nova/Desktop/alen-belen/results/genomad_bakta"
MINIMAP2_PATH="/home/nova/Desktop/alen-belen/results/minimap2"
OUTPUT_PATH="/home/nova/Desktop/alen-belen/results/mosdepth_mge_genes/plasmid"
LOG_PATH="/home/nova/Desktop/alen-belen/scripts/mge/mosdepth_logs"

# Crear carpetas de salida si no existen
mkdir -p "$OUTPUT_PATH"
mkdir -p "$LOG_PATH"

# Definir los ambientes de conda

# Iterar sobre cada subcarpeta en genomad_bakta
for barcode_dir in "$GENOMAD_BAKTA_PATH"/*; do
    barcode=$(basename "$barcode_dir")

    # Definir paths para el archivo GFF3 y BAM
    gff3_file="${GENOMAD_BAKTA_PATH}/${barcode}/plasmid/plasmid_genes_bakta.gff3"
    bam_file="${MINIMAP2_PATH}/${barcode}/${barcode}.bam"
    sorted_bam_file="${MINIMAP2_PATH}/${barcode}/${barcode}_sorted.bam"
    bed_file="${GENOMAD_BAKTA_PATH}/${barcode}/plasmid/plasmid_genes_bakta.bed"
    
    # Crear subcarpeta específica para el barcode en el output y logs
    mkdir -p "${OUTPUT_PATH}/${barcode}"
    log_file="${LOG_PATH}/${barcode}_process.log"
    mosdepth_output="${OUTPUT_PATH}/${barcode}/${barcode}_mosdepth.out"
    output_gff3="${GENOMAD_BAKTA_PATH}/${barcode}/plasmid/plasmid_genes_bakta_mosdepth.gff3"

    # Asegurarse de que los archivos GFF3 y BAM existen
    if [[ ! -f "$gff3_file" || ! -f "$bam_file" ]]; then
        echo "Error: Archivos faltantes para $barcode" >> "$LOG_PATH/error.log"
        continue
    fi

    # Crear el archivo BED desde el archivo GFF3
    echo "Generando archivo BED para $barcode" | tee -a "$log_file"
    grep -v 'region' "$gff3_file" | awk '{split($9, a, ";"); print $1"\t"$4"\t"$5"\t"a[1]}' > "$bed_file"
    wait
    if [[ ! -f "$bed_file" ]]; then
        echo "Error: No se pudo generar el archivo BED para $barcode" >> "$log_file"
        continue
    fi

    # Activar ambiente para samtools y ordenar el BAM
    echo "Ordenando BAM para $barcode" | tee -a "$log_file"
    samtools sort -@ 10 "$bam_file" -o "$sorted_bam_file" >> "$log_file" 2>&1
    wait
    samtools index "$sorted_bam_file" >> "$log_file" 2>&1
    wait
    
    # Verificar que el archivo BAM ordenado fue creado
    if [[ ! -f "$sorted_bam_file" ]]; then
        echo "Error: No se pudo ordenar el archivo BAM para $barcode" >> "$log_file"
        continue
    fi

    # Activar ambiente para mosdepth y calcular la profundidad
    echo "Calculando profundidad con mosdepth para $barcode" | tee -a "$log_file"
    mosdepth --by "$bed_file" "$mosdepth_output" "$sorted_bam_file" >> "$log_file" 2>&1
    wait

    # Verificar que el archivo de mosdepth fue generado
    if [[ ! -f "${mosdepth_output}.regions.bed.gz" ]]; then
        echo "Error: No se pudo generar el archivo de mosdepth para $barcode" >> "$log_file"
        continue
    fi

    # Añadir en el gff3 como nueva columna, la quinta columna generada por mosdepth que es la profundidad
    echo "Añadiendo la profundidad al archivo GFF3 para $barcode" | tee -a "$log_file"

    # Descomprimir el archivo de salida de mosdepth
    zcat "${mosdepth_output}.regions.bed.gz" > "${OUTPUT_PATH}/${barcode}/${barcode}_mosdepth.regions.bed"
    wait
    # Combinar el archivo GFF3 con la profundidad
    awk 'NR==FNR {split($4, id, "="); depth[id[2]]=$5; next} {split($9, a, ";"); split(a[1], id, "="); if (id[2] in depth) print $0 ";depth=" depth[id[2]]; else print $0}' "${OUTPUT_PATH}/${barcode}/${barcode}_mosdepth.regions.bed" "$gff3_file" > "$output_gff3"
    wait
    echo "Proceso completado para $barcode. Archivo GFF3 actualizado: $output_gff3" | tee -a "$log_file"

done

echo "Script finalizado. Revisa los logs en $LOG_PATH"
