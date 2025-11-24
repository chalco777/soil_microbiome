#!/bin/bash

# Directorio base donde se encuentran las muestras
WORKDIR="/media/crowfoot2/DATOS/alen_belen_241109/"


BASE_DIR="$WORKDIR/samtools_extract_human_filtered"
SUBSAMPLING_BASE_DIR="$WORKDIR/Downsampling_1.6/fastq_concat"

cd /home/crowfoot2/bbmap
chmod +x reformat.sh

for SAMPLE_DIR in ${BASE_DIR}/barcode33; do
    # Extraer el nombre de la muestra (por ejemplo, barcode01) del directorio
    SAMPLE_NAME=$(basename $SAMPLE_DIR)
    
    # Ruta del archivo fastq.gz de entrada
    FASTQ_PATH="${SAMPLE_DIR}/${SAMPLE_NAME}.fastq.gz"
  
    # Verificar si el archivo fastq existe
    if [[ -f $FASTQ_PATH ]]; then
        # Crear el directorio de salida para la muestra si no existe
        OUTPUT_DIR="${SUBSAMPLING_BASE_DIR}/${SAMPLE_NAME}"
        mkdir -p $OUTPUT_DIR
        
        # Ejecutar reformat.sh para la muestra
        ./reformat.sh in=$FASTQ_PATH out=${OUTPUT_DIR}/${SAMPLE_NAME}.fastq.gz samplebasestarget=1600000000
    else
        echo "ADVERTENCIA: ${FASTQ_PATH} no existe. Omitiendo ${SAMPLE_NAME}..."
    fi
done

echo "Procesamiento de subsampling completado para todas las muestras."
