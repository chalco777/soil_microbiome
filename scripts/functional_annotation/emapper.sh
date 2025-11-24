#!/bin/bash

# Ruta base donde están los directorios de los barcodes
BASE_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/bakta_2"

# Opciones comunes de eggNOG-mapper
EMAPPER_OPTS="-m diamond --cpu 10 --report_orthologs --excel --dmnd_iterate yes --go_evidence all"

# Iterar sobre cada directorio de barcode
for BARCODE_DIR in "$BASE_DIR"/barcode*; do
    echo "Procesando $BARCODE_DIR"

    # Definir rutas absolutas para los archivos
    GFF_FILE="$BARCODE_DIR/consensus.gff3"
    FAA_FILE="$BARCODE_DIR/consensus.faa"
    OUTPUT_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/eggnog/$(basename "$BARCODE_DIR")"

    # Crear el directorio de salida si no existe
    mkdir -p "$OUTPUT_DIR"

    # Procesar consensus.gff3
    PROCESSED_GFF="$OUTPUT_DIR/consensus_processed.gff3"
    # Eliminar todo a partir de la línea que contiene '##FASTA'
    sed '/^##FASTA/,$d' "$GFF_FILE" | \
    # Eliminar líneas que comienzan con '#'
    sed '/^#/d' > "$PROCESSED_GFF"

    # Procesar consensus.faa
    PROCESSED_FAA="$OUTPUT_DIR/consensus_processed.faa"
    awk '/^>/ {print $1; next} {print}' "$FAA_FILE" > "$PROCESSED_FAA"

    # Ejecutar eggNOG-mapper
    OUTPUT_PREFIX="$OUTPUT_DIR/resultados_emapper_$(basename "$BARCODE_DIR")"
    emapper.py -i "$PROCESSED_FAA" -o "$OUTPUT_PREFIX" $EMAPPER_OPTS \
    --output_dir "$OUTPUT_DIR" \
    --decorate_gff "$PROCESSED_GFF"

done
