#!/bin/bash

# Directorios base
BAKTA_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/bakta_2"              # Directorio base de bakta
GENOMAD_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/genomad_2"          # Directorio base de genomad
GENOMAD_BAKTA_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/genomad_bakta"  # Directorio base para los resultados

# Recorre cada barcode en el directorio BAKTA
for barcode_path in "$BAKTA_DIR"/*/; do
    # Extrae el nombre del barcode
    barcode=$(basename "$barcode_path")
    
    echo "Procesando barcode: $barcode"
    
    # Crea el directorio correspondiente en GENOMAD_BAKTA si no existe
    mkdir -p "$GENOMAD_BAKTA_DIR/$barcode/virus"

    # Ejecuta el primer comando (cut + grep + sort) usando rutas absolutas
    cut "$GENOMAD_DIR/$barcode/consensus_summary/consensus_virus_summary.tsv" -f1 | grep -v 'seq_name' | sort > "$GENOMAD_BAKTA_DIR/$barcode/virus/contigs_virus.txt"
    
    # Ejecuta el comando awk para filtrar el GFF3 sin cambiar de directorio
    awk '/##FASTA/ {exit} FNR==NR {contigs[$1]; next} $1 in contigs' "$GENOMAD_BAKTA_DIR/$barcode/virus/contigs_virus.txt" "$BAKTA_DIR/$barcode/consensus.gff3" > "$GENOMAD_BAKTA_DIR/$barcode/virus/virus_genes_bakta.gff3"
    
    # Copia el archivo plasmid/virus_genes_genomad.tsv al directorio correspondiente en GENOMAD_BAKTA
    cp "$GENOMAD_DIR/$barcode/consensus_summary/consensus_virus_genes.tsv" "$GENOMAD_BAKTA_DIR/$barcode/virus/virus_genes_genomad.tsv"

    echo "Procesamiento completado para $barcode"
done

echo "Script completado exitosamente."
