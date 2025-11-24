#!/bin/bash

# Directorio donde están tus archivos GFF3
DIRECTORIO="/home/nova/Desktop/alen-belen/results/genomad_bakta/plasmid"

# Directorio de salida
DIRECTORIO_SALIDA="/home/nova/Desktop/alen-belen/results/genomad_bakta/plasmid/output"

# Crear el directorio de salida si no existe
mkdir -p "$DIRECTORIO_SALIDA"

# Archivo de salida combinado
ARCHIVO_SALIDA="$DIRECTORIO_SALIDA/depth_genes_all_barcodes_output.txt"

# Limpiar el archivo de salida si ya existe
> "$ARCHIVO_SALIDA"

# Iterar sobre cada archivo que coincide con el patrón
for FILE in "$DIRECTORIO"/barcode*_plasmid_genes_bakta_mosdepth.gff3; do
    # Extraer el nombre base del archivo
    BASENAME=$(basename "$FILE")
    
    # Extraer el nombre del barcode (antes del primer guión bajo)
    BARCODE=$(echo "$BASENAME" | cut -d'_' -f1)
    
    # Procesar el archivo y añadir la columna del barcode
    grep -P -v $'\tregion\t' "$FILE" | awk -F '\t' -v barcode="$BARCODE" '{
        contig = $1;
        tipo = $3;
        inicio = $4;
        fin = $5;
        name = "NA";
        depth = "NA";
        split($9, atributos, ";");
        for (i in atributos) {
            n = split(atributos[i], kv, "=");
            if (n == 2) {
                key = kv[1];
                value = kv[2];
                if (key == "Name") {
                    name = value;
                } else if (key == "depth") {
                    depth = value;
                }
            }
        }
        print barcode "\t" contig "\t" tipo "\t" inicio "\t" fin "\t" name "\t" depth;
    }' >> "$ARCHIVO_SALIDA"
done
