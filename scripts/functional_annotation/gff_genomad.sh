#!/bin/bash

# Definir las rutas absolutas de los directorios
BASE_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2"
RESULTS_DIR="$BASE_DIR/results"
GFF_DIR="$RESULTS_DIR/count_bedtools_3"
CONSENSUS_DIR="$RESULTS_DIR/genomad_2"

# Directorio de salida para los archivos GFF modificados y secuencias no asignadas
OUTPUT_DIR="$RESULTS_DIR/gff_final_3"
UNASSIGNED_DIR="$RESULTS_DIR/unassigned_sequences"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$UNASSIGNED_DIR"

# Iterar sobre todos los barcodes existentes
for BARCODE in $(ls "$CONSENSUS_DIR")
do
    # Verificar que el directorio es un barcode válido
    if [[ $BARCODE == barcode* ]]; then
        echo "Procesando $BARCODE..."

        # Definir las rutas absolutas para los archivos
        GFF_FILE="$GFF_DIR/${BARCODE}_decorated_with_coverage.gff"
        VIRUS_SUMMARY="$CONSENSUS_DIR/$BARCODE/consensus_summary/consensus_virus_summary.tsv"
        PLASMID_SUMMARY="$CONSENSUS_DIR/$BARCODE/consensus_summary/consensus_plasmid_summary.tsv"

        # Verificar si el archivo GFF existe
        if [ ! -f "$GFF_FILE" ]; then
            echo "El archivo GFF para $BARCODE no existe: $GFF_FILE"
            continue
        fi

        # Verificar si los archivos de resumen existen
        if [ ! -f "$VIRUS_SUMMARY" ]; then
            echo "El archivo consensus_virus_summary.tsv para $BARCODE no existe: $VIRUS_SUMMARY"
            continue
        fi
        if [ ! -f "$PLASMID_SUMMARY" ]; then
            echo "El archivo consensus_plasmid_summary.tsv para $BARCODE no existe: $PLASMID_SUMMARY"
            continue
        fi

        # Archivo de salida
        OUTPUT_GFF="$OUTPUT_DIR/${BARCODE}_decorated_with_coverage_mge.gff"
        UNASSIGNED_FILE="$UNASSIGNED_DIR/${BARCODE}_unassigned_sequences.txt"

        # Procesar el archivo GFF y añadir la nueva columna
        awk -v VS="$VIRUS_SUMMARY" -v PS="$PLASMID_SUMMARY" -v UNASSIGNED="$UNASSIGNED_FILE" '
        BEGIN {
            FS=OFS="\t"

            # Leer el archivo GFF para obtener los seq_names presentes
            while ((getline < "'"$GFF_FILE"'") > 0) {
                if ($1 ~ /^#/) continue  # Saltar líneas de comentarios
                gff_seq_names[$1] = 1
            }
            close("'"$GFF_FILE"'")

            # Leer el archivo de resumen de virus y crear una tabla hash
            while ((getline < VS) > 0) {
                if (NR_VS == 0) { NR_VS++; continue }  # Saltar el encabezado
                seq_name_full = $1
                seq_length = $2
                topology = $3
                coordinates = $4

                # Extraer el nombre base de la secuencia (antes de "|")
                split(seq_name_full, arr, "\\|")
                seq_name = arr[1]

                # Si es un provirus, almacenamos también las coordenadas
                if (topology == "Provirus") {
                    virus_info[seq_name] = virus_info[seq_name] ? virus_info[seq_name] ";" "VIRUS_PROVIRUS;" seq_length ";" coordinates : "VIRUS_PROVIRUS;" seq_length ";" coordinates
                } else {
                    virus_info[seq_name] = virus_info[seq_name] ? virus_info[seq_name] ";" "VIRUS;" seq_length : "VIRUS;" seq_length
                }

                # Marcar que este seq_name está en el resumen de virus
                virus_seq_names[seq_name_full] = seq_name
            }
            close(VS)

            # Leer el archivo de resumen de plásmidos y crear una tabla hash
            while ((getline < PS) > 0) {
                if (NR_PS == 0) { NR_PS++; continue }  # Saltar el encabezado
                seq_name = $1
                seq_length = $2
                plasmid_info[seq_name] = "PLASMID;" seq_length

                # Marcar que este seq_name está en el resumen de plásmidos
                plasmid_seq_names[seq_name] = 1
            }
            close(PS)
        }

        {
            # Procesar el archivo GFF
            if ($1 ~ /^#/) {
                print $0, "MGE_TYPE"  # Añadir encabezado si es necesario
                next
            }

            seq_name = $1
            new_column = ""

            if (seq_name in virus_info) {
                new_column = virus_info[seq_name]
                assigned_virus[seq_name] = 1
            } else if (seq_name in plasmid_info) {
                new_column = plasmid_info[seq_name]
                assigned_plasmid[seq_name] = 1
            } else {
                new_column = "NOT_MGE"
            }

            print $0, new_column
        }
        END {
            # Guardar las secuencias de virus no asignadas
            for (seq_name_full in virus_seq_names) {
                seq_name = virus_seq_names[seq_name_full]
                if (!(seq_name in assigned_virus)) {
                    print seq_name_full > UNASSIGNED
                }
            }
            # Guardar las secuencias de plásmidos no asignadas
            for (seq_name in plasmid_seq_names) {
                if (!(seq_name in assigned_plasmid)) {
                    print seq_name > UNASSIGNED
                }
            }
        }
        ' "$GFF_FILE" > "$OUTPUT_GFF"

        echo "Procesado $BARCODE: archivo modificado guardado en $OUTPUT_GFF"
        echo "Secuencias no asignadas guardadas en $UNASSIGNED_FILE"
    fi
done
