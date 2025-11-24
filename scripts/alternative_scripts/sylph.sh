#!/bin/bash

# Activar el entorno de conda para Sylph
echo "Activating conda environment for Sylph..."
source ~/anaconda3/etc/profile.d/conda.sh
conda activate sylph || { echo "ERROR: Failed to activate Conda environment"; exit 1; }

# Limpiar cualquier variable GROUPS preexistente
unset GROUPS

# Definir el directorio de trabajo y el número de hilos
WORK_DIR="/media/crowfoot2/DATOS/alen_belen_241109"  # Ajusta esta ruta según sea necesario
THREADS=15

# Moverse al directorio de trabajo
cd "$WORK_DIR" || { echo "ERROR: Failed to access working directory $WORK_DIR"; exit 1; }

# Definir la lista de grupos
GROUPS=("samtools_extract_human_filtered")  # Añade los nombres de tus grupos aquí

# Definir la ruta a la base de datos GTDB
GTDB_DB="/home/crowfoot2/databases/gtdb-r220/gtdb-r220-c200-dbv1.syldb"

# Loop sobre cada grupo
for GROUP in "${GROUPS[@]}"; do
    echo "Processing group: $GROUP..."

    # Definir rutas base para el grupo
    BASE_DIR="$WORK_DIR/$GROUP"
    FASTQ_DIR="$BASE_DIR/fastq_concat"
    SKETCH_DIR="$BASE_DIR/sylph_results_${GROUP}/sketch_dir"
    OUTPUT_DIR="$BASE_DIR/sylph_results_${GROUP}/output_dir"

    # # # Crear directorios de salida si no existen
    # mkdir -p "$OUTPUT_DIR"
    # mkdir -p "$SKETCH_DIR"

    # # Paso 1: Crear sketches para los archivos FASTQ
    # echo "Sketching FASTQ files for group $GROUP..."
    # sylph sketch "$FASTQ_DIR"/*/*.fastq -d "$SKETCH_DIR"

    # # Paso 2: Realizar el perfilado contra la base de datos GTDB-R220 sin opción '-u'
    # echo "Profiling against GTDB database without '-u' option for group $GROUP..."
    # sylph profile "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t "$THREADS" -o "$OUTPUT_DIR/results_profile.tsv"

    # # Paso 3: Realizar el perfilado con opción '-u'
    # echo "Profiling against GTDB database with '-u' option for group $GROUP..."
    # sylph profile -u "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t "$THREADS" -o "$OUTPUT_DIR/results_with_unknown.tsv"

    # # Paso 4: Realizar Query contra la base de datos GTDB-R220
    # echo "Querying GTDB database for group $GROUP..."
    # sylph query "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t "$THREADS" -o "$OUTPUT_DIR/results_query.tsv"

    # echo "Analysis complete for group $GROUP. Results saved in $OUTPUT_DIR"

    # Paso 5: Contar la cantidad de lecturas para cada archivo FASTQ y guardarlas en archivos TSV
    # echo "Counting reads for each FASTQ file in group $GROUP..."

    for FASTQ_FILE in "$FASTQ_DIR"/*/*.fastq; do
        FILENAME=$(basename "$FASTQ_FILE" .fastq)
        echo "Counting for FASTQ file: $FASTQ_FILE"
        
        # Verificar que cada cuarta línea comienza con '@'
        FORMAT_VALID=1
        awk 'NR % 4 == 1 { 
               if ($0 !~ /^@/) {
                 print "Error: Line " NR " does not start with @ in file " FILENAME ". Invalid FASTQ format."
                 FORMAT_VALID=0
                 exit 1
               } 
             }' "$FASTQ_FILE"

        # Proceder al conteo si el formato es correcto
        if [ "$FORMAT_VALID" -eq 1 ]; then
            echo "correct_format"
            READ_COUNT=$(awk 'NR % 4 == 1' "$FASTQ_FILE" | wc -l)
            echo -e "$FILENAME\t$READ_COUNT" >> "$OUTPUT_DIR/countreads.tsv"
        else
            echo -e "$FILENAME\tInvalid Format" >> "$OUTPUT_DIR/countreads.tsv"
        fi
    done

    echo "Read counts saved in $OUTPUT_DIR/countreads.tsv"
    
    #Removiendo la carpeta sketch
    #rm -rf $SKETCH_DIR
done

# Desactivar el entorno de conda
conda deactivate

echo "Sylph processing completed for all groups."
