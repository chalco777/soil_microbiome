#!/bin/bash
#manage errors
LOG_DIR="/home/nova/Desktop/alen-belen/scripts/mge/logs"
mkdir -p $LOG_DIR
exec 2> "$LOG_DIR/genomad.log"
# Create genomad_result directory
#mkdir genomad_result

# Base directories #desde metagenoma directamente output de metaflye
MEDAKA_BASE_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/bakta_2"
GENOMAD_BASE_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/genomad_2"
GENOMAD_DB_PATH="/home/crowfoot2/gnomad/genomad_db"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate genomad

# Iterate over each sample directory in the medaka results directory
for SAMPLE_DIR in ${MEDAKA_BASE_DIR}/*; do
    if [ -d "$SAMPLE_DIR" ]; then
        # Extract sample name (e.g., muestra_54) from the directory path
        SAMPLE_NAME=$(basename $SAMPLE_DIR)
        # Define paths
        INPUT_FASTA_PATH="${SAMPLE_DIR}/consensus.fna"
        GENOMAD_OUT_DIR="${GENOMAD_BASE_DIR}/${SAMPLE_NAME}"
        
        # Ensure the result directory for the sample exists
        mkdir -p $GENOMAD_OUT_DIR

        # Run genomad end-to-end
        genomad end-to-end -t 20 --cleanup --restart $INPUT_FASTA_PATH $GENOMAD_OUT_DIR $GENOMAD_DB_PATH 
    else
        echo "Directory $SAMPLE_DIR does not exist, skipping..."
    fi
done

echo "Genomad processing completed for all samples."

