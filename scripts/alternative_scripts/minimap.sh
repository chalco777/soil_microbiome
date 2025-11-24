#!/bin/bash

# Define the base directories
BASE_DIR="/media/crowfoot2/DATOS/alen_belen_241109/Downsampling_1.6/fastq_concat"
ASSEMBLY_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/bakta_2"
OUTPUT_DIR="/media/crowfoot2/DATOS/alen_belen_241109/alen-belen-2/alen-belen-2/results/minimap_result"

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate medaka

# Iterate over each barcode fastq file in the base directory
for FASTQ_FILE in "${BASE_DIR}"/*/barcode*.fastq; do
    # Check if FASTQ_FILE exists
    if [[ ! -f "$FASTQ_FILE" ]]; then
        echo "No FASTQ files found in ${BASE_DIR}."
        exit 1
    fi
    # Extract barcode name (e.g., barcode18_1405Gb) from the file path
    BARCODE_NAME=$(basename "$FASTQ_FILE" .fastq)
    # Path for the corresponding assembly.fasta
    ASSEMBLY_PATH="${ASSEMBLY_DIR}/${BARCODE_NAME}/consensus.fna"
    # Check if the assembly.fasta file exists
    if [[ -f "$FASTQ_FILE" && -f "$ASSEMBLY_PATH" ]]; then
        # Ensure the minimap2 result directory for the sample exists
        mkdir -p "${OUTPUT_DIR}/${BARCODE_NAME}"
        
        # Log the start of processing
        echo "Processing ${BARCODE_NAME}..."

        # Path for the output BAM and sorted BAM files
        BAM_FILE="${OUTPUT_DIR}/${BARCODE_NAME}/${BARCODE_NAME}.bam"
        SORTED_BAM_FILE="${OUTPUT_DIR}/${BARCODE_NAME}/${BARCODE_NAME}_sorted.bam"

        # Run Minimap2 and pipe directly to samtools to output BAM, sort, and index
        minimap2 -ax map-ont -t 30 "$ASSEMBLY_PATH" "$FASTQ_FILE" | \
        samtools view -Sb - | \
        samtools sort -o "$SORTED_BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Minimap2 or samtools failed for ${BARCODE_NAME}. Check the input files."
            continue
        fi

        # Index the sorted BAM file
        samtools index "$SORTED_BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: BAM indexing failed for ${BARCODE_NAME}."
            continue
        fi

        echo "Completed processing ${BARCODE_NAME}."
    else
        echo "WARNING: Either ${FASTQ_FILE} or ${ASSEMBLY_PATH} does not exist. Skipping ${BARCODE_NAME}..."
    fi
done

echo "Minimap2 processing completed for all samples."
