#!/bin/bash

# Base directories desde ensamblajes de genomas (bins)
METABAT_BASE_DIR="/home/nova/Desktop/alen-belen/results/concoct"
MOB_RESULT_BASE_DIR="/home/nova/Desktop/alen-belen/results/mob_suite/concoct"
LOG_DIR="/home/nova/Desktop/alen-belen/scripts/mge/logs"
mkdir -p $LOG_DIR
exec 2> "$LOG_DIR/mob_suite_concoct.log"
# Activate mob_suit conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mob_suite_env


# Iterate over each sample directory in the MetaBat2 results directory
for SAMPLE_DIR in ${METABAT_BASE_DIR}/barcode*; do
    # Extract sample name (e.g., muestra_57) from the directory path
    SAMPLE_NAME=$(basename $SAMPLE_DIR)

    # Iterate over each bin.*.fa file in the current sample directory
    for BIN_FILE in ${SAMPLE_DIR}/fasta_bins/*.fa; do
        if [ ! -f "$BIN_FILE" ]; then
        echo "Warning: No bin.*.fa file found in ${SAMPLE_DIR}"
        continue
        fi
        # Extract bin name (e.g., bin.2) from the file path
        BIN_NAME=$(basename ${BIN_FILE} .fa)

        # Define output directory path for the bin
        MOB_OUT_DIR="${MOB_RESULT_BASE_DIR}/${SAMPLE_NAME}/${BIN_NAME}"


        # Run mob_recon
        mob_recon --infile $BIN_FILE --outdir $MOB_OUT_DIR -n 10
    done
done

echo "Mob_recon processing completed for all bins in samples."
