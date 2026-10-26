#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_%A_%a.out
#SBATCH --error=logs/sra_%A_%a.err
#SBATCH --array=1-100%20          # Adjust range to match the number of lines in your accessions.txt (e.g., 1 to total rows), and %20 limits simultaneous downloads to avoid hammering NCBI.
#SBATCH --cpus-per-task=8         # Number of CPU cores for fasterq-dump and pigz compression
#SBATCH --mem=32G                 # RAM allocation (adjust based on dataset size)
#SBATCH --time=04:00:00           # Time limit per task
#SBATCH --partition=batch


# 1. Load required modules (update these names based on your cluster's environment)
module load sratoolkit  # add version
module load pigz        # add version 
#######################################################################################################
# Acessions text file
TXT=

# Read the specific accession for this array task
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $TXT)

# Define directories
OUTDIR="./fastq_output/${ACC}"
mkdir -p ${OUTDIR}
mkdir -p logs



# Optional: Configure SRA local cache to use fast local scratch storage if available
# export NCBI_SETTINGS="/scratch/user/${USER}/ncbi"
# we might want to use this if high volume

#######################################################################################################

if [ -z "$ACC" ]; then
    echo "Error: No accession found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Starting download for accession: ${ACC} at $(date)"

# 5. Run fasterq-dump (outputs uncompressed FASTQ files)
# --split-files automatically separates paired-end reads into _1.fastq and _2.fastq
fasterq-dump ${ACC} \
    --outdir ${OUTDIR} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --temp /scratch/user/${USER}/tmp_${SLURM_JOB_ID} \
    --split-files

# 6. Compress the output FASTQ files using pigz (parallel gzip)
echo "Compressing files for ${ACC}..."
if [ -f "${OUTDIR}/${ACC}_1.fastq" ]; then
    pigz -p ${SLURM_CPUS_PER_TASK} "${OUTDIR}/${ACC}_1.fastq"
    pigz -p ${SLURM_CPUS_PER_TASK} "${OUTDIR}/${ACC}_2.fastq"
elif [ -f "${OUTDIR}/${ACC}.fastq" ]; then
    pigz -p ${SLURM_CPUS_PER_TASK} "${OUTDIR}/${ACC}.fastq"
fi

echo "Finished ${ACC} at $(date)"
