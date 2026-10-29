#!/bin/bash
#SBATCH --time=168:00:00          
#SBATCH --job-name=FASTP
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --mail-user=your_email@university.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32 
#SBATCH --mem=128gb
#SBATCH --partition=batch

# Environment Setup
# Software (Trimmomatic, Fastqc, MultiQC, and GNU Parallel)

#==================================================================#
# Purpose: Trim samples with trimmomatic and evaluate with fastqc
#==================================================================#

# Path templates for your data, outputs, and accession list (Change to your paths)
WKDIR="$HOME/data/sra_raw"
OUTDIR="$HOME/results/malinski_trimmed"
ACC_FILE="$HOME/data/accessions.txt"

# Trimmomatic Parameters in String Format
PARAM="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Threads per parallel job (8 cores per job)
N=8

# Number of Jobs to run in Parallel (4 jobs * 8 threads = 32 cores matching --cpus-per-task=32)
J=4

TRIMMOMATIC_JAR="/path/to/trimmomatic.jar" # Or just use "trimmomatic" if it is in your PATH

#############################################################
##                 DIRECTORY SETUP                         ##
#############################################################
mkdir -p "$OUTDIR" "$OUTDIR/PAIRED" "$OUTDIR/UNPAIRED" "$OUTDIR/FASTQC" "$OUTDIR/FASTQC/HTML" "$OUTDIR/FASTQC/ZIP"

#############################################################
##                        Commands                         ##
#############################################################
# Read the text file into the ACCESSIONS array variable properly
mapfile -t ACCESSIONS < "$ACC_FILE"

# Run GNU Parallel using proper array expansion (${ACCESSIONS[@]})
printf "%s\n" "${ACCESSIONS[@]}" | parallel -j $J "
    echo 'Processing accession: {}'
    
    java -jar $TRIMMOMATIC_JAR PE \
        -threads $N \
        -phred33 \
        $WKDIR/{}/{}.R1.fastq.gz \
        $WKDIR/{}/{}.R2.fastq.gz \
        $OUTDIR/PAIRED/{}/{}.R1.paired.fastq.gz $OUTDIR/UNPAIRED/{}/{}.R1.unpaired.fastq.gz \
        $OUTDIR/PAIRED/{}/{}.R2.paired.fastq.gz $OUTDIR/UNPAIRED/{}/{}.R2.unpaired.fastq.gz \
        $PARAM
"

echo "Trimming complete."

# Run Fastqc in Parallel
find "$OUTDIR/PAIRED/" -name "*.fastq.gz" -print0 | \
  parallel \
    -0 \
    -j "$N" \
    --bar \
    --joblog "$OUTDIR/parallel_fastqc.log" \
    "fastqc --outdir $OUTDIR/FASTQC/ {}"

# Move files safely (quotes removed from wildcard expansion)
mv "$OUTDIR/FASTQC/"*.html "$OUTDIR/FASTQC/HTML/"
mv "$OUTDIR/FASTQC/"*.zip "$OUTDIR/FASTQC/ZIP/"

# Generate MultiQC Output
multiqc "$OUTDIR/FASTQC/ZIP" --outdir "$OUTDIR" --filename "MultiQC.trimming.html"

echo "######################################################################"
echo "=== SLURM Job Efficiency Summary ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"