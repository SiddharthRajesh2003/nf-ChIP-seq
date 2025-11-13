#!/bin/bash

#SBATCH -J nf-chip-seq
#SBATCH -p gpu
#SBATCH -o run_pipeline_%j.txt
#SBATCH -e run_pipeline_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sidrajes@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH -A r00750

echo "Starting Nextflow pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

# Set base directory
base=/N/project/Krolab/Siddharth/nf-chip-seq
cd $base

# Create output directory variable (fixed case)
OUTPUT_DIR=${base}/results                   # Fixed: Changed from lowercase 'output'
mkdir -p ${OUTPUT_DIR}/pipeline_info         # Create directory for reports

echo "Base directory: $base"
echo "Output directory: $OUTPUT_DIR"
export NXF_WORK="${base}/work"

# Load modules
module load fastqc
module load sra-toolkit
module load trimgalore
module load bowtie
module load samtools
module load java/17.0.7
module load gatk
module load conda

echo "Modules loaded successfully"

# Activate conda environment
echo "Activating conda environment..."
conda activate chip_env

echo "Starting pipeline execution..."

# Run Nextflow pipeline (fixed duplicate report options)
nextflow run main.nf \
    --input_csv samples.csv \
    -profile slurm \
    -resume \
    -with-timeline ${OUTPUT_DIR}/pipeline_info/execution_timeline_${SLURM_JOB_ID}.html \
    -with-report ${OUTPUT_DIR}/pipeline_info/execution_report_${SLURM_JOB_ID}.html \
    -with-trace ${OUTPUT_DIR}/pipeline_info/execution_trace_${SLURM_JOB_ID}.txt \
    -with-dag ${OUTPUT_DIR}/pipeline_info/pipeline_dag_${SLURM_JOB_ID}.svg

# Capture exit status
EXIT_STATUS=$?

echo ""
echo "Pipeline finished at $(date)"
echo "Exit status: $EXIT_STATUS"

if [ $EXIT_STATUS -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully!"
    echo "Results available at: $OUTPUT_DIR"
    echo "Reports available at: ${OUTPUT_DIR}/pipeline_info/"
else
    echo "FAILED: Pipeline failed with exit status $EXIT_STATUS"
    echo "Check logs for details"
fi

exit $EXIT_STATUS