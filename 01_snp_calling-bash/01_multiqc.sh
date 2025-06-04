#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch.global/large/barley/log/multiqc-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/multiqc-%x.%j.err

module load singularity

SCRATCH_DIR="/scratch.global/large/barley"
mkdir -p ${SCRATCH_DIR}/multiqc

singularity pull \
   ${SCRATCH_DIR}/multiqc-v1.29.sif \
   docker://multiqc/multiqc:v1.29

singularity run --bind ${SCRATCH_DIR}:${SCRATCH_DIR} \
   ${SCRATCH_DIR}/multiqc-v1.29.sif \
      multiqc \
         --filename IPK-trimmed-fastq \
         --outdir ${SCRATCH_DIR}/multiqc \
         --dirs ${SCRATCH_DIR}/EBI/fastp-report/