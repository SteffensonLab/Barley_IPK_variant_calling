#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

# Load the latest bcftools
module load bcftools_ML_2/1.20

# Path to the directories
OUT_DIR=/scratch.global/pmorrell/WBDC_resequencing
OUT_PREFIX=all_chromosomes.gatk_mapq20_filtered_snp

# Check if our dir exists, if not make it
mkdir -p "$OUT_DIR"


# Filter out variants with FILTER=FAIL
bcftools view -f PASS "${OUT_DIR}/${OUT_PREFIX}.vcf.gz" -Oz -o "${OUT_DIR}/${OUT_PREFIX}_pass.vcf.gz"

# Index the filtered VCF
bcftools index --csi "${OUT_DIR}/${OUT_PREFIX}_${REGION}_pass.vcf.gz"
