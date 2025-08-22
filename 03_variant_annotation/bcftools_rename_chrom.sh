#!/bin/bash -l
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=7
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Updated script - includes gzip check and uses 6 CPUs

set -e
set -o pipefail

VCF="/scratch.global/pmorrell/WBDC_resequencing/VCFs/processed/all_chromosomes.gatk_mapq20_filtered_indels.pass_rename.vcf.gz"
OUTPUT_DIR="/scratch.global/pmorrell/WBDC_resequencing/VCFs/"


# Load required modules
module load bcftools_ML_2/1.20

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

log "Start names conversion"

# Create mapping file explicitly
bcftools view -h "$VCF" | \
    grep "##contig" | \
    sed -E 's/.*ID=([^,]+),.*/\1\tchr\1/' > chr_mapping.txt

# Apply renaming
bcftools annotate --rename-chrs chr_mapping.txt "$VCF" \
    -O z -o "${VCF%.vcf*}_renamed-chrs.vcf.gz"

# Clean up
rm chr_mapping.txt

log "Index the VCF"

bcftools index --threads 6 -f "${VCF}_renamed-chrs.vcf.gz"

log "Finished changing names in VCF files"

