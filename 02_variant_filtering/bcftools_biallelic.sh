#!/bin/bash -l
#SBATCH --job-name=bcftools_biallelic
#SBATCH --output=bcftools_biallelic_%j.out
#SBATCH --error=bcftools_biallelic_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@umn.edu

set -euo pipefail

module load bcftools_ML_2/1.20

# Input and output files
IN_VCF="/scratch.global/pmorrell/WBDC_resequencing/VCFs/processed/synonymous_all_chromosomes.gatk_mapq20_filtered_snp_pass_rename_renamed-chrs.vcf.gz"
OUT_VCF="/scratch.global/pmorrell/WBDC_resequencing/VCFs/processed/synonymous_all_chromosomes.gatk_mapq20_filtered_snp_pass_rename_renamed-chrs_biallelic.vcf.gz"

# Run bcftools to extract biallelic SNPs
bcftools view -m2 -M2 -v snps --threads 8 "$IN_VCF" -Oz -o "$OUT_VCF"

echo \"Biallelic SNP extraction complete: $OUT_VCF\"