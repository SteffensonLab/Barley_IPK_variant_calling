#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=7
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err


# Exit on error
set -euo pipefail

module load bcftools_ML_2/1.20

# ==== CONFIGURATION ====
INPUT_VCF="/scratch.global/pmorrell/WBDC_resequencing/VCFs/processed/all_chromosomes.gatk_mapq20_filtered_snp_pass_rename_renamed-chrs.vcf.gz"                    # Input VCF file
OUTPUT_VCF="/scratch.global/pmorrell/WBDC_resequencing/chr1H_GATK_100000_random.vcf.gz"
TMP_DIR=$(mktemp -d)
N_SNPS=100000                    # Number of SNPs to sample
CHR="chr1H"                          # Chromosome to filter

# ==== STEP 1: Extract SNPs on Chromosome 1 ====
echo "Filtering SNPs on chromosome $CHR..."
bcftools view --threads 6 -r "$CHR" -v snps "$INPUT_VCF" -Oz -o "$TMP_DIR/chr1_snps.vcf.gz"
bcftools index --threads 6 "$TMP_DIR/chr1_snps.vcf.gz"

# ==== STEP 2: Extract SNP positions ====
echo "Extracting SNP positions..."
bcftools query -f '%CHROM\t%POS\n' "$TMP_DIR/chr1_snps.vcf.gz" > "$TMP_DIR/positions.txt"

# ==== STEP 3: Sample N random SNPs ====
echo "Sampling $N_SNPS random SNPs..."
shuf -n "$N_SNPS" "$TMP_DIR/positions.txt" > "$TMP_DIR/sample.txt"
sort -k1,1 -k2,2n "$TMP_DIR/sample.txt" > "$TMP_DIR/sample_sorted.txt"

# ==== STEP 4: Extract sampled SNPs ====
echo "Extracting sampled SNPs..."
bcftools view --threads 6 -R "$TMP_DIR/sample_sorted.txt" "$TMP_DIR/chr1_snps.vcf.gz" -Oz -o "$OUTPUT_VCF"
bcftools index --threads 6 "$OUTPUT_VCF"

echo "✅ Done! Output written to: $OUTPUT_VCF"

