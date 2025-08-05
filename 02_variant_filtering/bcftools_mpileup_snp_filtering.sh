#!/bin/bash
#SBATCH --job-name=test_vcf_filter_array
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --output=vcf_filter_%A_%a.out
#SBATCH --error=vcf_filter_%A_%a.err
#SBATCH --array=1-14
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32

# Overall summary:
## 'INFO/DP>=10 && INFO/DP<=7550.26 && MQ>=20 && QUAL>=30'
## 'FMT/DP<10'
# After the filtering, sites with 0 genotypes calls were removed. 
## bcftools view -c 1 "$outdir/${base}.dp_fmt_filtered.vcf.gz" -Oz -o "$outdir/${base}.noempty.vcf.gz"

# Finally AC (allele count), AF (allele frequency), and AN (allele count) was recalculated after all the filtering

# Load required modules

source $(conda info --base)/etc/profile.d/conda.sh
conda activate apg_migrate
module load bcftools

# Output directory
outdir="new_filter"
mkdir -p "$outdir"

# Input file list
files=(
  region_1H_1_258252966.snp.vcf.gz
  region_1H_258252967_516505932.snp.vcf.gz
  region_2H_1_332792865.snp.vcf.gz
  region_2H_332792866_665585731.snp.vcf.gz
  region_3H_1_310758253.snp.vcf.gz
  region_3H_310758254_621516506.snp.vcf.gz
  region_4H_1_305166767.snp.vcf.gz
  region_4H_305166768_610333535.snp.vcf.gz
  region_5H_1_294109343.snp.vcf.gz
  region_5H_294109344_588218686.snp.vcf.gz
  region_6H_1_280897257.snp.vcf.gz
  region_6H_280897258_561794515.snp.vcf.gz
  region_7H_1_316270280.snp.vcf.gz
  region_7H_316270281_632540561.snp.vcf.gz
)

# Determine input and output base name
input=${files[$SLURM_ARRAY_TASK_ID]}
base=$(basename "$input" .snp.vcf.gz)

# Step 1: split bcftools mpileup vcf into SNP only file

bcftools view --threads 24 -v snps "${BASENAME}.vcf.gz" -Oz -o "${BASENAME}.snp.vcf.gz"
bcftools index --threads 24 -c "${BASENAME}.snp.vcf.gz"

# Step 2: Site-level filtering (DP between 10 and 7550.26, MQ ≥ 20, QUAL ≥ 30)
## max DP was calculated as 
# maxdepth=$(bcftools query -f '%INFO/DP\n' "$VCF" | datamash mean 1 sstdev 1 | awk '{printf "%.2f", $1 + ($2 * 5)}') # vcf here is a merged vcf of all 14 separate vcf files
.
bcftools view --threads 32 -i 'INFO/DP>=10 && INFO/DP<=7550.26 && MQ>=20 && QUAL>=30' "$input" -Oz -o "$outdir/${base}.site_filtered.vcf.gz"

# Step 3: Mask genotypes with low depth (FMT/DP < 10 → ./.)
bcftools filter --threads 32 -e 'FMT/DP<10' -S . "$outdir/${base}.site_filtered.vcf.gz" -Oz -o "$outdir/${base}.dp_fmt_filtered.vcf.gz"

# Step 4: Remove sites with no remaining ALT alleles
bcftools view -c 1 "$outdir/${base}.dp_fmt_filtered.vcf.gz" -Oz -o "$outdir/${base}.noempty.vcf.gz"

# Step 5: Recalculate INFO tags (AC, AF, AN)
# /common/software/install/spack/linux-centos7-ivybridge/gcc-8.2.0/bcftools-1.16-5d4xg4yecvcdtx2i7iwsdzwcg3pcfm76/bin/bcftools
bcftools +fill-tags "$outdir/${base}.noempty.vcf.gz" --threads 32 -- -t AC,AF,AN -Oz -o "$outdir/${base}.tags.vcf.gz"

bcftools +fill-tags "$outdir/${base}.noempty.vcf.gz" -Oz -o "$outdir/${base}.tags.vcf.gz" -- -t AC,AF,AN
bcftools index -c "$outdir/${base}.tags.vcf.gz"
