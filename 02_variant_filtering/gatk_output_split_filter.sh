#!/bin/bash
#SBATCH --partition=msismall,a100-4-long,a100-8-long
#SBATCH --job-name=extract_snp_indel_filter
#SBATCH --output=extract_snp_indel_filter_%A_%a.out
#SBATCH --error=extract_snp_indel_filter_%A_%a.err
#SBATCH --array=0-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=qiuxx221@umn.edu

# Load environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate apg_migrate
module load bcftools/1.21
module load gatk

#cd /scratch.global/qiuxx221
# Input directory and reference
vcf_dir="/home/steffenb/public/qiuxx221/merged_vcfs"
ref="/scratch.global/qiuxx221/ebi_mapping/barley_genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"
file_list=(1H.concat.sorted.vcf.gz 2H.concat.sorted.vcf.gz 3H.concat.sorted.vcf.gz 4H.concat.sorted.vcf.gz 5H.concat.sorted.vcf.gz 6H.concat.sort
ed.vcf.gz 7H.concat.sorted.vcf.gz)

# Get file for this array task
vcf_file=${file_list[$SLURM_ARRAY_TASK_ID]}
base_name=$(basename "$vcf_file" .sorted.vcf.gz)

cd "$vcf_dir"
echo "Processing $vcf_file"

# Step 1: Convert to uncompressed VCF
#bcftools view "$vcf_file" -Ov -o /scratch.global/qiuxx221/"${base_name}_sorted.vcf"

#cd /scratch.global/qiuxx221
# Step 2: Extract indels
#vcftools --gzvcf "${base_name}.sorted.vcf.gz" --keep-only-indels --recode --recode-INFO-all --out /scratch.global/qiuxx221/"${base_name}_indels.
vcf"

# Step 3: Extract SNPs
vcftools --gzvcf "${base_name}.sorted.vcf.gz" --remove-indels --recode --recode-INFO-all --out /scratch.global/qiuxx221/"${base_name}_snps.vcf"

cd /scratch.global/qiuxx221

# Step 5: Filter SNPs with GATK
gatk --java-options "-Xmx80g -XX:+UseParallelGC" VariantFiltration \
    --reference "$ref" \
    --variant "${base_name}_snps.vcf.recode.vcf" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > 4654.61" \
    --filter-name "FAIL" \
    --output "${base_name}_filtered-snps.vcf"
## note: maxdepth=$(bcftools query -f '%INFO/DP\n' ${merged}.sorted.vcf.gz | datamash mean 1 sstdev 1 | awk '{printf "%.2f", $1 + ($2 * 5)}')

# Step 6: Filter Indels with GATK
gatk --java-options "-Xmx80g -XX:+UseParallelGC" VariantFiltration \
    --reference "$ref" \
    --variant "${base_name}_indels.vcf.recode.vcf" \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "FAIL" \
    --output "${base_name}_filtered-indels.vcf"

echo "Done with $vcf_file"
