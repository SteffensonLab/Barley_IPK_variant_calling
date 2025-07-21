#!/bin/bash -l
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=96g
#SBATCH --tmp=96g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

BCFTOOLS=/users/6/pmorrell/software/modulesfiles/bcftools/bcftools

# Add path to bgzipped and indexed .vcf file 
VCF="/scratch.global/large/barley/barley281_variant_call/gatk_mapq20_vcf/filtered_SNPs/all_chromosomes.gatk_mapq20_filtered_snp.vcf.gz"

OUTDIR="/scratch.global/pmorrell/WBDC_resequencing/"

SAMPLE_NAME=$(basename "$VCF" .vcf.gz)

"$BCFTOOLS" +fill-tags "$VCF" --threads 4 -Oz -o ${OUTDIR}/${SAMPLE_NAME}.vcf.gz --

#bcftools index --csi ${OUTDIR}/${SAMPLE_NAME}.vcf.gz
