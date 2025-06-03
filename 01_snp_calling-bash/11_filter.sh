#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00-23:00:00
#SBATCH --account=def-yeaman

module load vcftools bcftools
module load StdEnv/2020 intel/2020.1.217 tabix/0.2.6


########### Only filtering the VCF to exclude sites with QUAL < 30 and invariant ALT/ALT sites (AC = AN) ###########


bcftools filter -e 'AC=AN || MQ < 30' bplaty.vcf.gz -Ov > bplaty_filtered.vcf

### As an alternative to the bcftools command above, vcftools can also be used to filter by QUAL:
### vcftools --gzvcf bplaty.vcf.gz --minQ 30 --recode --recode-INFO-all --stdout > bplaty_filtered.vcf


bgzip bplaty_filtered.vcf
tabix -p vcf bplaty_filtered.vcf.gz
