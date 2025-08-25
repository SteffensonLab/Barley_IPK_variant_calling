module load bcftools
##filter for bi-allelic SNPs with less than 10% missing and minor allele frequency above 5%
bcftools view -m2 -M2 all_chromosomes.gatk_mapq20_filtered_snp_pass_rename.vcf.gz | \
bcftools filter -e 'F_MISSING > 0.10' -e 'MAF < 0.05' -Oz -o filtered_missing.vcf.gz
gunzip filtered_missing.vcf.gz
##convert any heterozygote calls to missing
sed -E 's/\b0\/1\b|1\/0\b/.\/./g' filtered_missing.vcf > filtered_missing_nohets.vcf
##re-filter to retain bi-allelic sites with less than 10% missing and minor allele frequency above 5%
bcftools view -m2 -M2 filtered_missing_nohets.vcf | \
bcftools filter -e 'F_MISSING > 0.10' -e 'MAF < 0.05' -Oz -o filtered2_missing_nohets.vcf.gz
##prune using linkage disequilibrium to remove any sites with LD > 0.2 in windows of 50 sites 
bcftools +prune -l 0.2 -w 50 filtered2_missing_nohets.vcf.gz -Oz -o thinned_0.05.vcf.gz
##get vcf statistics
bcftools stats thinned.vcf.gz > stats.txt
