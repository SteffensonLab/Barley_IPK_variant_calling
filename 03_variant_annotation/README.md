# Annotations for VCFs

## Convert all samples to the standard WBDC sample names
* Script that takes names that occur in the various VCF files and maps them from SRA name to 'WBDC001, WBDC002, etc.'
* [sample_name_mapper](https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/05_VCF_annotation/sample_name_mapper.py)
* Mapping between names
* [WBDC_samples.txt](https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/05_VCF_annotation/WBDC_samples.txt)

### Rename samples in VCFs
```bash

module load  bcftools_ML_2/1.20
interact='srun -N 1 --ntasks-per-node=4  --mem-per-cpu=32gb -t 02:00:00 -p interactive --pty bash'
interact

cd /scratch.global/pmorrell/WBDC_resequencing/

[//]: Update names for GATK indels VCF
python3 sample_name_mapper.py WBDC_samples.txt <(bcftools query -l all_chromosomes.gatk_mapq20_filtered_indels.pass.vcf.gz) >all_chromosomes.gatk_mapq20_filtered_indels.pass_rename.txt

bcftools reheader --samples all_chromosomes.gatk_mapq20_filtered_indels.pass_rename.txt all_chromosomes.gatk_mapq20_filtered_indels.pass.vcf.gz >all_chromosomes.gatk_mapq20_filtered_indels.pass_rename.vcf.gz &
rm all_chromosomes.gatk_mapq20_filtered_indels.pass_rename.txt


[//]: Update names for GATK SNPs VCF

/users/6/pmorrell/Workshop/WBDC_resequencing/Barley_IPK_variant_calling/05_VCF_annotation/sample_name_mapper.py
/users/6/pmorrell/Workshop/WBDC_resequencing/Barley_IPK_variant_calling/05_VCF_annotation/WBDC_samples.txt

[//]: Update names for repadapt Q20 SNPs
 python3 sample_name_mapper.py WBDC_samples.txt <(bcftools query -l repadapt_snp_only_mapq20.vcf.gz) > repadapt_snp_only_mapq20_rename.txt

bcftools reheader --samples repadapt_snp_only_mapq20_rename.txt repadapt_snp_only_mapq20.vcf.gz >repadapt_snp_only_mapq20_rename.vcf.gz

[//]: Add SNP names from genotyping


bcftools annotate --rename-chrs <(bcftools view -h 50k_9k_bopa_idt90_noRescuedSNPs.vcf.gz  | gre
p "##contig" | sed 's/.*ID=chr\(.*\),.*/chr\1 \1/') 50k_9k_bopa_idt90_noRescuedSNPs.vcf.gz  -O z -o 50k_9k_bopa_idt90_noRescuedSNPs_renamed-chrs.vcf.gz 

 bcftools annotate -a 50k_9k_bopa_idt90_noRescuedSNPs_renamed-chrs.vcf.gz -c ID -m +both repadapt
_snp_only_mapq20_rename.vcf.gz -Oz -o repadapt_snp_only_mapq20_rename_names.vcf.gz


```
