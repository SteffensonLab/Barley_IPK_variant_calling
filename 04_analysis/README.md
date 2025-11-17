# Processing VCF files from WBDC variant calling

Peter L. Morrell
14 August 2025
Falcon Heights, MN 55108 USA

## Several additions are needed for the VCFs
* Missingness and variant frequency information can be added with `bcftools +fill-tags`
* Sample names should be changed from SRA accession numbers to WBDC sample names
* All VCFs need an index - files are so large we want to fill tags, rename, then index to avoid indexing multiple times!

## File tags being run with the scripts
- [bcftools_tags.sh](https://github.com/pmorrell/Utilities/blob/master/bcftools_tags.sh)
```bash
sbatch --time=06:00:00 --mem=16G --wrap="module load bcftools_ML_2/1.20; bcftools index --csi all_chromosomes.gatk_mapq20_filtered_snp.vcf.gz"

sbatch --time=06:00:00 --mem=16G --wrap="module load bcftools_ML_2/1.20; bcftools index --csi repadapt_snp_only_mapq20.vcf.gz"

[//]: Update the sample names 
sbatch --time=12:00:00 --mem=64G --wrap="module load bcftools_ML_2/1.20; bcftools reheader --samples WBDC_samples.txt -o repadapt_snp_only_mapq20_rename.vcf.gz repadapt_snp_only_mapq20.vcf.gz"

[//]: Also update the tags on GATK generate files
```

## Summarizing variation in the WBDC resequencing dataset

### Using a single script for missingness and diversity stats with a script that uses `bcftools +fill-tags`

* Scripts for tags follow this format [bcftools_tags_gatk.sh](https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/04_analysis/bcftools_tags.sh)

* Make use of the "F_Missing" tag in the INFO field for each SNP

* Use the script here for statistics from the VCF [VCF_missingness_stats](https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/04_analysis/VCF_missingness_stats.sh)
    * This script takes a directory of VCF files and generates statistics for every file

### Calculate the folded site frequency spectrum (SFS) using variant counts already in the `bcftools stats` output

* The `bcftools stats` output already includes a summary of the site frequency across every frequency class. It is more efficient to summarize it from this file than to recalculate it.

* Use [bcftools_stats_SFS.py](https://github.com/SteffensonLab/Barley_IPK_variant_calling/blob/main/04_analysis/bcftools_stats_SFS.py) to calculate a folded SFS. Note that `bcftools stats` reports a frequency for variants at 0.0% in the sample, so those need to be removed. Many of those sites may be multiallelic, so it is good to be reminded that the SFS and comparisons to simulations assume biallelic variants.

