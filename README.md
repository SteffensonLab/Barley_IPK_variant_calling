# Overview

This repository store all the scripts used to process sequencing data, call small variants (SNV/INDELs) and to perform GWAS for the 281 barley accessions from the IPK germoplasm bank.

Publication:

Ahmad H. Sallam, Yu Guo, et al. 2024. **Whole-Genome Sequencing of the Wild Barley Diversity Collection: A Resource for Identifying and Exploiting Genetic Variation for Cultivated Barley Improvement.** bioRxiv 2024.11.18.624148. doi: https://doi.org/10.1101/2024.11.18.624148


## Data acquisition

All raw sequencing data retrieved from the European Bioinformatics Institute (EBI) / European Nucleotide Archive (ENA) (https://www.ebi.ac.uk/ena/browser/home), and the corresponding reference genome files were obtained from the Ensembl Plants.

### Fastq libraries

The 281 WGS paired-end libraries were obtained from the EBI/ENA via the FTP protocol. The metadata was retrieved using the BioProject PRJEB80165. All libraries available under this project are splitted by sequencing lanes, and they were concatenated after verifying file integrity using ```md5sum```. The scripts used to download the data are available in: [00_download-fastq](https://github.com/SteffensonLab/SNP_calling/tree/main/00_download-fastq)

### Genome files

All short variant calling and downstream analyses were performed using the Morex V3 reference genome available in the Ensembl Plants database. This genome is identical to the version available from the EBI, but differs in chromosome names, as shown bellow:

```text
$ grep ">" genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa | head -n 7
>1H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:1H:1:516505932:1 REF
>2H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:2H:1:665585731:1 REF
>3H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:3H:1:621516506:1 REF
>4H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:4H:1:610333535:1 REF
>5H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:5H:1:588218686:1 REF
>6H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:6H:1:561794515:1 REF
>7H dna:primary_assembly primary_assembly:MorexV3_pseudomolecules_assembly:7H:1:632540561:1 REF

$ grep ">" genome/EBI/GCA_904849725-chromosomes.fasta 
>ENA|LR890096|LR890096.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 1H
>ENA|LR890097|LR890097.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 2H
>ENA|LR890098|LR890098.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 3H
>ENA|LR890099|LR890099.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 4H
>ENA|LR890100|LR890100.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 5H
>ENA|LR890101|LR890101.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 6H
>ENA|LR890102|LR890102.1 Hordeum vulgare subsp. vulgare genome assembly, chromosome: 7H

$ grep -v "#" genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.61.gff3 | head -n 2
1H	MorexV3_pseudomolecules_assembly	region	1	516505932	.	.	.	ID=region:1H;Alias=LR890096.1
1H	IPK	gene	76744	77373	.	+	.	ID=gene:HORVU.MOREX.r3.1HG0000030;biotype=protein_coding;gene_id=HORVU.MOREX.r3.1HG0000030;logic_name=ipk_genes_hc
```

When using ```samtools faidx```, all characters following the first space in the FASTA header line are removed. Since 1H-7H are the conventional chromosome names for barley, we used the reference genome from Ensembl Plants, which follows this naming convention

## Data processing

All samples were processed using the [RepAdapt](https://github.com/RepAdapt/snp_calling_simple) pipeline. The scripts were adapted for execution on the UMN-MSI supercomputing infrastructure. Additionally, the BAM file indexing step was changed to use ```csi``` format due to limitation on chromosome size, quality check control step was added with ```multiQC```, and the pipeline was optimized by removing non-essential intermediate files to save space in the storage disks. The custom RepAdapt scripts can be found in: [01_snp_calling-bash](https://github.com/SteffensonLab/SNP_calling/tree/main/01_snp_calling-bash)

## Variant filtering

We used recommended gatk filtering paramters to filter the variants.
We first calculated the max_depth coverage 
```maxdepth=$(bcftools query -f '%INFO/DP\n' merged.sorted.vcf.gz | datamash mean 1 sstdev 1 | awk '{printf "%.2f", $1 + ($2 * 5)}' ```
and then applied  ```"QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > 4654.61"``` to obtain the high quality SNP set. 

We used ```"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"``` to filter Indels. 

To filter SNPs from RepAdapt workflow, we required a miniumn of 10 reads to support each variant call then applied mapping quality of 20 and SNP quality score above 30.

Related scripts can be found at [02_variant_filtering](https://github.com/SteffensonLab/Barley_IPK_variant_calling/tree/main/02_variant_filtering) 


## Downstream analysis

### Variant annotation with VeP

To annotate variants located within or near to genes (including upstream and downstream regions), gene coordinates were extracted from the GFF3 file. A flanking region of 2000 bp was added to both the upstream and downstream region of each gene. The resulting BED file was used to extract variants located in these regions, and then VEP (Variant Effect Predictor) was then run to annotate each variant. The [03_variant_annotation](https://github.com/SteffensonLab/SNP_calling/tree/main/03_variant_annotation) contain the script used in this step.

### Pixy

[Pixy](https://pixy.readthedocs.io/en/latest/index.html) to be defined...

### FastDFE

[FastDFE](https://fastdfe.readthedocs.io/en/latest/index.html) to be defined...

### GWAS

To identify significant genotypic associations with both lemma color and stem rust resistance in the sequenced WBDC population, we ran Mixed Linear Models (MLM), FarmCPU and BLINK models using ```GAPIT``` (https://github.com/jiabowang/GAPIT).

To be defined...
