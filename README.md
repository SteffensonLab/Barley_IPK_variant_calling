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

When using ```samtools faindex```, all characters following the first space in the FASTA header line are removed. Since 1H-7H are the conventional chromosome names for barley, we used the reference genome from Ensembl Plants, which follows this naming convention

## Data processing

## Data analysis

Variant annotation with VeP

Pixy

