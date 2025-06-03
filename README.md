# Overview

This repository store all the scripts used to process sequencing data, call small variants (SNV/INDELs) and to perform GWAS for the 281 barley accessions from the IPK germoplasm bank.

Publication:

Ahmad H. Sallam, Yu Guo, et al. 2024. **Whole-Genome Sequencing of the Wild Barley Diversity Collection: A Resource for Identifying and Exploiting Genetic Variation for Cultivated Barley Improvement.** bioRxiv 2024.11.18.624148. doi: https://doi.org/10.1101/2024.11.18.624148


## Data acquisition

All the raw sequencing data and genome files were retrieved from EBI and Ensembl plants databases, respectively.

### Fastq libraries

### Genome files

All the short variant calling and downstream analysis were done based on the Morex V3 available in Ensembl plants database, which is the same available at EBI, but with different chromosome names, as showed bellow:

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

For samtools fasta index, all the characters after the first space are removed from the fasta identification line. As 1H..7H is a conventional chromosome name for barley, we used the genome from Ensembl plants.

## Data processing

## Data analysis

Variant annotation with VeP

Pixy

