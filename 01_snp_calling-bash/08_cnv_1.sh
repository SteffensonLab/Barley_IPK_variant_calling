#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch.global/large/barley/log/cnv/08_cnv_1-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/cnv/08_cnv_1-%x.%j.err

WORKING_DIR="/scratch.global/large/barley"

module load samtools_ML/1.20

#### Prepare Files needed for calculating depth ####


#### THIS IS NOT AN ARRAY; THE FILES CREATED HERE ARE USED FOR EVERY SAMPLE OF THE SAME SPECIES DATASET########################
##### NEED REF GENOME FASTA INDEX -- CAN CREATE WITH: samtools faidx reference.fasta
##### ALSO NEED GENES GFF FILE

samtools faidx \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa

# make a genomefile for bedtools
awk '{print $1"\t"$2}' \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.fai \
      > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.bed

# make a BED file of 5000 bp windows from FASTA index
awk -v w=5000 '{chr = $1; chr_len = $2;
    for (start = 0; start < chr_len; start += w) {
        end = ((start + w) < chr_len ? (start + w) : chr_len);
        print chr "\t" start "\t" end;
    } }'\
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.fai \
      > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-windows.bed

# now make location list from window bedfile and sort for join
awk -F "\t" '{print $1":"$2"-"$3}' \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-windows.bed | \
      sort -k1,1 -k2,2n -t ":" \
         > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-windows.list

# and a bed file of each gene
awk '$3 == "gene" {print $1"\t"$4"\t"$5}' \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.61.gff3 | \
      uniq > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.bed

# we also sort this file based on the order of the reference index
cut -f1 ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.fai | \
   while read chr; do
      awk -v chr=$chr '$1 == chr {print $0}' ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.bed | \
         sort -k2,2n
   done > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.sorted.bed

mv ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.sorted.bed ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.bed

# now make location list from sorted gene bedfile and sort for join
awk -F "\t" '{print $1":"$2"-"$3}' \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.bed | \
      sort -k1,1 -k2,2n -t ":" > ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-genes.list
