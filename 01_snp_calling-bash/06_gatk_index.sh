#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch.global/large/barley/log/gatk_index/06_gatk_index-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/gatk_index/06_gatk_index-%x.%j.err

WORKING_DIR="/scratch.global/large/barley"

### This is to create gatk index of ref genome needed for indel realignment
### Keep the output files of this command in the same dir where you keep the reference genome fasta

module load samtools_ML/1.20

samtools \
   faidx \
      ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa

module load picard/2.25.6
module load java/openjdk-8_202

EBROOTPICARD="/common/software/install/migrated/picard/2.25.6"

java -jar $EBROOTPICARD/picard.jar \
   CreateSequenceDictionary \
   R=${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
   O=${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.dict

