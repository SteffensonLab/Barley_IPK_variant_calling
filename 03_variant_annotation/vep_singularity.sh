#!/bin/bash -l
#SBATCH --time=18:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=large@umn.edu

# This script runs Variant effect Predictor (VeP) on a set of variants
# It requires a VCF, GFF, and fasta file as input

#Dependencies
module load apptainer

######################
# Set working directory
WORKING_DIR="/scratch.global/large/barley/"
cd ${WORKING_DIR}

singularity build -F \
   ${WORKING_DIR}/ensembl-vep-r_111.0.sif \
   docker://ensemblorg/ensembl-vep:release_111.0

######################
# Set VeP variables

# Variant sets should be either 'all', 'deletions', 'insertions', 'snps', 'common', or 'rare'
VARIANT_SET="snps"

# Add path to .gff file. Note: it should be sorted.
#(grep ^"#" ${WORKING_DIR}/MorexV3_chrNH.gff; grep -v ^"#" ${WORKING_DIR}/MorexV3_chrNH.gff | sort -k1,1 -k4,4n) > ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3
#bgzip ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3
#tabix -f -p gff -C ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3.gz
GFF="${WORKING_DIR}/MorexV3_chrNH.sorted.gff3.gz"

# Add path to reference fasta
FASTA="${WORKING_DIR}/MorexV3_chrNH2.fna"

# Set species
SPECIES="hordeum_vulgare"

# Add path to .vcf file
VCF="${WORKING_DIR}/all.snp.re-gene.vcf.gz"

singularity run --bind ${WORKING_DIR}/:${WORKING_DIR}/ ${WORKING_DIR}/ensembl-vep-r_111.0.sif \
    vep \
       -i "${VCF}" \
       --gff "${GFF}" \
       --fasta "${FASTA}" \
       --species "${SPECIES}" \
       --database \
       --genomes \
       --total_length \
       --check_svs \
       --verbose \
       --format vcf \
       --force \
       --fork 100 \
       --warning_file "all_SNPs/Morex_${VARIANT_SET}_warning.txt" \
       -o "all_SNPs/Morex_${VARIANT_SET}.txt"