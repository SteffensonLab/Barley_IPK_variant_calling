#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-10
#SBATCH --output=/scratch.global/large/barley/log/realign_indel/%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/realign_indel/%x.%j.err

# Set here the file number to be processed: 0..28
fnumber="0"

WORKING_DIR="/scratch.global/large/barley"

### Keep the lists below with the same order
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/BAM-${fnumber}.txt) ### list of bam files with read groups (output of script 05)
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)  ### list of output names (just remove .bam suffix from input list)
OUTDIR="${WORKING_DIR}/BAM_dedup"

mkdir -p ${WORKING_DIR}/BAM_intervals

module load samtools_ML/1.20

samtools index \
   --csi \
   ${OUTDIR}/$INPUT\_RG.bam

module load java/openjdk-8_322
module load gatk_ML/3.8.0

EBROOTGATK="/panfs/jay/groups/9/morrellp/public/Software/GATK_ML_3.8.0"

#### Here we generate a indel realigned bam file for each sample/library

java -jar \
   $EBROOTGATK/GenomeAnalysisTK.jar \
      -T RealignerTargetCreator \
      -R ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
      -I ${OUTDIR}/$INPUT\_RG.bam \
      -o ${WORKING_DIR}/BAM_intervals/$INPUT\.intervals

java -jar \
   $EBROOTGATK/GenomeAnalysisTK.jar \
      -T IndelRealigner \
      -R ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
      -I ${OUTDIR}/$INPUT\_RG.bam \
      -targetIntervals ${WORKING_DIR}/BAM_intervals/$INPUT\.intervals \
      --consensusDeterminationModel USE_READS \
      -o ${OUTDIR}/$INPUT\_realigned.bam
