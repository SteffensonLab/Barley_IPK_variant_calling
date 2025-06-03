#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-10
#SBATCH --output=/scratch.global/large/barley/log/add_RG/05_add_RG-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/add_RG/05_add_RG-%x.%j.err

# Set here the file number to be processed: 0..28
fnumber="0"

WORKING_DIR="/scratch.global/large/barley"

### Keep the lists below with the same order
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/BAM-${fnumber}.txt) ### List of input deduplicated bam files
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt) ### List of output names (just remove .bam) from the inputs
#NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list3.txt)   ### Here you want to extract the sample name from the input name, which is used to set the read IDs. Can be the same as list2.txt
NAME=$(echo ${INPUT} | sed -E 's/\/.*\/|_dedup.bam//g')
OUTDIR="${WORKING_DIR}/BAM_dedup"

module load picard/2.25.6
module load java/openjdk-8_202

EBROOTPICARD="/common/software/install/migrated/picard/2.25.6"

### Here we add read groups, we start with our deduplicated bam files and we get a deduplicated bam with read groups assigned per sample/library

java -jar \
   $EBROOTPICARD/picard.jar \
      AddOrReplaceReadGroups \
      I=${OUTDIR}/$INPUT\_dedup.bam \
      O=${OUTDIR}/$INPUT\_RG.bam \
      RGID=$NAME \
      RGLB=$NAME\_LB \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=$NAME && rm -f ${OUTDIR}/$INPUT\_dedup.bam
