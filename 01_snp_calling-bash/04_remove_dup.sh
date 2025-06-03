#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-10
#SBATCH --output=/scratch.global/large/barley/log/dedup/04_remove_dup-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/dedup/04_remove_dup-%x.%j.err

# Set here the file number to be processed: 0..28
fnumber="0"

WORKING_DIR="/scratch.global/large/barley"

### The lists below need to follow the same order
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/BAM-${fnumber}.txt) ### list of input bam files (output of script 03)
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/BAM_dedup.txt) ### list of output names --  I usually just remove the suffix (.bam) to extract the name from the input, and then the command below will add _dedup.bam
OUTDIR="${WORKING_DIR}/BAM_dedup"

mkdir -p ${OUTDIR}
mkdir -p ${WORKING_DIR}/dup_metrics

module load picard/2.25.6
module load java/openjdk-8_202

EBROOTPICARD="/common/software/install/migrated/picard/2.25.6"

### Here we remove duplicates. We feed it a bam, and we get a deduplicated bam per sample/library

java -jar $EBROOTPICARD/picard.jar \
   MarkDuplicates \
      INPUT=${WORKING_DIR}/bwa_output/$INPUT\_sorted.bam \
      OUTPUT=${OUTDIR}/$INPUT\_dedup.bam \
      METRICS_FILE=${WORKING_DIR}/dup_metrics/$INPUT\_DUP_metrics.txt \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=true && rm -f ${WORKING_DIR}/bwa_output/$INPUT\_sorted.bam
