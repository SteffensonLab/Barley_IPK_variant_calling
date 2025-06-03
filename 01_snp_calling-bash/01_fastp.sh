#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-281

WORKING_DIR="/scratch.global/large/barley"

### THESE LISTS NEED TO FOLLOW THE SAME ORDER
INPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/R1-merged.txt)  ### A list of your R1 fastq files
INPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/R2-merged.txt)  ### A list of your R2 fastq files
OUTPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/R1-trimmed.txt) ### A list of R1 output names -- just capture the meaningful part of the fastq names including R1 (remove the fq.gz or fq suffix)
OUTPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/R2-trimmed.txt) ### A list of R2 output names -- just capture the meaningful part of the fastq names including R2 (remove the fq.gz or fq suffix)
JSON_HTML=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/fastp-report.txt) ### Fastp report file name

module load fastp_ML

### Trimming --  for each sample pair of raw fastq reads or for each library, we produce a pair of trimmed output files.
### Remember to keep R1 and R2 in the output names created in lists 3 and 4 above
### This below uses 4 cores

fastp \
   -w 2 \
   -i $INPUT1 \
   -I $INPUT2 \
   -o $OUTPUT1\_trimmed.fastq.gz \
   -O $OUTPUT2\_trimmed.fastq.gz \
   -j $JSON_HTML\.json \
   -h $JSON_HTML\.html
