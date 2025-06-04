#!/bin/bash

WORKING_DIR="/scratch.global/large/barley"
cd ${WORKING_DIR}/
mkdir -p ${WORKING_DIR}/text_lists
mkdir -p ${WORKING_DIR}/log/{map,remove_dup,add_RG,gatk_index,dedup,realign_indel}

# R1 fastq list - input
find ${WORKING_DIR}/EBI/merged-fastq/ -name *.gz | \
   grep "R1" > ${WORKING_DIR}/text_lists/R1-merged.txt

# R2 fastq list - input
find ${WORKING_DIR}/EBI/merged-fastq/ -name *.gz | \
   grep "R2" > ${WORKING_DIR}/text_lists/R2-merged.txt

mkdir -p ${WORKING_DIR}/trimmed-fastq
# R1 trimmed list - output
sed -e 's/merged-fastq/trimmed-fastq/' \
    -e 's/.fastq.gz//' \
       ${WORKING_DIR}/text_lists/R1-merged.txt \
          > ${WORKING_DIR}/text_lists/R1-trimmed.txt
# R2 trimmed list - output
sed -e 's/merged-fastq/trimmed-fastq/' \
    -e 's/.fastq.gz//' \
       ${WORKING_DIR}/text_lists/R2-merged.txt \
          > ${WORKING_DIR}/text_lists/R2-trimmed.txt
# Output file for fastp reports - json and html
sed -e 's/trimmed-fastq/fastp-report/' \
    -e 's/R1/fastp/' \
       ${WORKING_DIR}/text_lists/R1-trimmed.txt \
       > ${WORKING_DIR}/text_lists/fastp-report.txt

# Split BAM files by 10 samples
split \
   ${WORKING_DIR}/text_lists/BAM.txt \
   -l 10 \
   -d ${WORKING_DIR}/text_lists/splitted-files/BAM- \
   --additional-suffix=.txt

for i in 1 2; do
   split \
      ${WORKING_DIR}/text_lists/R${i}-trimmed.txt \
      -l 10 \
      -d ${WORKING_DIR}/text_lists/splitted-files/R${i}- \
      --additional-suffix=-trimmed.txt
done

# Remove zero from the file name
for i in $(seq 0 9); do
   mv ${WORKING_DIR}/text_lists/splitted-files/BAM-0${i}.txt ${WORKING_DIR}/text_lists/splitted-files/BAM-${i}.txt
   mv ${WORKING_DIR}/text_lists/splitted-files/splitted-files/R1-0${i}-trimmed.txt ${WORKING_DIR}/text_lists/splitted-files/splitted-files/R1-${i}-trimmed.txt
   mv ${WORKING_DIR}/text_lists/splitted-files/splitted-files/R2-0${i}-trimmed.txt ${WORKING_DIR}/text_lists/splitted-files/splitted-files/R2-${i}-trimmed.txt
done
