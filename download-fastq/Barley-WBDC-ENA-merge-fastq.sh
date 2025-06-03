#!/bin/bash -l
#SBATCH --time=32:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=large@umn.edu

# Developer: Luis Arge - large@umn.edu

####################################################
# Load modules

####################################################
# Set working dir
WORKING_DIR="/scratch.global/large/barley/"
cd ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/EBI/merged-fastq

####################################################
# Concatenate files

N=10
for i in $(cut -f1 ${WORKING_DIR}/EBI/fastq_files.txt | uniq); do
   (
       R1=$(grep ${i} ${WORKING_DIR}/EBI/fastq_files.txt | cut -f2 | grep "_R1_" | awk -v path="${WORKING_DIR}EBI/fastq/" '{printf path$1"\n"}')
       R2=$(grep ${i} ${WORKING_DIR}/EBI/fastq_files.txt | cut -f2 | grep "_R2_" | awk -v path="${WORKING_DIR}EBI/fastq/" '{printf path$1"\n"}')
       cat ${R1} > ${WORKING_DIR}EBI/merged-fastq/${i}_R1.fastq.gz 2>> ${WORKING_DIR}EBI/cat.log
       cat ${R2} > ${WORKING_DIR}EBI/merged-fastq/${i}_R2.fastq.gz 2>> ${WORKING_DIR}EBI/cat.log
       echo ${i}
   ) &
   if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait -n 
   fi
done
