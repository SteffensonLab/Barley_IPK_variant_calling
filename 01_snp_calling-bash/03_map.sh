#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-2
#SBATCH --output=/scratch.global/large/barley/log/map/03_map-%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/map/03_map-%x.%j.err

# Set here the file number to be processed: 0..28
fnumber="0"

WORKING_DIR="/scratch.global/large/barley"
mkdir -p ${WORKING_DIR}/bwa_output

#### All these lists below need to follow same order
INPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/R1-${fnumber}-trimmed.txt) ### List of trimmed R1 fastq reads (produced with script 01)
INPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/R2-${fnumber}-trimmed.txt) ### List of trimmed R2 fastq reads (produced with script 01)
OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/text_lists/splitted-files/BAM-${fnumber}.txt) ### List of output names -- extract the meaningful part of the name from the trimmed reads names. There is a single output file for each pair of reads (each sample or library)

module load bwa/0.7.17
module load samtools_ML/1.20

#### Here we map the trimmed reads to the ref genome and we produce a sam, then we convert it to bam, we sort it and finally we index it
#### We end up with 1 bam file per sample after this. If you had multiple libraries per sample, you'd end up with 1 bam per library


bwa mem \
   -t 50 \
   ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
   $INPUT1\_trimmed.fastq.gz \
   $INPUT2\_trimmed.fastq.gz  \
      > ${WORKING_DIR}/bwa_output/$OUTPUT\.sam && rm -f $INPUT1\_trimmed.fastq.gz $INPUT2\_trimmed.fastq.gz

samtools view \
   -Sb \
   -q 10 \
   ${WORKING_DIR}/bwa_output/$OUTPUT\.sam \
      > ${WORKING_DIR}/bwa_output/$OUTPUT\.bam && rm ${WORKING_DIR}/bwa_output/$OUTPUT\.sam

samtools sort \
   --threads 50 \
   ${WORKING_DIR}/bwa_output/$OUTPUT\.bam \
      > ${WORKING_DIR}/bwa_output/$OUTPUT\_sorted.bam && rm ${WORKING_DIR}/bwa_output/$OUTPUT\.bam

samtools index \
   --csi \
   ${WORKING_DIR}/bwa_output/$OUTPUT\_sorted.bam

