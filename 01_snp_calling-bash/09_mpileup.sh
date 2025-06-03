#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nstasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=large@umn.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch.global/large/barley/log/variant_calling/%x.%j.out
#SBATCH --error=/scratch.global/large/barley/log/variant_calling/%x.%j.err

# Set here the file number to be processed: 0..28
fnumber="0"

module load bcftools_ML_3/1.21

WORKING_DIR="/scratch.global/large/barley"
mkdir -p ${WORKING_DIR}/raw-VCF/

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare-chromosomes.txt) ### list of chromosomes (can be found in the FASTA index file of reference genome -- .fai file). This species has 14, that's why array number is 14. We parallelize by chromosome.

### If you have a very fragmented reference with 1000s of scaffolds, rather than sending 1000s of jobs, you can feed a list of chromosomes to the command rather than a single chromosome/scaffold name.

### So you can for example split 1000 scaffolds names into 10 lists of 100 scaffolds each, and feed those 10 lists to the command using a job array of size 10 (1-10) -- one list of scaffolds per job.


### Here we call SNPs.

### list.txt is a list of ALL the realigned bam files of ALL samples.

### If you have multiple bam per samples due to multiple libraries, merge them into one bam per sample before this step.

### ploidymap.txt is a list that keeps the same order of samples in list.txt and should use the sample ID names you set in script 05.
### For each sample ID name has the ploidy -- for a diploid species just create a tab separated file with 2 columns (first column: sample ID name set in script 05, second columnd: 2)


bcftools mpileup \
   -Ou \
   --threads 64 \
   -f ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa \
   --bam-list ${WORKING_DIR}/text_lists/splitted-files/BAM_for_variant_calling-${fname}.txt \
   -q 5 \
   -r $CHROM \
   -I \
   -a FMT/AD,FMT/DP | \
      bcftools call \
         --threads 64 \
         -S ${WORKING_DIR}/text_lists/splitted-files/ploidy-${fname}.txt \
         -G - \
         -f GQ \
         -mv \
         -Ov \
            > ${WORKING_DIR}/raw-VCF/${CHROM}-subpop-${fname}.vcf
