#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-user=large@umn.edu

WORKING_DIR="/scratch.global/large/barley"

### Here we just create the bwa index of the reference genome, needed for mapping
### Keep the output of this command in the same dir where you keep the reference genome fasta

module load bwa/0.7.17

bwa \
   index \
      -a bwtsw \
      ${WORKING_DIR}/genome/ensembl-plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa
