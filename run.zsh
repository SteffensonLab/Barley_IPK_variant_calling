#!/bin/bash

# e.g. SNP calling on chr1H from 0Mb to 20Mb

# --excl-flags: flag could be found in https://broadinstitute.github.io/picard/explain-flags.html.3332 means exclude reads with (1) unmapped (2) not primary alignment (3) PCR or optical duplicate (4) supplementary alignment 

# -R: bed file for target region

# -f: genome reference

# -b: indexed bam or cram file list (e.g. for 5 samples: s1, s2, s3, s4, s5)

# -V indels: exclude indels, only generate SNP in output file

bcftools mpileup --threads 5 -a DP,AD -q 20 -Q 20 --excl-flags 3332 -R chr1H.part1.bed -f MorexV3.fa -b bam.list | bcftools call -V indels --threads 5 -mv -O z -o chr1H.part1.snp.vcf.gz

echo finish
