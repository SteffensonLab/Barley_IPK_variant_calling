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
module load parallel

####################################################
# Set working dir
WORKING_DIR="/scratch.global/large/barley/"
cd ${WORKING_DIR}/
mkdir -p ${WORKING_DIR}/EBI/{Accession_TSVs,fastq}

####################################################
# Get project and its accessions
wget -O ${WORKING_DIR}/EBI/PRJEB80165.tsv \
   'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB80165&result=analysis&fields=study_accession,sample_accession,analysis_accession,analysis_type,tax_id,scientific_name,submitted_ftp&format=tsv&download=true&limit=0'

####################################################
# Get FTP links for each FASTQ and MD5 files to be obtained
for i in $(cut -f2 ${WORKING_DIR}/EBI/PRJEB80165.tsv | sed '/^[[:space:]]*$/d' | grep -v 'sample_accession'| sed 's/;/ /g'); do
   accession=$(echo ${i})
   wget -O ${WORKING_DIR}/EBI/Accession_TSVs/${accession}.tsv \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=study_accession,sample_accession,tax_id,scientific_name,fastq_md5,fastq_ftp,submitted_md5,submitted_ftp,bam_ftp,bam_md5&format=tsv&download=true&limit=0"
done

# Look for empty file and try again for each one
for i in $(find ${WORKING_DIR}/EBI/Accession_TSVs/ -type f -empty); do
   accession=$(basename ${i} | sed 's/.tsv//')
   wget -O ${WORKING_DIR}/EBI/Accession_TSVs/${accession}.tsv \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=study_accession,sample_accession,tax_id,scientific_name,fastq_md5,fastq_ftp,submitted_md5,submitted_ftp,bam_ftp,bam_md5&format=tsv&download=true&limit=0"
done

####################################################
# Get fastq.gz files and check md5 hashes
cat ${WORKING_DIR}/EBI/Accession_TSVs/SAMEA*.tsv | \
   grep -v "submitted_md5" | \
   cut -f2,7,8 | \
   sed 's/;/\t/g' | \
   awk '{printf $1"\t"$4"\t"$2"\n"$1"\t"$5"\t"$3"\n"}' | \
   awk '{printf $1"\t"$2"\t"$3"\t"$2"\n"}' | \
   sed 's/ftp.sra.ebi.ac.uk\/vol1\/run\/ERR...\/ERR........\///' \
      > ${WORKING_DIR}/EBI/fastq_files.txt

####################################################
# Download fastq.gz files
parallel -j 100 wget -q -P ${WORKING_DIR}/EBI/fastq/ {} ::: $(cut -f4 ${WORKING_DIR}/EBI/fastq_files.txt)

####################################################
# Get md5 hashes for each downloade file
cd ${WORKING_DIR}/EBI/fastq/

rm -f ../md5-hashes_downloaded.txt
M=10
for i in *.fastq.gz; do
   (
      md5sum ${i} >> ../md5-hashes_downloaded.txt
   ) &
   if [[ $(jobs -r -p | wc -l) -ge $M ]]; then
      wait -n 
   fi
done

####################################################
# Compare md5 hashes
cd ../../

paste <(cut -f2,3 ${WORKING_DIR}/EBI/fastq_files.txt | sort -k1 | awk '{printf $1" "$2"\n"}') \
     <(sort -k2 ${WORKING_DIR}/EBI/md5-hashes_downloaded.txt | awk '{printf $2" "$1"\n"}') | \
        awk '{if ($2!=$4) printf $0"\n"}' > fastq_diff.txt
