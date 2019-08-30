#PBS -N extract_snps
#PBS -o ukb_all_snps
#PBS -e ukb_all_snps
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -t 1-22

# Script for extracting individual SNPs in UK Biobank on BC
# The only input is a text file of rsID

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

module add languages/R-3.0.2
module load apps/qctool-2.0

# Gen files on BC3. You should not need to change this.
genfile="ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded"

# Change this to your directory on Blue Crystal
data="uk_biobank/"
# Point this towards your text file with RSids in it
snpfile="snplists/snps.txt"

echo ${i}

#This runs QCTOOL and extracts the SNPs
qctool  -g ${genfile}/ukb_imp_chr${i}_v2.bgen \
		-incl-rsids ${data}/${snpfile} \
		-og ${data}/extracted_snps/final_full_snps_${i}.gen

#This converts the extracted files in dta format
Rscript --vanilla \
			${data}/scripts/R-transpose.r \
			${data}/extracted_snps/final_full_snps_${i}.gen