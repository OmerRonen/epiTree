#!/bin/bash
#
#SBATCH --job-name=plinkSubs
#SBATCH --time=300:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
ml biology
ml plink/1.90b5.3


chr=$1
echo $1
BIN="path_to_plink_files/name_of_plink_files"$chr
OUT="path_to_output_plink_file/name_of_output_plink_file"$chr
FAM="path_to_fam_file/name_of_fam_file.fam"
SNP="path_to_snp_list/name_to_snp_list.snplist"

srun plink --bfile $BIN --indiv-sort f $FAM --keep $FAM --extract $SNP --make-bed --out $OUT
