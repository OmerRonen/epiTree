#!/bin/bash
#SBATCH --job-name=dos2predix
#SBATCH --time=2880:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
ml python/2.7.13
ml py-numpy
ml biology
ml plink/1.90b5.3



GENOP="path_to_dosage_files"
GENO='name_of_dosage_files_chr'
WEI="../data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db"
OUT="../data/name_of_output_prediXcan_txt_file"
FAM='name_of_fam_file.fam'
EXE="./PrediXcan.py"


srun $EXE --predict --dosages $GENOP --dosages_prefix $GENO --samples $FAM --weights $WEI --output_prefix $OUT
