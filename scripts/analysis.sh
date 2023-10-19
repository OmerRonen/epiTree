#!/bin/bash
#SBATCH --job-name=epitree
#SBATCH --mail-type=ALL
#SBATCH --mail-user=omer_ronen@berkeley.edu
#SBATCH -o epitree.out #File to which standard out will be written
#SBATCH -p yugroup
PY="/scratch/users/omer_ronen/mutemb/bin/python"
WEI="/accounts/campus/omer_ronen/projects/epiTree/data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db"
FAM='/accounts/campus/omer_ronen/projects/epiTree/epitree/pick.fam'
GENOP="/scratch/users/omer_ronen"
GENO="1"
DOS="$GENOP/$GENO"
OUT="$GENOP/pick"

#$PY scripts/convert_plink_to_dosage.py -b $1 -o $DOS
#$PY scripts/PrediXcan.py --predict --dosages $GENOP --dosages_prefix $GENO --samples $FAM --weights $WEI --output_prefix $OUT

$PY scripts/PrediXcan.py --predict --dosages /scratch/users/omer_ronen/ --dosages_prefix 1 --samples /accounts/campus/omer_ronen/projects/epiTree/epitree/pick.fam --weights /accounts/campus/omer_ronen/projects/epiTree/data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db --output_prefix $OUT
