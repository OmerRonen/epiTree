#!/bin/bash
#SBATCH --job-name=epitree
#SBATCH --mail-type=ALL
#SBATCH --mail-user=omer_ronen@berkeley.edu
#SBATCH -o epitree.out #File to which standard out will be written
#SBATCH -p yugroup


PY="/scratch/users/omer_ronen/mutemb/bin/python"
PTH="/accounts/campus/omer_ronen/projects/epiTree"
WEI="$PTH/data/ctimp_Brain_Cortex.db"
DOS="$PTH/results/dosage"
EXPR="$PTH/results/expression"
B_FILE="data/Multiple_sclerosis/ukbb/Multiple_sclerosis_chr"
FAM="/accounts/campus/omer_ronen/projects/epiTree/data/Multiple_sclerosis/ukbb/Multiple_sclerosis_chr"
i=$1
PRE_FIX="chr_$i"
B_FILE_i="$B_FILE$i"
FAM_i="$FAM$i.fam"
DOS_i="$DOS/chr_$i"
EXPR_i="$EXPR/chr_$i"
PREDX_i="$EXPR/chr_"$i"_predicted_expression.txt"
RESULTS_i="$PTH/results/ms/chr_"$i
# create Results directory if it doesn't exist
mkdir -p $RESULTS_i
# Execute the commands for each iteration
$PY scripts/convert_plink_to_dosage.py -b "$B_FILE_i" -o "$DOS_i"
$PY scripts/PrediXcan.py --predict --dosages "$DOS" --dosages_prefix "$PRE_FIX" --samples "$FAM_i" --weights "$WEI" --output_prefix "$EXPR_i"
Rscript scripts/analysis_iRF_gene.R --predx $PREDX_i --db $WEI --pheno data/Multiple_sclerosis/pheno.csv --out $RESULTS_i
#for i in {1..23}
#do
#  PRE_FIX="chr_$i"
#  B_FILE_i="$B_FILE$i"
#  FAM_i="$FAM$i.fam"
#  DOS_i="$DOS/chr_$i"
#  EXPR_i="$EXPR/chr_$i"
#  PREDX_i="$EXPR/chr_"$i"_predicted_expression.txt"
#  RESULTS_i="$PTH/results/ms/chr_"$i
#  # create Results directory if it doesn't exist
#  mkdir -p $RESULTS_i
#  # Execute the commands for each iteration
#  $PY scripts/convert_plink_to_dosage.py -b "$B_FILE_i" -o "$DOS_i"
#  $PY scripts/PrediXcan.py --predict --dosages "$DOS_i" --dosages_prefix "$PRE_FIX" --samples "$FAM_i" --weights "$WEI" --output_prefix "$EXPR_i"
#  Rscript scripts/analysis_iRF_gene.R --predx $PREDX_i --db $WEI --pheno data/Multiple_sclerosis/pheno.csv --out $RESULTS_i
#done