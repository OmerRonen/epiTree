#!/bin/bash
#SBATCH --job-name=epitree_irf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=omer_ronen@berkeley.edu
#SBATCH -o epitree_irf.out #File to which standard out will be written
#SBATCH -p yugroup

PHENO="Multiple_sclerosis"
TISSUE="Brain_Cortex"
PY="/scratch/users/omer_ronen/mutemb/bin/python"
PTH="/accounts/campus/omer_ronen/projects/epiTree"
WEI="$PTH/data/ctimp_$TISSUE.db"
RESULTS=$PTH"/results/"$PHENO"/"$TISSUE
EXPR=$PTH"/results/expression/"$TISSUE"_"$PHENO


Rscript scripts/analysis_iRF_gene.R --predx $EXPR --db $WEI --pheno data/$PHENO/pheno.csv --out $RESULTS
