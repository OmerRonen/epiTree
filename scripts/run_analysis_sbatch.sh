#!/bin/bash
for i in {1..22}
do
sbatch /accounts/campus/omer_ronen/projects/epiTree/scripts/analysis.sh Systemic_lupus_erythematosus Whole_Blood $i
done
