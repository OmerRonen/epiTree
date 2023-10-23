#!/bin/bash
for i in {1..22}
do
sbatch /accounts/campus/omer_ronen/projects/epiTree/scripts/analysis.sh Multiple_sclerosis Brain_Cortex $i
done
