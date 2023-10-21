#!/bin/bash
for i in {1..23}
do
sbatch /accounts/campus/omer_ronen/projects/epiTree/scripts/analysis.sh $i
done
