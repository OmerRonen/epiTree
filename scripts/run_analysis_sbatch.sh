#!/bin/bash
for i in {1..22}
do
sbatch /accounts/campus/omer_ronen/projects/epiTree/scripts/analysis.sh Rheumatoid_arthritis Whole_Blood $i
done
