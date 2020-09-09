#!/bin/bash
#
#SBATCH --job-name=plink2dos
#SBATCH --time=900:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
ml python
ml biology
ml plink/1.90b5.3

chr=$1
echo $chr
EXE="./convert_plink_to_dosage.py"
BIN="path_to_plink_files/name_of_plink_files"$chr
OUT="path_to_output_plink_files/name_of_output_dosage_files_chr"$chr

srun python $EXE -b $BIN -o $OUT -p plink 
 
