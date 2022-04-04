#!/bin/bash
 
# get count of sample files in this directory

NUMFILES=$(cat /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/sample_list.txt | wc -l)

# subtract 1 as we have to use zero-based indexing (first element is 0)
ZBNUMFILES=$(($NUMFILES - 1))

# submit array of jobs to SLURM
if [ $ZBNUMFILES -ge 0 ]; then
  sbatch --array=0-$ZBNUMFILES rna_array.sh
else
  echo "No jobs to submit, since no input files in this directory."
fi
