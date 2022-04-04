#!/bin/bash -l
#SBATCH --nodes=1  --ntasks-per-node=1 --cpus-per-task=8 --mem-per-cpu=3500M
#SBATCH --time=3-00:00:00
#SBATCH --job-name rna_align_$SLURM_ARRAY_TASK_ID
#SBATCH --output=%A.%a.out              # The output file name: <job_name>.<job_id>.out
#SBATCH --error=%A.%a.err               # The error file name: <job_name>.<job_id>.err
#SBATCH --account=cpu-s2-matocqlab-0            # The account to charge
#SBATCH --partition=cpu-s2-core-0               # The parition
 


module load unr-rc

cd /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/alignments

echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
 
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# grab our filename from a directory listing
NAMES=($(cat /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/sample_names.txt))
FILENAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
echo "My sample is ${FILENAME}"

# make new directory, change into it, and run
mkdir ${FILENAME}_align
cd ${FILENAME}_align
hisat2 -q --phred33 --no-temp-splicesite --no-mixed --no-discordant --max-intronlen 150000 \
--rna-strandness RF --no-unal -p 8 \
-x /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/genome/Nbryanti_db \
-1 /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/trimmed/${FILENAME}_R1_001_val_1.fq.gz \
-2 /data/gpfs/assoc/matocqlab/Neotoma_transcriptomics/trimmed/${FILENAME}_R2_001_val_2.fq.gz \
-S ${FILENAME}_align.sam

/data/gpfs/home/mholding/software/samtools-1.11/samtools view -bh -@ 8 -f 3 -F 256 -q 40 ${FILENAME}_align.sam | \
/data/gpfs/home/mholding/software/samtools-1.11/samtools sort -@ 8 > ${FILENAME}_align_sorted.bam
/data/gpfs/home/mholding/software/samtools-1.11/samtools index -b ${FILENAME}_align_sorted.bam
rm ${FILENAME}_align.sam
