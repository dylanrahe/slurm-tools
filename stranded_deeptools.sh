
#!/bin/bash

###############################################################################
###  This will be the script to create normalized, stranded bigwigs from    ###
###  aligned bams					                    ###
###############################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=RNAseq
#SBATCH --output=%x_%j_%a.log
##SBATCH --array=189-205,279,282,285,288

#load the necessary modules on hcp

module purge
module load samtools/intel/1.12
module load deeptools/3.5.0

set -e

bam=$1
bamname=$(basename $bam)

#forward strand only bigwig
bamCoverage -b $bam -o ${bamname%bam}fwd.bw -bs 1 --normalizeUsing CPM --samFlagExclude 16 --ignoreDuplicates;

#reverse strand only bigwig
bamCoverage -b $bam -o ${bamname%bam}rev.bw -bs 1 --normalizeUsing CPM --samFlagInclude 16 --ignoreDuplicates
