#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=80GB
#SBATCH --job-name=bedcov
#SBATCH --output=%x_%a.log
#SBATCH --array=395-404

#load the necessary modules on hcp

module purge
module load bedtools/intel/2.29.2
module load samtools/intel/1.12

#exit upon error

set -e

#set flags

aflag=0
bflag=0
lflag=0

#Add "help" section

usage="
Usage:

Use bedtools coverage to generate coverage information for specific bed regions
$(basename "$0") [-h] -a <bed file with regions> -b <bam file path prefix> -l <library prefix>.

where:
  -h show this help text
  -a bed file containing regions/bins to compare
  -b bam file name prefix, with path
  -l library prefix

***NOTE: Make sure to edit script to change array numbers before running!***

For example: for a bam file in /path/to/bams/ with the name HVTDGEWW1_l01_n01_mb395.fastq.sorted.bam,
  -b would be \"/path/to/bams/HVTDGEWW1_l01_n01_\" and -l would be \"mb\"
  "

#collect options

while getopts :ha:b:l: opt; do
  case ${opt} in
    h )
      echo "$usage"
      exit
      ;;
    a )
      bed_regions=$OPTARG >&2 ; aflag=1
      ;;
    b )
      bam_prefix=$OPTARG >&2 ; bflag=1
      ;;
    l )
      lib_prefix=$OPTARG >&2 ; lflag=1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    \? )
      printf "invalid option: -%s\n" "$OPTARG" >&2
      exit 1
        ;;
  esac
done

#check option flags

if [ $aflag -ne 1 ];
then
  echo "must use -a to specify bed file" >&2
  echo "$usage"
  exit 1
fi
if [ $bflag -ne 1 ];
then
  echo "must use -b to specify bam files with path prefix" >&2
  echo "$usage"
  exit 1
fi
if [ $lflag -ne 1 ];
then
  echo "must use -l to specify library prefix" >&2
  echo "$usage"
  exit 1
fi

#print prefix for log
echo "Starting processing for sample ${lib_prefix}${SLURM_ARRAY_TASK_ID}..."

#execute and generate outputs

bedtools coverage -a $bed_regions -b ${bam_prefix}${lib_prefix}${SLURM_ARRAY_TASK_ID}.fastq.sorted.bam -mean > ${lib_prefix}${SLURM_ARRAY_TASK_ID}_full.cov.txt
samtools flagstat -O tsv ${bam_prefix}${lib_prefix}${SLURM_ARRAY_TASK_ID}.fastq.sorted.bam > ${lib_prefix}${SLURM_ARRAY_TASK_ID}.stats.tsv
