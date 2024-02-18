#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=8:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=ChIPseq
#SBATCH --output=%x_%a.log

set -e

module load bwa-mem2/2.1
module load samtools/intel/1.14
module load deeptools/3.5.0

#get options

usage="
Usage:

Use bwa-mem2 to align, samtools to process, and deeptools to create bigwigs for single-end read ChIPseq data.
$0 [-h] -i <indexed fasta file> -r <reads> -o <output directory>.

where:
  -h show this help text
  -i bwa-mem2 index prefix
  -r fastq read input
  -o output directory
  "

while getopts "hi:r:o:" options; do
  case ${options} in
    h)
      echo $usage
      exit 1
      ;;
    i)
      index=${OPTARG}
      ;;
    r)
      reads=${OPTARG}
      ;;
    o)
      outdir=${OPTARG}
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    \?)
      echo "invalid option: -%s\n" "$OPTARG" >&2
      exit 1
      ;;
  esac
done





filename=$(basename $reads);

bwa-mem2 mem -t 12 $index $reads > ${outdir}/${filename%fastq}.sam;

samtools sort ${outdir}/${filename%fastq}.sam -@ 12 -o ${outdir}/${filename%fastq}sort.bam -O BAM;

samtools index ${outdir}/${filename%fastq}sort.bam;

bamCoverage -b ${outdir}/${filename%fastq}sort.bam -o ${outdir}/${filename%fastq}_iD_.bw -p 12 -bs 1 --normalizeUsing CPM --ignoreDuplicates;

exit 1
