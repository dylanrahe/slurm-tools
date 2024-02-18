#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=60GB
#SBATCH --job-name=STAR_generate_index
#SBATCH --output=%x%j.log

module purge
module load star/intel/2.7.6a

set -e

dflag=0
fflag=0
gflag=0

usage="
Usage:

Generate index using STAR aligner.
$(basename "$0") [-h] -d <index directory> -f <fasta files> -g <gtf annotations>.

where:
  -h show this help text

  -d directory to generate indexes in
  -f fasta file(s) from which to generate indexes
  -g gtf annotations file"

while getopts :hd:f:g: opt; do
  case ${opt} in
    h )
      echo "$usage"
      exit
      ;;
    d )
      gdir=$OPTARG >&2 ; dflag=1
      ;;
    f )
      fasta=$OPTARG >&2 ; fflag=1
      ;;
    g )
      gtf=$OPTARG >&2 ; gflag=1
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
if [ $dflag -ne 1 ];
then
  echo "must use -d to specify genome directory" >&2
  echo "$usage"
  exit 1
fi

if [ $fflag -ne 1 ];
then
        echo "must use -f to specify fasta file(s)" >&2
        echo "$usage"
        exit 1
fi
if [ $gflag -ne 1 ];
then
       echo "must use -g to specify gtf file" >&2
       echo "$usage"
       exit 1
fi

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $gdir --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 74

exit 0

