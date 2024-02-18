#!/bin/bash

###############################################################################
### This will download fastq files and gzip them from GEO using a text      ###
### file containing the accession numbers.                                  ###
###############################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=fasterq-dump
#SBATCH --output=%x.log

fflag=0

usage="
Usage:

This will download fastq files and gzip them from GEO using a text file containing the accession numbers.

$(basename "$0") [-h] -f <accession file>

where:
  -h show this help text
  -f file containing one accession number per line"

while getopts :h:f: opt; do
  case ${opt} in
    h )
      echo "$usage"
      exit
      ;;
    f )
      file=$OPTARG ;
      fflag=1
      echo "Processing Accession numbers in file: $file" >&2
      ;;
    \? )
      printf "invalid option: -%s\n" "$OPTARG" >&2
      exit 1
        ;;
  esac
done

if [ $fflag -ne 1 ];
then
  echo "must use -f to specify accesion list" >&2
  echo "$usage"
  exit 1
fi

#Execute fasterq-dump to retreive files

/home/dpr274/bin/sratoolkit.2.10.8-centos_linux64/bin/fasterq-dump $(cat $file) -S -p

for i in *.fastq;
  do gzip $i --verbose;
done
