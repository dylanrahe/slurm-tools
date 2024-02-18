#!/bin/bash

###############################################################################
###This will be the script to align create normalized bigwigs from          ###
###single-end or paired-end RNAseq fast files using the STAR aligner. Make  ###
###sure to EDIT ARRAY NUMBERS IN THIS SCRIPT BEFORE USING!!!                ###
###############################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=RNAseq
#SBATCH --output=%x_%j_%a.log
#SBATCH --array=189-205,279,282,285,288

#load the necessary modules on hcp

module purge
module load samtools/intel/1.11
module load star/intel/2.7.6a
module load deeptools/3.5.0

set -e


#setup arguments for the script. Include usage and an argument for help,
#an option to force starting the script at the begninning regardless of whether
#files have been created in a previous run, and options for the required file
#prefix and index files.

iflag=0
pflag=0
dflag=0
oflag=0
peflag=0
libpflag=0
skip=1

usage="
Usage:

Align fastq.gz reads (single or paired end) to specified genome and generate read-count normalized bigwig:
$(basename "$0") [-hF] -i <index directory> -r <1|2> -p <run prefix> -l <library prefix> -d <run directory> -o <output directory>.

*Note* Make sure to edit array numbers appropriately before running this script

where:
  -h show this help text
  -i index directory
  -f force reprocessing of reads even if intermediate files exist
  -r 1 for single end, 2 for paired end
  -p run prefix
  -l library name prefix
  -d run directory
  -o output directory

  Note: fastq names should be <run prefix>_n01_<lib_prefix>.fastq.gz, where "n01" represents the read number.
  For example: for the files, \"HGLJYAFX2_n01_dpr167.fastq.gz\" and \"HGLJYAFX2_n02_dpr167.fastq.gz\", the run prefix will be \"HGLJYAFX2\" and the lib prefix will be \"dpr\""

while getopts :hfi:r:p:l:d:o: opt; do
  case ${opt} in
    h )
      echo "$usage"
      exit
      ;;
    f )
      echo "Forcing re-processing of all steps" >&2
      skip=0
      ;;
    i )
      index=$OPTARG >&2 ; iflag=1
      ;;
    r )
      pe=$OPTARG >&2 ; peflag=1
      ;;
    p )
      prefix=$OPTARG >&2 ; pflag=1
      ;;
    l )
      lib_prefix=$OPTARG >&2 ; libpflag=1
      ;;
    d )
      dir=$OPTARG >&2 ; dflag=1
      ;;
    o )
      out_dir=$OPTARG >&2 ; oflag=1
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
if [ $iflag -ne 1 ];
then
  echo "must use -i to specify index directory" >&2
  echo "$usage"
  exit 1
fi
if [ $oflag -ne 1 ];
then
  echo "must use -o to specify output directory" >&2
  echo "$usage"
  exit 1
fi

if [ $pflag -ne 1 ];
then
        echo "must use -p to specify run prefix" >&2
        echo "$usage"
        exit 1
fi
if [ $dflag -ne 1 ];
then
       echo "must use -d to specify run directory" >&2
       echo "$usage"
       exit 1
fi
if [ $peflag -ne 1 ];
then
  echo "must use -r to specify index directory" >&2
  echo "$usage"
  exit 1
fi
if [ $libpflag -ne 1 ];
then
  echo "must use -l to specify library prefix name" >&2
  echo "$usage"
  exit 1
fi

#print prefix for log
echo "Starting processing for sample ${lib_prefix}${SLURM_ARRAY_TASK_ID}...";

#By default, the script will check to see if some part of the script has been run
#and skip unneccesary steps.


#Check if the alignment had been run previously, and had generated a file.
#If so, skip this step so as not to realign. If not, align with STAR using
#specified index file. If -f flag is set, force alignment
#regardless of whether alignment file already exists (and overwrite)

if [[ -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam ]] && [[ $skip = "1" ]] ; then
  echo "File ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam already exists, skipping alignment" ;
elif [[ ! -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam ]] || [[ $skip = "0" ]] ; then
  if [[ $pe -eq 1 ]] ; then
    echo "Starting single-end alignment with STAR..." ;
    STAR --runThreadN 10 --genomeDir $index --readFilesIn ${dir}/${prefix}_n01_${lib_prefix}${SLURM_ARRAY_TASK_ID}.fastq.gz --outFileNamePrefix ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat;
  elif [[ $pe -eq 2 ]] ; then
    echo "Starting paired-end alignment with STAR...";
    STAR --runThreadN 10 --genomeDir $index --readFilesIn ${dir}/${prefix}_n01_${lib_prefix}${SLURM_ARRAY_TASK_ID}.fastq.gz ${dir}/${prefix}_n02_${lib_prefix}${SLURM_ARRAY_TASK_ID}.fastq.gz --outFileNamePrefix ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat;
  fi
fi

#now check for previously generated index as before. If doesn't exist, index the bam file,
#or if -f flag is set, force reindexing (and overwrite)

if [[ -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam.bai ]] && [[ $skip = "1" ]] ; then
  echo "File ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam.bai already exists, skipping indexing";
elif [[ ! -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam.bai ]] || [[ $skip = "0" ]] ; then
  echo "Starting indexing with samtools...";
  samtools index ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam ;
fi

#Use Deeptools to generate normalized bigwigs instead of bedtools

if [[ -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.bw ]] && [[ $skip = "1" ]] ; then
  echo "File ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.bw already exists, nothing left to do!";
elif [[ ! -f ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.bw ]] || [[ $skip = "0" ]] ; then
  if [[ $pe -eq 2 ]] ; then
    echo "Creating CPM-normalized bigwig for paired-end reads with deeptools...";
    bamCoverage -b ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam -o ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}_eRiD_.bw -bs 1 --normalizeUsing CPM --extendReads --ignoreDuplicates;
  elif [[ $pe -eq 1 ]] ; then
    echo "Creating CPM-normalized bigwig for single-end reads with deeptools...";
    bamCoverage -b ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam -o ${out_dir}/${lib_prefix}${SLURM_ARRAY_TASK_ID}_iD_.bw -bs 1 --normalizeUsing CPM --ignoreDuplicates;
  fi
fi

echo "Done!";
exit 0
