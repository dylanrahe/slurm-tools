**Slurm Tools**

This series of tools was created to analyze genomic data using the Slurm workload manager in conjunction with the Lmod environment modules system

These tools were written in 2019-2021 and may need to be updated to refelct current versions, etc.

Tools:

**ChIPseq.sh**

This tool is used for aligning, sorting, indexing, and creating normilex bigWig files for single-end ChIP-seq reads.

Usage:

Use bwa-mem2 to align, samtools to process, and deeptools to create bigwigs for single-end read ChIPseq data.
$0 [-h] -i <indexed fasta file> -r <reads> -o <output directory>.

where:
  -h show this help text
  -i bwa-mem2 index prefix
  -r fastq read input
  -o output directory

**RNAseq.sh**

This tool is used for analyzing RNAseq data, whether single or paired end. For each run, array indexes need to be changed to be able to use this script.

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
  For example: for the files, \"HGLJYAFX2_n01_dpr167.fastq.gz\" and \"HGLJYAFX2_n02_dpr167.fastq.gz\", the run prefix will be \"HGLJYAFX2\" and the lib prefix will be \"dpr\"


**STAR_generate_index**

This tool is used to generate index files for using the STAR aligner, which is used in the RNAseq script. 

Usage:

Generate index using STAR aligner.
$(basename "$0") [-h] -d <index directory> -f <fasta files> -g <gtf annotations>.

where:
  -h show this help text

  -d directory to generate indexes in
  -f fasta file(s) from which to generate indexes
  -g gtf annotations file


**bed_cov.sh**

For a given bed file containing genomic regions, this script will use bedtools to determine the coverage of specified regions in a bam file

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


**fasterq-dump.sh**

This script uses the sratoolkit to download fastq files and gzip them from GEO

Usage:

This will download fastq files and gzip them from GEO using a text file containing the accession numbers.

$(basename "$0") [-h] -f <accession file>

where:
  -h show this help text
  -f file containing one accession number per line



**stranded_deeptools**

This script takes a single argument, a bam file, and outputs two bigWigs, one for each strand


