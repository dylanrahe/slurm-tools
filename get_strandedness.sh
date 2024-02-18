#!/bin/bash

#start with input file containing paths to target bams

#This begins a csv file, samples x strandedness
echo ",Forward,Reverse" > stranded_table.csv

for i in $(cat file_list.txt); do
  #Get only mapped reads
  samtools view -bF 4 $i > ./temp.bam;
  #Get sample name (wihtout path)
  echo "${i##*/}" > temp1;
  #Get the numnber of reads where the first pair (64) is on the forward strand (no 16 flag)
  samtools view -bcf 64 -F 16 ./temp.bam > temp2;
  #Get the number of reads where the first pair (64) is on the reverse strand (16): 64+16=80
  samtools view -bcf 80 ./temp.bam > temp3;
  #put the values together, separate by a comma, and append to table
  paste -d "," temp1 temp2 temp3 >> stranded_table.csv;
  #remove temporary files
  rm -f temp{.bam,1,2,3};
done
