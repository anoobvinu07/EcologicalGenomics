#!/bin/bash

# moving the place the script will be run
cd ~/EcologicalGenomics/myresults/

# making a new directory to store results
mkdir fastqcT1

for file in /data/project_data/RS_RNASeq/fastq/XBM*D*fastq.gz

do 

fastqc ${file} -o fastqcT1/

done



for file in /data/project_data/RS_RNASeq/fastq/XBM*H*fastq.gz

do

fastqc ${file} -o fastqcT1/

done

