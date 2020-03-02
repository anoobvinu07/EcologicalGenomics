#!/bin/bash

# moving the place the script will be run
cd ~/EcologicalGenomics/myresults/

# making a new directory to store results
mkdir fastqcT1_cleaned

for file in /data/project_data/RS_RNASeq/fastq/cleanreads/XBM*D*cl.fq

do 

fastqc ${file} -o fastqcT1_cleaned/

done



for file in /data/project_data/RS_RNASeq/fastq/cleanreads/XBM*H*cl.fq

do

fastqc ${file} -o fastqcT1_cleaned/ #output folder location

done

