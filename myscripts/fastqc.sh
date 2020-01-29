#!/bin/bash

# moving the place the script will be run
cd ~/EcologicalGenomics/myresults/

# making a new directory to store results
mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XGL*fastq.gz

do 

fastqc ${file} -o fastqc/

done



