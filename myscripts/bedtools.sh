#!/bin/bash

/data/popgen/bedtools2/bin/bedtools closest -a /users/a/p/ap1/EcologicalGenomics/myresults/bedfiles_copepod/diffmeth.bed\
      -b /data/project_data/epigenetics/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.bed \
      -D b | \
      awk '!($10=="-")' > hits.bed # save filename is hits.bed