#!/bin/bash

/users/a/p/ap1/EcologicalGenomics/myresults/bedfiles_copepod/diffmeth.bed closest -a diffmeth.bed \
      -b /data/project_data/epigenetics/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.bed \
      -D b | \
      awk '!($10=="-")' > hits.bed 