#!/bin/bash

# We add the sleep command because the script writes to a file
# and needs some time in order to write to it!
while read -r p; do
  sbatch 2.Anvio_Phylogenomics_Create_database.sh $p
  sleep 2
done <"${HOME}"/Pseudomonas/Pseudomonas_database/Strains_Downloaded_from_NCBI.txt
