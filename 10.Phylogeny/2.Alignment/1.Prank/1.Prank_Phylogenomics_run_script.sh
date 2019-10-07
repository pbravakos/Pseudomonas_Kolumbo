#!/bin/bash

<< ////
        This script calls another one, which runs prank for each gene.
        
        The genes names should be checked every time before 
        calling the script!! 
////


# First, we declare an array with all the gene names.
declare -a arr=("POG09090273" "POG090903FS"
		"POG090902FY" "POG0909018I"
		"POG090901ND" "POG090900GE"
		"POG090901RJ" "POG090902JS"
		"POG090901XJ" "POG0909006H"
		"POG090903B5" "POG090901S2"
		)




for i in "${arr[@]}"
do
    sbatch 2.Prank_Phylogenomics.sh $i
    sleep 1
done

