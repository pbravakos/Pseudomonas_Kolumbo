#!/bin/bash
<< ////
	We want to change the fasta headers of each fasta single copy gene in order to make them identical between all the genes.
	That way we will be able to concatenate the genes later on.
	Moreover, we want to remove all gaps from the alignment in order to repeat it again with other programs of our choosing.
	Anvio does the alignment with muscle by default, but we want to try other MSA programs. 

	The bellow sed command removes all gaps from the alignment and in the fasta header only the genome name of the strain remains.
	Everything else is stripped from the header.
	The first part of the pipe removes everything preceding the strain names. 
	106 is a constant number in all fasta headers and is the exact number of characters that should be removed before the strain name!
	If going to use this sed command please check that 106 applies to the new fasta headers!!!!
	The second part removes everything after the strain name till the next line character (but not the next line character itself!!)
	The third part removes all gaps from the alignment (and the '-' from some headers that contain it!).
	The final fasta file is ready for alignment by any program we want e.g. Prank.
////

# INITIAL PARAMETERS
AlignedDir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Anvio_aligned_genes
OutputDir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Genes_fasta_headers_and_gaps_stripped

# Check that the directories needed to run have already been created.
[[ ! -d $AlignedDir ]] && { echo "Directory ${AlignedDir} cannot be not found!" >&2; exit 1; }

[[ ! -d $OutputDir ]] && { echo "Directory ${OutputDir} cannot be not found!" >&2; exit 1; }


# First, we declare an array with all the gene names.

declare -a arr=("POG090903HI" "POG090902MV" 
                "POG09090004" "POG090900O9" 
                "POG09090188" "POG090901ZM" 
                "POG090900YW" "POG090901TG" 
                "POG0909039Q" "POG0909032E" 
                "POG090901TJ" "POG090901ZB"
                "POG09090095" "POG090900UZ"
                "POG0909008M" "POG090900IU"
		"POG090903GW" "POG0909008P"
		"POG09090288" "POG090901JA"
		"POG0909006H" "POG0909001N"
		"POG090900YH" "POG090903CH"
		"POG090901HE" "POG090901H1"
		"POG090900CF" "POG090900MN"
		"POG090901OY" "POG090900GC"
		"POG090903B5" "POG0909024W"
		"POG090900OV" "POG0909001O"
		"POG09090325" "POG090900ZF"
		"POG0909008D" "POG0909006C"
		"POG0909002C" "POG090901XG"
		"POG090900YI" "POG090900CU"
		"POG090900OV" "POG0909001O"
		"POG09090325" "POG090900ZF"
		"POG0909008D" "POG0909006C"
		"POG0909002C" "POG090901XG"
		"POG090900YI" "POG09090384"
		"POG0909028G" "POG090901S2"
		"POG090900RW" "POG0909000G" 
		"POG090900FV" "POG090903FS"
		"POG09090273" "POG090903FS"
		"POG090902FY" "POG0909018I"
		"POG090901ND" "POG090900GE"
		"POG090901RJ" "POG090902JS"
		"POG0909000T" "POG090903DD" 
		"POG090901L6" "POG090903LP"
		"POG090901XJ"
) 



# Next, we run the sed command inside a for loop, for each gene of the array.
for i in "${arr[@]}"
do
   sed -E 's/>.{106}/>/g' ${AlignedDir}/Pseudomonas_${i}.fasta | sed -E 's/\|\S*//g' | sed 's/-//g' > ${OutputDir}/Pseudomonas_${i}.fas
done

