#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="SIFT"
#SBATCH --output=SIFT_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

## Some job specific info
echo "Job ID is = " "$SLURM_JOBID"
echo "SLURM cluster name = " "$SLURM_CLUSTER_NAME"
echo "SLURM partition = " "$SLURM_JOB_PARTITION"
echo "SLURM node list = " "$SLURM_JOB_NODELIST"
echo "SLURM num of nodes = " "$SLURM_JOB_NUM_NODES"
echo "SLURM number of tasks = " "$SLURM_NTASKS"
echo "SLURM memory per node = " "$SLURM_MEM_PER_NODE"
echo "SLURM memory per cpu = " "$SLURM_MEM_PER_CPU"
echo "working directory = " "$SLURM_SUBMIT_DIR"
echo "=================================================="
echo "SBATCΗ job started " "$(date)"
echo "=================================================="
echo

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a

generalInfo () {
    cat << ////
    	
    	This script allows us to predict whether a protein sequence variation affects protein function. 
    	It should run after we have found SNPs with Mummer, in a previous step of this pipeline.
	
	To run this script two arguments (two Strains) must be given. The first argument will be used as the Reference strain and the second as the Query strain e.g.:
	sbatch $0 Strain02 Strain03
	
	IMPORTANT:
	We have prepared the RefSeq complete non redundant blast database, before running this script, and have named it "RefseqCompleteNrDB"!
	
	IMPORTANT:
	We have modified manually the file "SIFT_for_submitting_fasta_seq.csh" following installation instructions, as:
	
		#       Location of blastpgp
			setenv NCBI /home1/pbravakos/Software/ncbi-blast-2.8.1+/bin

		#       Location of SIFT
			setenv SIFT_DIR /home1/pbravakos/Software/sift6.2.1

		#       SIFT's output files are written here
			setenv tmpdir /home1/pbravakos/Pseudomonas/SIFT
	
	IMPORTANT:
	We have modified the file "seqs_chosen_via_median_info.csh" to run psiblast with the number of cores of the current script ($SLURM_NTASKS)
		
////
}

# Check that an argument has been given in the correct form.
if [[ $# -ne 2 ]]; then
    echo "Exactly two arguments should be given!" >&2
    generalInfo >&2
    exit 1
fi

if [[ ! $1 =~ ^Strain[0-9]{2}$ ]] || [[ ! $2 =~ ^Strain[0-9]{2}$ ]] ; then
    echo "Exactly two arguments should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
    generalInfo >&2
    exit 1
fi



# INITIAL PARAMETERS
StrainX1=$1
StrX1=${StrainX1/ain/}
StrainX2=$2
StrX2=${StrainX2/ain/}

SiftDir=${HOME}/Software/sift6.2.1/bin
RefSeqDb=/mnt/dbs/RefSeq-Complete/nonredundant_data/RefseqCompleteNrDB

VarTxt=variation.txt

AnvioDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX1}
anviGetSeq="anvi-get-sequences-for-gene-calls"

SedTxt=AA_3letter_to1letter

SNPDir=${HOME}/Pseudomonas/Mummer/SNP_discovery/${StrX1}_vs_${StrX2}
SNPFile=${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv


#---------------------------------------------------------------------------------------------------------------------------------
# Check for the prerequisites 
if [[ ! -e ${SNPDir}/${SNPFile} ]]; then
    echo "${SNPFile} cannot be found in ${SNPDir}. Please run the Mummer SNP analysis first!" >&2
    generalInfo >&2
    exit 1
fi
    
#------------------------------------------------------------------------------------------------------------------------------------

# Start a job dependency to move the Sdtout file to the correct folders
sbatch --dependency=afterany:"$SLURM_JOB_ID" Move_Slurm_output_files.sh "$SLURM_JOB_NAME" "$SLURM_JOB_ID" "${StrX1}_vs_${StrX2}"

echo
echo " 					START THE ANALYSIS!!!!!"
echo
echo


mkdir -p ${StrX1}_vs_${StrX2}
cd ${StrX1}_vs_${StrX2}
ln -s ${SNPDir}/${SNPFile} .

# First, create the sed file to turn all three letter amino acid codes to one letter codes
cat > ${SedTxt} <<"EOF"
s/Thr/T/g
s/Ala/A/g
s/Tyr/Y/g
s/His/H/g
s/Gln/Q/g
s/Asn/N/g
s/Lys/K/g
s/Asp/D/g
s/Glu/E/g
s/Cys/C/g
s/Trp/W/g
s/Arg/R/g
s/Ser/S/g
s/Gly/G/g
s/Phe/F/g
s/Leu/L/g
s/Ile/I/g
s/Val/V/g
s/Pro/P/g
s/Met/M/g
s/Stop/*/g
EOF

# Next, we want to extract from the SNP file we created in a previous part of this pipeline, only the non redundant and useful information.
# We add 0.99 to second column of the first awk command, because we want to round up, whenever there is a decimal number.
awk 'BEGIN{FS="\t"} !seen[$1,$5,$10]++ && NR!=1 {printf("%d\t%d\t%s\t%s\n", $1,$5/3+0.99,$7,$10)}' ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv | awk 'BEGIN{FS=OFS="\t"; print "Gene caller ID", "AA position", "AA", "AA in query"} !seen[$1,$2,$3,$4]++ && $3!=$4 {print}' | sed -f ${SedTxt} > AA_${StrX1}_vs_${StrX2}_SNPs.tsv


# Next, we try to extract all the amino acid fasta files, for each gene call with a SNP.
while read GeneID
do
    $anviGetSeq -c ${AnvioDir}/${StrX1}_Pilon_contigs.db --gene-caller-ids ${GeneID} --get-aa-sequences -o ${StrX1}_${GeneID}.faa
done < <(awk 'BEGIN{FS="\t"} NR>1 {print $1}' AA_${StrX1}_vs_${StrX2}_SNPs.tsv | sort -nu)


# First delcare an associative array (like the dictionary in Python!)
declare -A GeneIDArray

# Next, we will fill the associative array with the Gene call IDs as keys, and values the corresponding variations. (Not to be confused with the array of the while command!)
# while IFS=$'\t' read -r -a Array
while read -r -a Array
do
   if [[ -v "GeneIDArray[${Array[0]}]" ]]; then
       GeneIDArray["${Array[0]}"]+=" ${Array[2]}${Array[1]}${Array[3]}"
   else
       GeneIDArray["${Array[0]}"]="${Array[2]}${Array[1]}${Array[3]}"
   fi    
done < <(awk 'BEGIN{FS="\t"} NR>1 {print $1,$2,$3,$4}' AA_${StrX1}_vs_${StrX2}_SNPs.tsv | sort -n)

# Finally, we will run SIFT, after creating the variation file.
for GenId in "${!GeneIDArray[@]}" ; do
    echo "${GeneIDArray[${GenId}]}" | sed 's/ /\n/g' > GeneID_${GenId}_${VarTxt}
    ${SiftDir}/SIFT_for_submitting_fasta_seq.csh ${StrX1}_${GenId}.faa ${RefSeqDb} GeneID_${GenId}_${VarTxt}  
done

# Move output file to working directory (this only works if we have set up in the SIFT configuration file the output file to be there!)
mv ../${StrX1}_*.SIFTprediction .
rename "s/${StrX1}/${StrX1}_vs_${StrX2}/g" *.SIFTprediction 
#mv ../${StrX1}_*.faa.query.out
rm ../${StrX1}_*.* ${StrX1}_*.faa AA_3letter_to1letter



<<"HERE"

	NEXT STEPS:
	1) When we have run all Strain combinations (and have created the respective files), we can run the following commands in the working directory:
	mkdir Results && cp Str*/St*.SIFTprediction Results/ && cd Results/
	grep "DEL" * | sed -E 's/\./:/g;s/:/\t/g;s/([[:digit:]]{2})_([[:digit:]]+)/\1\t\2/g;s/^Str([[:digit:]]+)_vs_Str([[:digit:]]+)/Strain\1_vs_Strain\2/g' | awk 'BEGIN{FS=OFS="\t"; print "Strain_Comparison", "Gene_Caller_ID", "AA_Substitution", "SIFT_prediction"} {print $1, $2,$4, $5}' > Deleterious.tsv
	grep "WAR" * | sed -E 's/\./:/g;s/:/ /g;s/([[:digit:]]{2})_([[:digit:]]+)/\1 \2/g;s/^Str([[:digit:]]+)_vs_Str[[:digit:]]+/Strain\1/g' | awk 'BEGIN{OFS="\t";print "Strain", "Gene_Caller_ID", "AA_Position", "SIFT_prediction"} {print $1, $2,$5, $4" not allowed!"}' > Warnings.tsv
	rm Str*
	
	As a result of the above commands we will have two new files in the Results directory that summarize the output of SIFT.
	The most important file is "Deleterious.tsv" which summarizes all the Deleterious substitutions found by SIFT.
	
	2) Next,  manually add a new column named "SIFT_predictions" on each Mummer output file, based on the "Deleterious.tsv" file, which will add the word "DELETERIOUS" to each row with a deleterious amino acid variant.
	
	3) Then, we go to the Mummer Directory (where Mummer was run in the previous part of this pipeline) create the file with all the Deleterious :
	mkdir Results && cp Str*/*_final.tsv Results && cd Results
	grep "DELETER" Str* | sed -E 's/^Str([[:digit:]]+)_versus_Str([[:digit:]]+)_SNPs_AA_final.tsv:/Strain\1_vs_Strain\2\t/g' | awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4, $5, $8, $11, $NF}' | sort -t$'\t' -k3,5 | awk 'BEGIN{FS=OFS="\t"; print "Reference_vs_Query_Strains", "Reference_Gene_Caller_ID", "Source", "Accession", "Annotation", "Reference_Amino_Acid", "Query_Amino_Acid", "SIFT_prediction"} {print}' > Deleterious_AA_sub.tsv
	rm Str*	
	
	4) Then we will try to find which of these annotations belong to the core of the pangenome (This requires to have finished with the pangenome analysis first!)
	
	CoreDir=/mnt/283983cd-07b2-44de-b352-e9f8ccf9668e/Analysis/8.Anvio/2.Pangenome/1.Pangenome_Prokka_internal_genomes_ONLY/Taxonomy_groups_Functional_enrichment/CORE
	CoreAeruDir=/mnt/283983cd-07b2-44de-b352-e9f8ccf9668e/Analysis/8.Anvio/2.Pangenome/1.Pangenome_Prokka_internal_genomes_ONLY/Taxonomy_groups_Functional_enrichment/Core_aeruginosa
	CoreStutDir=/mnt/283983cd-07b2-44de-b352-e9f8ccf9668e/Analysis/8.Anvio/2.Pangenome/1.Pangenome_Prokka_internal_genomes_ONLY/Taxonomy_groups_Functional_enrichment/Core_stutzeri
	while IFS='@' read -r line
	do
	    printf "%s" "$line" > line
	    while read -r annot
            do
                grep -F "$annot" $CoreDir/* | sed -E 's/:/ /g' | awk '{print $1}' | head -n1 | awk -F"/" '{print $NF}' | sed -E 's/Functions_//g;s/\.tsv//g' > Core.txt
                if [[ -s Core.txt ]]; then
                    echo "yes" > In_Core.txt
                else
                    echo "no" > In_Core.txt
                fi
                
                grep -F "$annot" $CoreAeruDir/* | sed -E 's/:/ /g' | awk '{print $1}' | head -n1 | awk -F"/" '{print $NF}' | sed -E 's/Functions_//g;s/\.tsv//g' > CoreAeru.txt
                if [[ -s CoreAeru.txt ]]; then
                    echo "yes" > In_Core_Aeru.txt
                else
                    echo "no" > In_Core_Aeru.txt
                fi
                
                grep -F "$annot" $CoreStutDir/* | sed -E 's/:/ /g' | awk '{print $1}' | head -n1 | awk -F"/" '{print $NF}' | sed -E 's/Functions_//g;s/\.tsv//g' > CoreStut.txt
                if [[ -s CoreStut.txt ]]; then
                    echo "yes" > In_Core_Stut.txt
                else
                    echo "no" > In_Core_Stut.txt
                fi
                # Then paste together all newly created info with the original info we already had.
                paste -d "@" line In_Core.txt In_Core_Aeru.txt In_Core_Stut.txt | sed -e 's/\r//g;s/@/\t/g' >> Deleterious_temp.tsv
            done< <(awk -F'@' '{print $5}' line | sed 's/"//g') 
	done< <(sed -E 's/\t/@/g' Deleterious_AA_sub.tsv)
	
	awk 'BEGIN {FS=OFS="\t"; print "Reference_vs_Query_Strains", "Reference_Gene_Caller_ID","Source","Accession","Annotation","Reference_Amino_Acid","Query_Amino_Acid","SIFT_prediction","Annotation_in_core_pangenome","Annotation_in_aeruginosa_core_pangenome", "Annotation_in_stutzeri_core_pangenome"} NR>1 {print}' Deleterious_temp.tsv > Deleterious_final.tsv 
	
	rm line Core.txt In_Core.txt CoreAeru.txt In_Core_Aeru.txt CoreStut.txt In_Core_Stut.txt Deleterious_temp.tsv
	
	
	# If we want to find how many unique CDS we have in our results, type the following
	sed -E 's/(^Strain[[:digit:]][[:digit:]])_vs_Strain[[:digit:]][[:digit:]]/\1/g' Deleterious_final.tsv | awk '!seen[$1,$2]++' |sed -e "1d" |wc -l
	
	# If we want to fine the number of unique Strains we have in our results we do the following
	sed -E 's/(^Strain[[:digit:]][[:digit:]])_vs_Strain[[:digit:]][[:digit:]]/\1/g' Deleterious_final.tsv | awk '{print $1}' | sort -u -k1 | sed -e "1d" | wc -l
	
	# If we want to find the number of the unique CDS that belong to the core of the pangenome, we do the following:
	sed -E 's/(^Strain[[:digit:]][[:digit:]])_vs_Strain[[:digit:]][[:digit:]]/\1/g' Deleterious_final.tsv | awk '!seen[$1,$2]++' | awk '$NF=="yes"' | wc -l
		
HERE

echo
echo "==============================="
echo "SLURM job finished " "$(date)"
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
