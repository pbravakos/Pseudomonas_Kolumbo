#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="6th_Anvio_Prokka"
#SBATCH --output=6th_Anvio_Prokka_job_%j.out


# for calculating the amount of time the job takes
begin=$(date +%s)

# Some job specific info
echo "Job ID is = " $SLURM_JOBID
echo "SLURM cluster name = " $SLURM_CLUSTER_NAME
echo "SLURM partition = " $SLURM_JOB_PARTITION
echo "SLURM node list = " $SLURM_JOB_NODELIST
echo "SLURM num of nodes = " $SLURM_JOB_NUM_NODES
echo "SLURM number of tasks = " $SLURM_NTASKS
echo "SLURM memory per node = " $SLURM_MEM_PER_NODE
echo "SLURM memory per cpu = " $SLURM_MEM_PER_CPU
echo "working directory = " $SLURM_SUBMIT_DIR
echo "=================================================="
echo "SBATCΗ job started " `date`
echo "=================================================="
echo

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a

generalInfo () {
    cat <<END

	This is the Anvio pipeline. PART SIX.
	We will try to parse into anvio the BlastKOALA, Rast and Blast2Go annotations.
	Check the tutorial http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/ "Importing GhostKOALA/KEGG annotations into anvio".

	This script takes as input one argument. It runs from the master folder of Anvio and then we change directory to the directory of each specific Strain.
	
	For Strain 01 the correct usage is: 
	./$0 Strain01

	IMPORTANT!
	We have already run the Rast annotation online at http://rast.nmpdr.org/
	We have already run blastp on the anvio amino acid fasta files, exported the results to an xml file and then imported them into Blast2Go.
	We have also imported into Blast2Go the interproscan xml output files from the run we did earlier on this pipeline.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
    echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
    generalInfo >&2
    exit 1
fi

# Check that the working directory exists.
[[ ! -d $1 ]] && { echo "Directory $1 cannot be not found!" >&2; exit 1; }
echo
# Start a job dependency to move the Sdtout file to the correct folders
${HOME}/sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1
echo

# INITIAL PARAMETERS
StrainX=${1}
StrX=${StrainX/ain/}

#StrX=`echo ${StrainX} |cut -c1,2,3,7,8` # sh compatible command!

CntgDB=${StrX}_Pilon_contigs.db

BlastKOALAfile=${StrX}_user_ko.txt
Blast2GoDir=${HOME}/Pseudomonas/Blast2GO/${StrainX}
RastProkkaGeneCallsDir=${HOME}/Pseudomonas/Rast/Rast_with_Prokka_gene_calls/${StrainX}

Blast2GoTable=${StrX}_blast2go_go_table.txt



# Move to the working directory!
cd "${StrainX}"

#--------------------------SANITY CHECKS -------------------------------------------------------------------#
if [[ ! -e ${Blast2GoDir}/${Blast2GoTable} ]]; then
    echo "Blast2go output ${Blast2GoDir}/${Blast2GoTable} is missing!" >&2
    exit 1
fi

if [[ ! -e ${RastProkkaGeneCallsDir}/${StrX}_Rast.gff ]]; then
    echo "Rast gff file with Prokka gene calls cannot be found. Please upload it as ${RastProkkaGeneCallsDir}/${StrX}_Rast.gff" >&2
    exit 1
fi

if [[ ! -s ko00001.keg ]]; then 
    wget -N -O ko00001.keg "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext"
fi

# Check to see whether we have uploaded the BlastKOALA output file.
if [[ ! -r ${BlastKOALAfile} ]]; then
    echo "The file from BlastKOALa online analysis for $StrainX needs to be uploaded!" >&2
    exit 1
fi

#--------------------------SANITY CHECKS END! -------------------------------------------------------------------#


echo
echo '				Parse the results from BlastKOALA'
echo

# Create from the downloaded ko00001.keg file a tab delimited file 
# where the first column corresponds to the broadest classification, 
# the fifth corresponds to the gene itself

{
kegfile="ko00001.keg"
while read -r prefix content
do
    case "$prefix" in A) col1="$content";; 
                      B) col2="$content" ;; 
                      C) col3="$content";; 
                      D) echo -e "$col1\t$col2\t$col3\t$content";;
    esac 
done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > KO_Orthology_ko00001.txt
}

echo

# Next, we prepare the KEGG orthology file for comparison with the file we downloaded from BlastKOALA.
sed -E 's/K[0-9]{5};//g;s/(K[0-9]{5})\s*/\1\t/g' KO_Orthology_ko00001.txt | awk 'BEGIN{FS="\t";OFS="\t"} { print $4,$1"; "$2,$3,$5 }' > KEGG_orthology.tsv

# Next, we prepare the BlastKOALA file we downloaded from the internet, for comparison with KEGG orthology.
sed -E 's/^genecall_//g' "${BlastKOALAfile}" | awk 'BEGIN{FS="\t";OFS="\t"} {print $2,$1}' > "${StrX}"_BlastKoala.tsv


# Next, we compare the first column of each file, and when there is a match, we append the annotation of the first file to the second. Finally, we save to a file format that can be imported to Anvio.
awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$2;next}{print $2,$1,a[$1]}' KEGG_orthology.tsv "${StrX}"_BlastKoala.tsv | \
sed -E 's/\t([0-9]{5}) /\t\1;\t/g;s/; ([0-9]{5}) /; \t \1\t/g;s/Brite Hierarchies/-/g'  | \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1,"KEGG Brite/Pathway",$3$5,$4$6,"0"}' | \
awk -F"\t" '$3!=""' | awk '!a[$1$3]++' > AnvioKEGGBritePathwayImportable.tsv

# We repeat the same procedure for the Kegg Functions Annotations
awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$3;next}{print $2,$1,a[$1]}' KEGG_orthology.tsv "${StrX}"_BlastKoala.tsv | \
sed -E 's/\t([0-9]{5}) /\t\1\t/g' | \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1, "KEGG Function",$3,$4,"0"}' | \
awk -F"\t" '$3!=""' | awk '!a[$1$3]++' > AnvioKEGGFunctionImportable.tsv

# We repeat one last time the same procedure for Kegg Proteins Annotations.
awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$4;next}{print $2,$1,a[$1]}' KEGG_orthology.tsv "${StrX}"_BlastKoala.tsv | \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} {print $1, "KEGG Genes",$2,$3,"0"}' | \
awk -F"\t" '$3!=""' | awk '!a[$1$3]++' > AnvioKEGGGenesImportable.tsv

# Finally, we import these annotations to Anvio.
anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioKEGGBritePathwayImportable.tsv

anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioKEGGFunctionImportable.tsv

anvi-import-functions --contigs-db "${CntgDB}"  --input-files AnvioKEGGGenesImportable.tsv
echo
rm KEGG_orthology.tsv "${StrX}"_BlastKoala.tsv KO_Orthology_ko00001.txt
echo
echo '---------------------------BlastKOALA ANALYSIS COMPLETED---------------------------------------'	

echo
echo
echo
echo "			Parse Rast Annotations into Anvio"
echo


cp "${RastProkkaGeneCallsDir}"/"${StrX}"_Rast.gff .

sort -nk4,5 "${StrX}"_Rast.gff | \
grep -Ev "^#" | \
awk 'BEGIN{FS="\t";OFS="\t"} {print NR,$4,$5,$9}' | \
sed -E 's/ID=(fig\|[0-9.peg]+);/\1\t/g;s/Name=//g;s/Ontology_term=/ /g' > "${StrX}"_Rast_sorted.tsv



sort -nk4,5 "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff | \
grep -Ev "^#" | \
awk 'BEGIN{FS="\t";OFS="\t"} {print NR,$4,$5,$9}' > "${StrX}"_Prokka_sorted.tsv


awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$4;next} {print $4, $5,a[$1]}' "${StrX}"_Prokka_sorted.tsv "${StrX}"_Rast_sorted.tsv > "${StrX}"_Rast_Anvio_gene_caller_IDs.tsv


sed 's/ID=//g' "${StrX}"_Rast_Anvio_gene_caller_IDs.tsv | \
awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} $2!="hypothetical protein" {print $3,"Rast FigFam",$1,$2,"0"}' > AnvioRastAnnotationsImportable.tsv

anvi-import-functions --contigs-db "${CntgDB}"  --input-files AnvioRastAnnotationsImportable.tsv

rm "${StrX}"_Rast.gff "${StrX}"_Rast_sorted.tsv "${StrX}"_Prokka_sorted.tsv "${StrX}"_Rast_Anvio_gene_caller_IDs.tsv

echo
echo "-----------------------------------RAST ANNOTATIONS PARSED INTO ANVIO!!----------------------------------------------"
echo
echo
echo "			Parse Blast2Go Annotations into Anvio"
echo
## Start a job dependency for the next part of the Anvio pipeline to start as soon as this one has finished
#${HOME}/sbatch --dependency=afterok:$SLURM_JOB_ID 4.Anvio_Prokka_Jun18_PROFILING_AND_HMMs.sh ${StrainX}

cp "${Blast2GoDir}"/"${Blast2GoTable}" .

awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} NR>1 && $10!="" {print $3,"Blast2Go",$10,$11,$7}' "${Blast2GoTable}" > AnvioBlast2GoImportable.tsv

anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioBlast2GoImportable.tsv

awk 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","source", "accession","function","e_value"} NR>1 && $15!="no GO terms" && $15!="no IPS match" {print $3,"Interpo GoTerms",$15,$16,$7}' "${Blast2GoTable}" > AnvioInterProGoTermsImportable.tsv

anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioInterProGoTermsImportable.tsv

rm "${Blast2GoTable}"

echo
echo "----------------------------------------BLAST2GO ANALYSIS COMPLETED-------------------------------------------------------"
echo


'true'<< ////
	If we want to search a contigs database we can give a command like the following:
anvi-search-functions --contigs-db ${StrX}_Pilon_contigs.db --search-terms polyketide --output-file ${StrX}_polyketide_short_report.tsv --full-report ${StrX}_polyketide_full_report.tsv


	NOTE:
	The commands to run the blastp are the following: 
	blastp -query ${AnvioProkkaDir}/${StrX}_InterProScan_amino-acid-sequences.fasta \
        				-db nr \
        				-out ${StrX}_blastp_Anvio_Prokka.xml \
        				-task $Task \
        				-word_size $WordSize \
        				-evalue $Eval  \
        				-outfmt 14 \
        				-max_target_seqs $MaxSeqs \
        				-remote
////




echo
echo
echo "==============================================="
echo "SLURM job finished " "$(date)"
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
