#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="1st_Anvio_Prokka"
#SBATCH --output=1st_Anvio_Prokka_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

# Some job specific info
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

set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a

generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	./$0 Strain01
	This script runs from the master folder of Anvio and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	This is the Anvio pipeline. PART ONE
	Check the tutorial http://merenlab.org/2016/06/18/importing-functions/ "Importing functions into contigs database"

	For this script to run we should already have finished the Prokka annotation, because we use the output of the Prokka pipeline as input to the contigs database.
	We should also have done the genecalling with GeneMarkS2, online at http://exon.gatech.edu/GeneMark/genemarks2.cgi.

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
sbatch --dependency=afterany:"$SLURM_JOB_ID" Move_Slurm_output_files.sh "$SLURM_JOB_NAME" "$SLURM_JOB_ID" "$1"
echo
# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
CntgDB=${StrX}_Pilon_contigs.db

ContigPilonDir=${HOME}/Pseudomonas/Seqtk/${StrainX}
ProkkaDir=${HOME}/Pseudomonas/Prokka/${StrainX}/${StrX}_Prokka
GeneMarkS2Dir=${HOME}/Pseudomonas/GeneMarkS2/${StrainX}
GfftoGenbankDir=${HOME}/Software/Gff_to_genbank

ContigPilon=${StrX}_Pilon_CLA_Blast.fasta
ProkkaGff=${StrX}_Prokka.gff
GeneMarkGff=${StrX}_gms2.gff

# Variables for Anvio command options that can be altered
ContigDBSplitLength=20000 # This is mostly for better visualization. Not to be changed.
ContigDBKmerSize=4 # Tetranucleotides is the default and the most used one. Not to be changed.

# We substitute the command anvi-get-sequences-for-gene-calls with a parameter, but this is certainly not necessary!
anviGetSeq='anvi-get-sequences-for-gene-calls'

echo
## Start a job dependency for the next part of the Anvio pipeline to start as soon as this one has finished
#${HOME}/sbatch --dependency=afterok:$SLURM_JOB_ID 2.Anvio_Prokka_Interpro_Pfam_KEGG.sh $StrainX
echo

#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that the files we have access to the files needed to run this pipeline
if [[ ! -e $ContigPilonDir/$ContigPilon ]] || [[ ! -e $ProkkaDir/${ProkkaGff} ]] || [[ ! -e ${GeneMarkS2Dir}/${GeneMarkGff} ]]; then
    echo "Some or all of the files necessary to run this pipeline (Contigs, Prokka annotation, GeneMarkS2 annotation) are missing!" >&2
    exit 1
fi

if [[ ! -e ${GfftoGenbankDir}/gff_to_genbank.py ]]; then
    echo "Software to get genbank sequences is missing in ${GfftoGenbankDir}/gff_to_genbank.py" >&2
    exit 1
fi

echo
anvi-self-test --version
echo
echo
echo "			Anvio Analysis begins!!"
echo

cd "${StrainX}" || { echo "Directory ${StrainX} cannot be not found!" >&2; exit 1; }

echo
# Create a soft link of the contigs fasta file in the working directory
ln -s "${ContigPilonDir}"/"${ContigPilon}"
echo
echo
echo
echo '		Export Prokka gene calls to an importable format for Anvio '
echo

# Get the gene calls and annotations from a gff3 Prokka output in a format that anvio can import.

# We have to subtract one (1) from the start position column. Please check here http://merenlab.org/2016/06/22/anvio-tutorial-v2/#external-gene-calls

# First, export CDS gene calls found by prodigal, and save it to a temporary file. We will add the gene caller ids at the end, after we concatenate all the temporary files with the gene calls into one file.
awk '$3 ~ /CDS/{print $1,($4-1),$5,$7,$8,"prodigal","2.6.3"}' "${ProkkaDir}"/"${ProkkaGff}" > "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Second, we export Aragorn gene calls as partial genes (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS annotation.
awk '$2 ~ /Aragorn:1.2/{print $1,($4-1),$5,$7,"1","Aragorn","1.2"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Third, we export barnap gene calls as partial genes (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS and Aragorn annotation.
awk '$2 ~ /barrnap:0.9/{print $1,($4-1),$5,$7,"1","Barnap","0.9"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# Fourth, we export Infernal gene calls as partial genes (to prevent Anvio from exporting these genes as amino acids!) and append it to the same temporary file we save the CDS, Aragorn and barnap annotation.
awk '$2 ~ /Infernal:1.1/{print $1,($4-1),$5,$7,"1","Infernal","1.1"}' "${ProkkaDir}"/"${ProkkaGff}" >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

# NOTE
# We do NOT export SignalP from the Prokka gff file, as we will search and import these later on

#------------------------------------------------------------------------------------------------------------------------------------

# Next we will try to export gene calls from GenemarkS2, which we have already run online http://exon.gatech.edu/GeneMark/genemarks2.cgi
# We will export only gene calls that have not been found by Prokka.
# We MUST append the GeneMarkS2 gene calls last, after all the other gene calls have already been added to our temporary file!

# Prepare the two files (one from Prokka, one from GenemarkS2) for comparison
awk '$3 ~ /CDS/  {print $1,$4,$5,$7,$8}' "${ProkkaDir}"/"${ProkkaGff}" > "${StrX}"_Prokka_gff.out
awk '$3 ~ /CDS/  {print $1,$4,$5,$7,$8}' "${GeneMarkS2Dir}"/"${GeneMarkGff}" > "${StrX}"_GeneMark_gff.out

# Check if GeneMarkS2 gene calls can be found in Prokka gene calls. If not, then write these gene calls to a new file.
{
while read -r p; do
  if ! grep -q "$p" "${StrX}"_Prokka_gff.out; then
    echo "$p" >> GeneMarkS2_genes_NOT_in_Prokka.tsv
  fi
done < "${StrX}"_GeneMark_gff.out
}

# Prepare the temporary file (i.e. every column exept for the gene call ids). Please don't forget that we have to subtract 1 (one) from the start position!
awk 'OFS="\t" {print $1,($2-1),$3,$4,$5, "GeneMarkS2","3"}' GeneMarkS2_genes_NOT_in_Prokka.tsv >> "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv

#------------------------------------------------------------------------------------------------------------------------------------

# Create the necessary headers for the output file
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'contig' 'start' 'stop' 'direction' 'partial' 'source' 'version' > "${StrX}"_Prokka_Anvio_gene_calls.tsv

# Finally, we append all the gene calls from the temporary file we created earlier to the final file we are going to use to parse all the Prokka annotations into Anvio. 
# IMPORTANT
# We have to subsitute all the +,- with f,r and print it with the first column being the number of the row (which corresponds to the gene caller id!)
awk 'OFS="\t" {sub("+","f",$4);sub("-","r",$4)} {print NR,$0}' "${StrX}"_Prokka_Anvio_temporary_gene_calls.tsv >> "${StrX}"_Prokka_Anvio_gene_calls.tsv

echo "----------------------PROKKA GENE CALLS WERE EXPORTED!!!!---------------------------------------------------------"
echo
echo
echo '			Create Contigs DATABASE'
echo
# Create the Contig Database description file
echo "${StrX} contigs database created from the Pilon output scaffold fasta files, after filtering out contigs marked either as contamination or plasmids by blast search. The exact pipeline that created these scaffolds is described elswhere. Scaffolds were created only from one Illumina Miseq run performed on June 2018! Gene finding was done inside Prokka  v1.13.3 with Prodigal 2.6.3. Annotation was Completed with Prokka, using the genus Pseudomonas command in Prokka (after creating the appropriate database). We have also included gene calls from GeneMarkS2, but only those that have not already been found by Prokka (prodigal)" >  "${StrX}"_ContigDB_Anvio_decription.txt

# Finally, create the Contigs Database with the Prokka gene calls.
anvi-gen-contigs-database -f "${ContigPilon}" \
			  --project-name "${StrX}"_Prokka_Anvio \
			  --output-db-path "${CntgDB}" \
                          --description "${StrX}"_ContigDB_Anvio_decription.txt \
			  --split-length "${ContigDBSplitLength}"\
                          --kmer-size "${ContigDBKmerSize}" \
			  --external-gene-calls "${StrX}"_Prokka_Anvio_gene_calls.tsv \
			  --ignore-internal-stop-codons # We ignore stop codons in order to avoid to be forced to a stop by Anvio. But we should not forget to check the stdout file!
echo
echo
echo '---------------------------------CONTIGS DATABASE CREATED!!!---------------------------------------'
echo
echo
echo '		Parse Prokka Functional Annotations into Anvio'
echo

'true' << ////
	Now, we will export functional annotations from the ${ProkkaGff} file output of the Prokka pipeline to an importable tsv file
	
	IMPORTANT no1
	Gene caller ids should be the same in this file as the ones assigned by Anvio to the contigs database. 
	For this reason we use a similar awk command (i.e. the same .gff file, filtered by CDS/Aragorn/Barnap/Infernal) as the one used 
	to parse the gene calls into Anvio in the previous steps.
	
	IMPORTANT no2
	We MUST follow the same order of parsing these annotations as was the order that we parsed the gene calls.
	i.e. First MUST be parsed CDS, second MUST be parsed Aragorn, third MUST be parsed Barnap and last MUST be parsed Infernal!!
////

# Fisrt, export  CDS annotations to an importable file format, and save these annotations to a temporary file
awk -F'[\t=;,]' '$3 ~ /CDS/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
# We also remove non useful info (with sed) to have finally only the gene accession and function information
sed -E 's/ID=.+sequence://g; s/ID=.+motif://g; s/;locus_tag=.+;product=/\t/g' | \
awk -F"\t" '{print "Prokka""\t"$9"\t"$10}' > AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# At this point, we want to have the number of CDS annotations produced by Prokka, and save it to a variable for later use.
ProkkaCDSNum=$(wc -l AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv | cut -d' ' -f1)

# Second, we follow the same procedure to export the Aragorn annotation, and append to the same temporary file we saved the CDS annotations.
awk -F'[\t=;]' '$2 ~ /Aragorn:1.2/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Aragorn""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# Third, we repeat the same procedure to export the Barnap annotation, and append to the same temporary file we saved the CDS and Aragorn annotations.
awk -F'[\t=;]' '$2 ~ /barrnap:0.9/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Barnap""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

# Fourth, we repeat the same procedure to export the Infernal annotation, and append to the same temporary file we saved the CDS, Aragorn and Barnap annotations.
awk -F'[\t=;]' '$2 ~ /Infernal:1.1/ {print $0}' "${ProkkaDir}"/"${ProkkaGff}" | \
sed -E 's/ID=.+locus_tag=//g; s/;product=/\t/g' | \
awk -F"\t" '{print "Prokka:Infernal""\t"$9"\t"$10}' >> AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv

#------------------------------------------------------------------------------------------------------------------------------------------

# Create the file to be parsed into Anvio with the correct headers, in order to append later the Prokka functional annotation
printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioProkkaFunctionalAnnotationsImportable.tsv

# Finally parse the temporary file we created to Anvio, with the correct gene calls, and (after that!!) apply any further filtration, like removing hypothetical proteins etc. Finally we fix the e-value to 0 (zero).
awk '{print NR"\t"$0}' AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv | \
awk -F"\t" '$4!="hypothetical protein" {print $1"\t"$2"\t"$3"\t"$4"\t""0"}' >> AnvioProkkaFunctionalAnnotationsImportable.tsv


# Finally, parse functional annotations into Anvio!
anvi-import-functions --contigs-db "${CntgDB}" --input-files AnvioProkkaFunctionalAnnotationsImportable.tsv

echo
echo '-------------------------PROKKA FUNCTIONAL ANNOTATIONS WERE PARSED INTO ANVIO!-------------------------------'
echo
echo
echo '				Export Amino Acid fasta file with NO empty records. Also export genes to  a gff file'
echo
echo

# We are going to use the amino acid fasta file from the Anvio Prokka pipeline as input for all subsequent analyses
# For this reason, we want to have an amino acid sequence fasta file without any empty fasta records. These empty records could potentially cause problems in the analysis downstream.

# Export amino acids from the the contigs database into a fasta file
$anviGetSeq --contigs-db "${CntgDB}" --get-aa-sequences --output-file "${StrX}"_amino_acid_seqs_with_empty_fields.fasta

# Remove empty fasta records with awk.
# With RS, we declare that the records are separated by >. If the record contains a second line $2 (i.e. is not an empty fasta record), print the record (add a > in front).
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' "${StrX}"_amino_acid_seqs_with_empty_fields.fasta > "${StrX}"_amino_acid_seqs.faa
echo

# Export gff3 from anvi contig database
${anviGetSeq} -c "${StrX}"_Pilon_contigs.db -o "${StrX}"_Anvio.gff --export-gff3

echo
echo
echo "--------------------------------------AMINO ACID FASTA FILE READY-------------------------------------------------------------"
echo
echo "		Create new fasta file WITHOUT the GenemarkS2 gene calls, ONLY PROKKA CDS gene calls!. Also create a gff file with ONLY PROKKA CDS annotations"
echo


# Next, we create the file "patterns.txt" which will be used as input in grep. Patterns.txt contains one gene caller id per line, for all genes found by Prokka.
for i in $(seq 1 "$ProkkaCDSNum"); do echo "^>$i\b" >> patterns.txt; done

# Finally, we linearize the FASTA with awk, pipe to grep to filter for items of interest named in patterns.txt, then pipe to tr to delinearize:
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "${StrX}"_amino_acid_seqs.faa | \
grep -Ef patterns.txt - | \
tr "\t" "\n" >  "${StrX}"_ONLY_PROKKA_amino_acid_seqs.faa
echo

# Also, create a gff file with ONLY Prokka CDS gene calls.
awk -v var="$ProkkaCDSNum" -F"[\t=]" '(NR==1) || ($10 <= var )' "${StrX}"_Anvio.gff > "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff
# IMPORTANT!
# Only this gff file can actually be used because otherwise we get all the RNAs and non coding RNAs as CDS, which is a huge mistake!!! 

echo
echo
echo '--------------------------AMINO ACID FASTA SEQUENCES and GFF file ONLY FROM PROKKA CDS ANNOTATIONS ARE READY!----------------------------------'
echo
echo 
echo "				Create Genbank File only from Prokka gene calls"
echo

# We want to use this file for uploading to the Rast servers, and later parse the Rast annotations into Anvio, but besides this it is not a bad idea to have a genbank file.

sed 's/Scaffold_//g; s/_length//g; s/_pilon//g' "${ContigPilon}" > "${StrX}"_contigs_for_gff_to_genbank.fasta
sed 's/Scaffold_//g; s/_length//g; s/_pilon//g' "${StrX}"_ONLY_PROKKA_amino_acid_seqs.gff > "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff

"${GfftoGenbankDir}"/gff_to_genbank.py "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff "${StrX}"_contigs_for_gff_to_genbank.fasta

mv "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gb "${StrX}"_ONLY_PROKKA_aa.gbk

rm "${StrX}"_contigs_for_gff_to_genbank.fasta "${StrX}"_ONLY_PROKKA_aa_for_gff_to_genbank.gff


# Now we have to go to the Rast server http://rast.nmpdr.org/ and upload the ${StrX}_ONLY_PROKKA_aa.gbk. When we have the results we have to upload them to the working directory in order to import these annotations into Anvio!



echo
echo
echo '		Export amino acid sequences from Anvio for KEGG functional annotation'
echo
# Instructions can be found here http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
echo

# For safety reasons copy the amino acid fasta file to a new one!
cp "${StrX}"_amino_acid_seqs.faa "${StrX}"_KEGG_protein-sequences.fasta


# BlastKOALA does have a gene number limit of 5000 genes. 
# If your anvi’o database has more than that number of genes, you can split the file into multiple FASTA files, 
# and concatenate the BlastKOALA outputs you will be generating during this step.

# The fasta sequence headers start with a number (e.g. 0 or 1 etc) and the next command adds 'genecall' as a prefix of each header
sed -i 's/>/>genecall_/g' "${StrX}"_KEGG_protein-sequences.fasta

echo
echo '--------------------------------BlastKOALA/KEEG AA FASTA FILE READY FOR UPLOAD!-----------------------------------------------------'


# Remove files not needed downstream in our analysis.
rm ./*glob*.out ./*glob*temporary* patterns.txt AnvioProkkaTEMPFunctionalAnnotationsImportable.tsv "${StrX}"_amino_acid_seqs_with_empty_fields.fasta

# ATTENTION, FURTHER STEPS NEEDED!
# ${StrX}_KEGG_protein-sequences.fasta file is now ready to be submitted to BlastKOALA for annotation.
# To do this go to the BlastKOALA webserver https://www.kegg.jp/blastkoala/, 
# and click the “Choose File” button underneath the section that says Upload query amino acid sequences in FASTA format. 
# From the menu you will upload the ${StrX}_KEGG_protein-sequences.fasta, and will be asked to provide an email address.

# Also don't forget to run the Rast annotation online at http://rast.nmpdr.org/, selecting not Rasttk, but Rast while keeping the gene calls from the genbank file!

	
echo
echo "==============================="
echo "SLURM job finished " "$(date)"
echo

# finished commands

# getting end time to calculate time elapsed
end=$(date +%s)
elapsed=$(${$end - $begin})
echo Time taken: "$(printf '%dh:%dm:%ds\n' "${elapsed}"/3600 "${elapsed}"%3600/60 "$elapsed"%60)"

exit 0
