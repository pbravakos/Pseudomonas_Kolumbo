#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="2nd_Anvio_Prokka"
#SBATCH --output=2nd_Anvio_Prokka_job_%j.out

# for calculating the amount of time the job takes
begin=`date +%s`

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


generalInfo () {
    cat <<END

	This script takes as input one argument. 
	For Strain 01 that would be: 
	$0 Strain01
	This script runs from the master folder of Anvio and then we change directory to the folder of each specific Strain folder.

	This is the Anvio pipeline. PART TWO

	IMPORTANT!
	This script uses the busco gammaproteobacteria single copy genes. It should be changed if we do not have gammaproteobacteria!!

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

# Start a job dependency to move the Sdtout file to the correct folders
${HOME}/sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
CntgDB=${StrX}_Pilon_contigs.db
InterproOutputDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}/${StrX}_INTEPRO

InterProDir=${HOME}/Software/interproscan-5.33-72.0
BuscoGammaProteoDir=${HOME}/Software/busco_hmms_for_anvio/gammaproteobacteria_anvio    # This is a Database specifically for gammaproteobacteria

KaijuDB=${HOME}/Software/kaiju/kaijudb

anviGetSeq='anvi-get-sequences-for-gene-calls'
anviImpTax='anvi-import-taxonomy-for-genes'

#------------------------------------------------------------------------------------------------------------
cd ${StrainX}


#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that Kaiju Database is installed.
if [[ ! -e ${KaijuDB}/kaiju_db.fmi ]] && [[ ! -e ${KaijuDB}/names.dmp ]] && [[ ! -e ${KaijuDB}/nodes.dmp ]]; then
    echo "Kaiju Database was not installed properly. Please check the database installation again." >&2
    exit 1
fi

# Check whether interproscan executable exists or not
if [[ ! -r ${InterProDir}/interproscan.sh ]]; then
    echo "Interproscan executable is missing!" >&2
    exit 1
fi

# Check that the Busco Directory can be found
if [[ ! -d $BuscoGammaProteoDir ]]; then
    echo "Busco Directory is missing! It shloud be found here: $BuscoGammaProteoDir" >&2
    exit 1
fi

## Start a job dependency for the next part of the Anvio pipeline to start as soon as this one has finished
#${HOME}/sbatch --dependency=afterok:$SLURM_JOB_ID 4.Anvio_Prokka_Jun18_PROFILING_AND_HMMs.sh ${StrainX}

echo
echo
echo '		Start Interproscan analysis'
echo
# Instructions can be found here http://merenlab.org/2016/06/18/importing-functions/

# For safety reasons copy the amino acid fasta file to a new one!
cp ${StrX}_amino_acid_seqs.faa ${StrX}_InterProScan_amino-acid-sequences.fasta

<< ////
	IMPORTANT!
	The "cpu" option just overrides both of the following values in the INTERPROSCAN properties file:
	#Number of embedded workers at start time
	number.of.embedded.workers=1
	#Maximum number of embedded workers
	maxnumber.of.embedded.workers=35
	https://github.com/ebi-pf-team/interproscan/issues/41
////


# Run interproscan for functional annotation
LC_ALL=C ${InterProDir}/interproscan.sh --input ${StrX}_InterProScan_amino-acid-sequences.fasta \
                                        --output-file-base ${StrX}_interpro-output \
					--cpu $SLURM_NTASKS \
                                        --formats TSV,GFF3,HTML,JSON,XML \
					--goterms \
					--iprlookup \
					--pathways \
					--disable-precalc \
					--applications CDD,Coils,Gene3D,Hamap,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,ProSitePatterns,ProSiteProfiles,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius

echo
# Move the Interpro output to a new directory, in case we delete the files in the future, to be able to save this directory at least.
mkdir ${InterproOutputDir}
mv ${StrX}_interpro-output* ${InterproOutputDir}
cp ${InterproOutputDir}/${StrX}_interpro-output.tsv .


# From the above list we have excluded SignalP and TMHMM, which we are going to run downstream in our analysis.

echo
echo "---------------------------INTERPROSCAN ANALYSIS COMPLETED!!-------------------------------------------"
echo
echo 
echo "			Import Pfam database search results into Anvio"
echo

# The following command will search the Pfam database and add a source named "Pfam" into Anvio
# anvi-run-pfams --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS

echo
echo
echo "-----------------------------------PFAM ANALYSIS COMPLETED!----------------------------------------------------"
echo
echo

echo '		Import NCBI COGs annotations into anvio'
# Instructions can be found here http://merenlab.org/2016/06/22/anvio-tutorial-v2/ and http://merenlab.org/2016/10/25/cog-annotation/
# NCBI’s COG database is quite outdated and not maintained but still awesome. The release is of December 2014! 
echo

# Prerequisite for this command to run is to have completed the anvi-setup-ncbi-cogs first!
anvi-run-ncbi-cogs --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --search-with blastp

echo
echo '----------------------------COG ANNOTATION COMPLETED---------------------------------------'
echo
echo
echo
echo "		Import Kaiju taxonomy into contigs Database"

# Instructions can be found here http://merenlab.org/2016/06/18/importing-taxonomy/
echo
# Fist get the nucleotide gene calls
${anviGetSeq} -c ${CntgDB} -o ${StrX}_gene_calls.fna
# File ${StrX}_gene_calls.fna is nucleotide sequences and for this reason Anvio exports every gene sequence
# so the good news is that we can use this file without any further analysis!!!
echo
echo '---------------------Nucleotide gene calls fasta file was created------------------------------'
echo
echo
echo "			Kaiju analysis starts now!!!"
echo
echo
## Instructions on how to run kaiju can be found here https://github.com/bioinformatics-centre/kaiju/blob/master/README.md
## IMPORTANT!!
## The following commands need the Refseq database for kaiju.
# We have downloaded the whole NCBI Refseq on 23/09/2018. This procedure takes a lot of time and space!!!

# Run the Kaiju classifier
LC_ALL=C kaiju -t ${KaijuDB}/nodes.dmp \
      			-f ${KaijuDB}/kaiju_db.fmi \
      			-i ${StrX}_gene_calls.fna \
      			-o ${StrX}_gene_calls_RefSeq.tab \
      			-z $SLURM_NTASKS \
      			-v

# Next, add taxon names 
LC_ALL=C addTaxonNames -t ${KaijuDB}/nodes.dmp \
              		-n ${KaijuDB}/names.dmp \
              		-i ${StrX}_gene_calls_RefSeq.tab \
              		-o ${StrX}_gene_calls_RefSeq.names \
              		-r superkingdom,phylum,order,class,family,genus,species

# At this point its not a bad idea to make a copy of the contigs database –just in case
cp ${CntgDB} ${CntgDB}.bak

# Finally, run the anvio parser for Kaiju
${anviImpTax} --input-files ${StrX}_gene_calls_RefSeq.names \
                       --contigs-db ${CntgDB} \
                       --parser kaiju \
                       --just-do-it

echo
echo '------------------------KAIJU TAXONOMOMY ANALYSIS COMPLETED------------------------------------'
echo
echo
echo
echo '		Parse HMM hits to the contigs database'
echo
# Check this tutorial:
# http://merenlab.org/2015/06/25/screening-cultivars/ "Removing contaminants from cultivars with anvi'o"

# We can use the program busco_hmms_for_anvio (https://github.com/guyleonard/busco_hmms_for_anvio) to download 
# (we did it on 29/09/2018) and extract the gammaproteobacteria hmm busco dataset of single copy genes found here: 
# http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz
# This dataset can be used as input to the anvi-run-hmms command
anvi-run-hmms --num-threads $SLURM_NTASKS -c ${CntgDB} -H ${BuscoGammaProteoDir}
echo
echo
echo
# The following command will use the default anvio bacterial single copy gene database.
anvi-run-hmms --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS

echo
echo '-------------------------------HMM HITS ANALYSIS COMPLETE-------------------------------------'
echo
echo
echo
echo '		How many genomes?'
echo

anvi-script-gen_stats_for_single_copy_genes.py ${CntgDB}
anvi-script-gen_stats_for_single_copy_genes.R ${CntgDB}.hits ${CntgDB}.genes --output_prefix=${StrX}_single_copy_genes
# The result of the above command is a pdf for later viewing!
echo
echo "--------------------------SINGLE COPY GENES ANALYSIS IS COMPLETE---------------------------------"
echo
echo
echo
echo '		Contigs Statistics'
echo

# Shows simple stats of the contigs database to assess your assembly output, and also estimate the number of bacterial genomes to recover
anvi-display-contigs-stats --report-as-text --output-file ${StrX}_contig_stats.tab ${CntgDB}

echo
echo '--------------------------------CONTIGS STATISTICS ANALYSIS COMPLETE----------------------------------'
echo
echo
echo
echo "			Export PANTHER annotations for online analysis"
echo
<< ////
	NOTE
	The sequence identifiers (gene calls) have been mapped to PANTHER HMM IDs by InterproScan and can be used in the gene list analysis tools.
	Panther Generic Mapping File
	For IDs from organisms other than the 82 organisms in the PANTHER database, user-generated data containing mappings between those IDs and their corresponding PANTHER IDs can be used. The file must be tab-delimited and must contain the following columns: the first column can contain a list of unique IDs from the user; the second column should be the corresponding PANTHER family or subfamily ID (e.g., PTHR10078 or PTHR10078:SF6), and is used to look up the association with GO and PANTHER terms (molecular function, biological process and pathway).
	from: https://www.nature.com/articles/nprot.2013.092

////


# We export the Panther annotation from Interproscan results and we also make the Panther IDs unique.
awk -F'[\t]' '$4 ~ /PANTHER/ {print $1"\t"$5}' ${StrX}_interpro-output.tsv | sort -k2 -u > ${StrX}_Panther.tab


# FURTHER STEPS NEEDED!


echo
echo "----------------------------------------------PANTHER FILE READY FOR UPLOAD!-------------------------------------------------------" 
echo
echo
echo " 			Export PIRSF annotations for online analysis"
echo

# We first export the PIRSF annotation from Interproscan results and we also make the PIRSF IDs unique.
awk '$4 ~ /PIRSF/{print $5}' ${StrX}_interpro-output.tsv | sort -u > ${StrX}_PIRSF_InterproScan_IDS.txt

echo
echo "----------------------------------------------PIRSF FILE READY FOR UPLOAD!-------------------------------------------------------" 
echo
<< ////
	FURTHER STEPS NEEDED!

	FOR PANTHER
	Now it is time to go online to http://www.pantherdb.org/ and upload the file "${StrX}_Panther.tab" as Panther Generic Mapping File. Please, do not forget to change the Reference List from Human to Pseudomonas aeruginosa.

	FOR PIRSF
	Now it is time to go online to 	https://pir.georgetown.edu/pirwww/search/batch_sf.shtml and copy and paste the IDs found in the file "${StrX}_PIRSF_InterproScan_IDS.txt". Then, click on the "Save Results As: Table", download the output and rename the file to "${StrX}_PIRSF.tsv". Upload "${StrX}_PIRSF.tsv" to the correct folder on the server, and we are ready to go!
////
echo
echo
# Remove files no longer needed!
# rm ${StrX}_InterProScan_columns_12456_filtered.tsv ${StrX}_InterProScanVariousDatabasesImportable.tsv ${StrX}_InterProScan_columns_1213_filtered.tsv AnvioInterProAccessionsImportable.tsv ${StrX}_InterProScan_columns_14_filtered.tsv AnvioGoTermsImportable.tsv ${StrX}_InterProScan_columns_15_filtered.tsv AnvioInterProPathwaysImportable.tsv ko00001.keg KO_Orthology_ko00001.txt KeggOrthology_Table1.txt ${StrX}_KeggAnnotations-AnviImportable.txt


echo "==============================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
