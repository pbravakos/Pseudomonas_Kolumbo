#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="5th_Anvio_Prokka"
#SBATCH --output=5th_Anvio_Prokka_job_%j.out

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

	This is the 5th part of the Anvio Prokka pipeline. 
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01

	It runs from the master folder of Anvio and then we change directory to the folder of each specific Strain folder.

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
StrainX=$1
StrX=${StrainX/ain/}
CntgDB=${StrX}_Pilon_contigs.db

AntiSmashDir=${HOME}/Software/antismash-4.2.0
PfamDir=${HOME}/Software/Pfam32.0
OutputDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}/Antismash
AnvioDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}
Genbank2GffDir=${HOME}/Software/gfftools
EmapperDir=${HOME}/Software/eggnog-mapper

ContigPilon=${StrX}_Pilon_CLA_Blast.fasta

GProteoNOGDB=${HOME}/Software/gproNOG_hmm # This is the Proteobacteria Database for Eggnog mapper

EggNOGParseKEGGDir=${HOME}/Software/EggNogKEGGParser

AnvioGffONLYProkka=${StrX}_ONLY_PROKKA_amino_acid_seqs.gff

# Emapper parameters that can be changed
TaxScope=gproNOG # This is to limit the search on a specific database e.g. gamma Proteobacteria.
MaxSeqLen=5000
Eval=0.001


cd ${StrainX}

#--------------------------SANITY CHECKS -------------------------------------------------------------------#

# Check that the Contigs file already exists in the working folder.
if [[ ! -e $ContigPilon ]]; then
    echo "The contigs fasta file $ContigPilon is missng. " >&2
    exit 1
fi

if [[ ! -s ko00001.keg ]]; then 
    wget -N -O ko00001.keg "https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext"
fi

# Check that Antismash run file exists.
if [[ ! -e ${AntiSmashDir}/run_antismash.py && ! -e ${Genbank2GffDir}/genbank2gff.pl && ! -e ${EmapperDir}/emapper.py && ! -e ${EggNOGParseKEGGDir}/KEGG-to-anvio ]]; then 
    echo "Some or all of the following tools have not been installed properly: Antismash, gfftools, eggonog mapper, kegg-to-anvio (on ${EggNOGParseKEGGDir})." >&2
    exit 1
fi

# Check that the gff annotation file with only Prokka genes exists
if [[ ! -s $AnvioGffONLYProkka ]]; then
    echo "The Anvio gff file containing only prokka gene calls is missing!" >&2
    exit 1
fi

# Check that the Proteobacteria database for Eggnog mapper has been installed
if [[ ! -e $GProteoNOGDB ]]; then
   echo "The Proteobacteria database for Eggnog mapper has not been downloaded" >&2
   exit 1
fi

# Check that gfftools have been installed.
if [[ ! -e ${Genbank2GffDir}/genbank2gff.pl ]]; then
    echo "Gfftools have not been installed. Check the folder ${Genbank2GffDir}" >&2
    exit 1
fi

# Check that the Pfam hmm file exists.
if [[ ! -e $PfamDir/Pfam-A.hmm ]]; then
   echo "Pfam required files cannot be found. Please check the directory $PfamDir." >&2
   exit 1
fi



## Start a job dependency for the next (and final) part of the Anvio pipeline to start as soon as this one has finished
#${HOME}/sbatch --dependency=afterok:$SLURM_JOB_ID 4.Anvio_Prokka_Jun18_PROFILING_AND_HMMs.sh ${StrainX}


echo
echo "		Antismash Analysis begins"
echo

if [[ ! -d ${OutputDir} ]]; then
mkdir ${OutputDir}
fi

## Create a soft link of the contigs fasta file (This file should already exist, from previous parts of the Anvio pipeline!)
#ln -s ${HOME}/Pseudomonas/Seqtk/${StrainX}/${ContigPilon}

# IMPORTANT
# We are going to use the gff file which has the gene caller ids, only from Prokka CDS! We have already created this file in the first part of the pipeline.


echo
echo " 				Start Antismash Analysis!!!"
echo

${AntiSmashDir}/run_antismash.py --cpus $SLURM_NTASKS \
				 --verbose \
				 --debug \
				 --taxon bacteria \
				 --input-type nucl \
				 --transatpks_da \
				 --clusterblast \
				 --subclusterblast \
				 --knownclusterblast \
				 --smcogs \
				 --inclusive \
				 --borderpredict \
				 --full-hmmer \
				 --asf \
				 --tta \
				 --minlength 200 \
				 --nclusters 10 \
				 --cf_cdsnr 5 \
				 --cf_npfams 5 \
				 --fix-id-line \
				 --pfamdir ${PfamDir} \
				 --gff3 ${AnvioGffONLYProkka} \
				 --outputfolder ${OutputDir} \
				 --enable-BiosynML ${ContigPilon}
echo
echo
echo "---------------------------------ANTISMASH ANALYSIS COMPLETED!-------------------------------------------------"
echo
echo
echo
echo "				Import Antismash SmCOG Analysis into Anvio"
echo


cd ${OutputDir}

<< ////
	Secondary metabolism gene families (smCOGs) analysis attempts to allocate each gene in the detected gene clusters to a secondary 
	metabolism-specific gene family using profile hidden Markov models specific for the conserved sequence region characteristic of this family.
////

LC_ALL=C ${Genbank2GffDir}/genbank2gff.pl *...final.gbk > ${StrX}_Antismash.gff

sed -E -i "s/[\']//g;s/\(Score//g;s/\)\;/;/g;s/gene=/=/g;s/NADH\://g" ${StrX}_Antismash.gff

printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioSMCOGsImportable.tsv

awk -F'[\t=:;]' '$3 ~ /CDS/ && $12 ~ /smCOG/ {print $10"\t"$12"\t"$13"\t"$14"\t"0}' ${StrX}_Antismash.gff >> AnvioSMCOGsImportable.tsv

mv AnvioSMCOGsImportable.tsv ${AnvioDir}

cd ${AnvioDir}

anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioSMCOGsImportable.tsv

echo
echo "----------------------------------ANTISMASH SmCOG RESULTS PARSED INTO ANVIO!!--------------------------------------------"
echo
echo
echo
echo "		Eggnog Mapper Analysis begins!!"
echo


## hmmer search
#${EmapperDir}/emapper.py --database gproNOG \
#			 --hmm_maxseqlen ${MaxSeqLen} \
#			 --hmm_evalue ${Eval} \
#			 --override \
#			 -i ${StrX}_InterProScan_amino-acid-sequences.fasta \
#			 --output ${StrX}_hmmer \
#			 --temp_dir temp \
#			 --cpu $SLURM_NTASKS 
			

# diamond search
${EmapperDir}/emapper.py --database gproNOG \
			--override \
			-i ${StrX}_InterProScan_amino-acid-sequences.fasta \
			--output ${StrX}_diamond \
			--temp_dir temp \
			--cpu $SLURM_NTASKS \
			-m diamond \
			--matrix BLOSUM62 \
			--tax_scope $TaxScope

#emapper-1.0.3


echo
echo "------------------------------EGGNOG MAPPER ANALYSIS COMPLETED!---------------------------------"
echo
echo "			Import Eggnog Annotations into Anvio"
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


printf '%s\t%s\t%s\t%s\t%s\n' 'gene_callers_id' 'source' 'accession' 'function' 'e_value' > AnvioEGGNOGImportable.tsv


grep -v '^#' ${StrX}_diamond.emapper.annotations | awk  -F"\t" '{print $1"\t""EGGNOG""\t"$2"\t"$13"\t""0"}' | sed -E 's/\b[0-9]+\.//g' >> AnvioEGGNOGImportable.tsv 

anvi-import-functions --contigs-db ${CntgDB} --input-files AnvioEGGNOGImportable.tsv 



grep -v '^#' ${StrX}_diamond.emapper.annotations | awk -F"\t" '{print "genecall_"$1"\t"$7}' > ${StrX}_EggNOGKEGGParserInput.txt

# Fist format the KEGG results into a format that can be parsed into anvio
${EggNOGParseKEGGDir}/KEGG-to-anvio --KeggDB KO_Orthology_ko00001.txt \
					-i ${StrX}_EggNOGKEGGParserInput.txt \
					-o ${StrX}_EggNOGKEGGAnnotations-AnviImportable.tsv


# Second import the results to anvio
anvi-import-functions --contigs-db ${CntgDB} --input-files ${StrX}_EggNOGKEGGAnnotations-AnviImportable.tsv

rm KO_Orthology_ko00001.txt AnvioEGGNOGImportable.tsv ${StrX}_EggNOGKEGGParserInput.txt

echo
echo "----------------------------------EGGNOG ANNOTATIONS PARSED INTO ANVIO!------------------------------------------------------------"
echo
# Depracated!
#anvi-script-run-eggnog-mapper --contigs-db ${CntgDB} \
#				--annotation ${StrX}_diamond.emapper.annotations \
#				--use-version 1.0.3\
#				--num-threads $SLURM_NTASKS

<< ////
	For a manual download of the eggonog database please do the following:
	1) Download all HMM models for your target taxonomic level. 
	For instance, if Bacteria, download them from http://eggnogdb.embl.de/download/eggnog_4.5/data/bactNOG/
	2) Build a HMMER database using hmmpress. 
	For this, you will need to concatenate all models in a single file (i.e. cat bactNOG_hmm/*.hmm > bactDB.hmmer), 
	and run hmmpress bactDB.hmmer.
	3) Use hmmscan to query your sequences against the database.For instance: hmmscan bactDB.hmmer MyQueryFasta.fa
	
	More specifically:
	The following commands should run only once!
	Download the gammaproteobacteria database, concatenate the hmm files and run hmmpress!	
	wget http://eggnogdb.embl.de/download/eggnog_4.5/data/gproNOG/gproNOG.hmm.tar.gz
	tar xvf gproNOG.hmm.tar.gz
	cd gproNOG_hmm
	cat *.hmm > gproNOG_DB.hmmer
	hmmpress gproNOG_DB.hmmer
	hmmscan -o eggnog_${StrainX}_hmmscan.output --tblout eggnog_${StrainX}_hmmscan_per_sequence.table \
					--domtblout eggnog_${StrainX}_hmmscan_per_domain.table \
					--pfamtblout eggnog_${StrainX}_hmmscan_pfam-like.table --noali -E 0.0001 \
					--domE 0.0001 --cpu $SLURM_NTASKS \
					${GProteoNOGDB}/gproNOG_DB.hmmer \
					${StrX}_InterProScan_amino-acid-sequences.fasta

////

echo
echo "======================================================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
