#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Kraken"
#SBATCH --output=Kraken_job_%j.out

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

generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01
	This script runs from the master folder of Kraken and then we change directory to the folder of each specific Strain folder.

	NOTE:
	We are supposed to run this script after CLA, in order to filter out any scaffolds that is considered to be contamination.

	IMPORTANT:
	Please, read the report from the Kraken output and evaluate it. This is very important, because based on this report we try to filter out the scaffolds. Based on the report evaluation, make any adjustments to this script e.g. by changing the search pattern in the grep command.
	IMPORTANT:
	Here we take into account only the classified scaffolds by Kraken. Please check the Kraken output (StdOut) to check that ALL scaffolds have been classified. In case there are unclassified scaffolds, these should probably be included in the final output.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
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
StrNum=${StrainX/Strain/}

KrakenDir=/mnt/big/Metagenomics/kraken-2.0.7/installation
KrakenDB=/mnt/big/Metagenomics/kraken2-db
KrakenClassified=Kraken_Classified_CLA_scaffolds.fasta
KrakenOutput=Kraken_CLA_scaffolds.out
KrakenReport=Kraken_CLA_scaffolds.report

CLADir=/home1/pbravakos/Pseudomonas/CLA/${StrainX}/CLA-Results

#-----------------------------------------------------------------------------------------------------

# First, change directory
cd ${StrainX}

# Create a link of the CLA output file.
ln -s ${CLADir}/final_Scaffolds.fa

# Run Kraken
LC_ALL=C ${KrakenDir}/kraken2 --db ${KrakenDB} \
			--classified-out ${KrakenClassified} \
			--output ${KrakenOutput} \
			--use-names \
			--threads $SLURM_NTASKS \
			final_Scaffolds.fa

# Get the Kraken report. 
LC_ALL=C ${KrakenDir}/kraken2 --report ${KrakenReport} --threads $SLURM_NTASKS --db ${KrakenDB} $KrakenClassified

# Extract the fasta headers that have been classified as Pseudomonas by Kraken.
grep -P -o "[Sa-z_0-9]*\tPseudo" ${KrakenOutput} | grep -P -o "S[a-z_0-9]*" > Pseudomonas_fasta.headers

# Create the final fasta output. (We do this in order to make sure that it is empty before going to the while loop in the next step.)
cat /dev/null > ${StrX}_CLA_Kraken.fasta

# Extract (using awk) the Pseudomonas fasta from the CLA output into a new file. We also delete with sed all the empty lines (for a reason we get an empty line in the end of the file!)
while read line; do
    awk -v pattern="$line" 'BEGIN {RS=">"} $0 ~ pattern {print ">"$0}' final_Scaffolds.fa | sed -r '/^\s*$/d' >> ${StrX}_CLA_Kraken.fasta
 done < Pseudomonas_fasta.headers


echo "==============================="
echo "PBS job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
