#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="KmerGenie"
#SBATCH --output=KmerGenie_job_%j.out

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
	For Strain01 that would be: 
	$0 Strain01
	This script runs from the master folder of KmerGenie and then we change directory to the folder of each specific Strain folder.

	KmerGenie estimates the best k-mer length for genome de novo assembly. http://kmergenie.bx.psu.edu/
	
	NOTE:	
	We are working here with KmerGenie version 1.7016 because later versions do not run on the server.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ Strain[0-9]{2} ]]; then
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
KmergenieDir=${HOME}/Software/kmergenie-1.7016
DateRun=Nov18b

ReadsFileTemplate=Reads.template
ReadsFileStrainX=${StrainX}_${DateRun}_reads.txt

if [[ ! -e ${ReadsFileTemplate} ]]; then
    echo "$ReadsFileTemplate cannot be found in ${SLURM_SUBMIT_DIR}. Please upload it." >&2
    exit 1
fi

#==================================================================================================

cd ${StrainX}

# Create a new reads file with bash parameter substitution.
( echo "cat <<EOF >${ReadsFileStrainX}";
  cat $SLURM_SUBMIT_DIR/${ReadsFileTemplate};
  echo "EOF"; 
) > temp.txt
. temp.txt # Source the temp.txt file, to actually create the desired reads file!

#rm temp.txt # The temporary file can be deleted since we do not need it anymore!

# Declare an array with all the kmers that we want kmergenie to search for
kmers=( 65 85 97 101 119 125 133 155 187 )

# Run kmergenie for each item in the array. 
for i in "${kmers[@]}"
do
${KmergenieDir}/kmergenie ${ReadsFileStrainX} -k $i -t $SLURM_NTASKS -o ${StrainX}_${DateRun}_kmer_${i}.hist
done

rm *.histo *.pdf *.dat

echo
echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
