#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="BBmerge"
#SBATCH --output=BBmerge_job_%j.out


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
	$0 Strain01
	This script runs from the master folder of BBmerge and then we change directory to the folder of each specific Strain folder.

	BBmerge will try to merge paired end sequences. For more information look https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/

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

BBMapDir=${HOME}/Software/bbmap
ReadsDir=${HOME}/Pseudomonas/BWA/${StrainX}

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "05"|"08" ]]; then
    DateRun=Jun18 # This has to be changed to "Jun18" or "Nov18a", and run this script twice, one for each case.
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

R1=${StrX}_${DateRun}_PE1_ONLY_genome.fastq
R2=${StrX}_${DateRun}_PE2_ONLY_genome.fastq

# Parameters for BBmerge that can change
MinQual=14 # Ignore bases with quality below this.
Mismatch=2  # Do not allow more than this many mismatches.

#==================================================================================================

cd ${StrainX}

${BBMapDir}/bbmerge.sh in1=${ReadsDir}/${R1} in2=${ReadsDir}/${R2} out=${StrX}_${DateRun}_after_BWA_merged.fastq outu1=${StrX}_${DateRun}_after_BWA_Read1_unmerged.fastq outu2=${StrX}_${DateRun}_after_BWA_Read2_unmerged.fastq ihist=${StrX}_${DateRun}_after_BWA_length_hist.txt ordered=t mismatches=${Mismatch} ultrastrict=t minq=${MinQual} t=$SLURM_NTASKS

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
