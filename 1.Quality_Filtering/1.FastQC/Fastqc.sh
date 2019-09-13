#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="FastQC"
#SBATCH --output=FastQC_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END

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
	This script runs from the master folder of FastQC and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	Here we take the raw fastq files, taken from the illumina machine, and check their quality.

	We have already created soft links of each raw fastq paired reads inside the correct directory

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] && [[ ! $1 =~ Strain[0-9]{2} ]]; then
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

FastQCDir=${HOME}/Software/FastQC
Output=${HOME}/Pseudomonas/FastQC/Before_filtering/${StrainX}
RawReadsDir=${HOME}/RawReads_Nov2018b
DateRun=Nov18b
Pair1Fastq=${StrX}_R1_${DateRun}.fastq
Pair2Fastq=${StrX}_R2_${DateRun}.fastq


cd ${StrainX}
echo
mkdir tmp
echo
${FastQCDir}/fastqc --version
echo

LC_ALL=C ${FastQCDir}/fastqc --threads $SLURM_NTASKS --noextract --dir tmp -o ${Output} ${RawReadsDir}/${Pair1Fastq} ${RawReadsDir}/${Pair2Fastq} 


<< ////
	
	NOTE:
	-t --threads    Specifies the number of files which can be processed
                        simultaneously.  Each thread will be allocated 250MB of
                        memory so you shouldn't run more threads than your
                        available memory will cope with, and not more than
                        6 threads on a 32 bit machine

////


echo
echo "==========================================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
