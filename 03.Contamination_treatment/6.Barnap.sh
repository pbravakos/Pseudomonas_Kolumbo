#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Barnap"
#SBATCH --output=Barnap_job_%j.out

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
	./0 Strain01
	This script runs from the master folder of Barnap and then we change directory to the folder of each specific Strain folder.
	As input we use here the Spades output after BWA, scaffold fasta files.

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

BarrnapDir=${HOME}/Software/barrnap/bin

InputFasta=${HOME}/Pseudomonas/Spades/${StrainX}/2ndtrial/scaffolds.fasta


cd ${StrainX}
LC_ALL=C ${BarrnapDir}/barrnap --threads $SLURM_NTASKS --outseq ${StrX}_after_BWA_ribosomal_seqs_ONLY.fasta $InputFasta


echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
