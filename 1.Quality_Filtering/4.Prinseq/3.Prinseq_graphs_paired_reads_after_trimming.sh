#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="Prinseq_after_trimming"
#SBATCH --output=Prinseq_after_trimming_job_%j.out

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
	This script runs from the master folder of Prinseq and then we change directory to the folder of each specific Strain folder.

	
	This is the Prinseq pipeline part 3.
	Here we are trying to graph the filtered fastq reads (by Prinseq in the previous step of this pipeline). We have to manually inspect the graph to see whether we are happy or not with the filtering outcome.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] && [[ ! $1 =~ Strain[0-9]{2} ]]; then
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

PrinseqDir=${HOME}/Software/prinseq-lite-0.20.4
FastpDir=${HOME}/Pseudomonas/Fastp/${StrainX}
FilteredReadsDir=${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering

DateRun=Nov18b
PrinseqFilterdR1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
PrinseqFilterdR2=${StrX}_${DateRun}_prinseq_good_R_2.fastq

#==================================================================================================

cd ${StrainX}

LC_ALL=C ${PrinseqDir}/prinseq-lite.pl -fastq ${FilteredReadsDir}/${PrinseqFilterdR1} -fastq2 ${FilteredReadsDir}/${PrinseqFilterdR2} -graph_data ${StrX}_${DateRun}_prinseq_PE_after_trimming.gd -out_good null -out_bad null -log ${StrX}_${DateRun}_prinseq_after_trimming.log

LC_ALL=C ${PrinseqDir}/prinseq-graphs-noPCA.pl -i ${StrX}_${DateRun}_prinseq_PE_after_trimming.gd -html_all

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
