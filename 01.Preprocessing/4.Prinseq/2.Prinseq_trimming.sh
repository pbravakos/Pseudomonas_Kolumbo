#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="Prinseq_trimming"
#SBATCH --output=Prinseq_trimming_job_%j.out

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

	
	This is the Prinseq pipeline no 2.
	Here we are actually trimming the fastq reads, based on the results of the first part of this pipeline.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ Strain[0-9]{2} ]]; then
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

PrinseqDir=${HOME}/Software/prinseq-lite-0.20.4
FastpDir=${HOME}/Pseudomonas/Fastp/${StrainX}

DateRun=Nov18b
FastpFastqR1=${StrX}_${DateRun}_R1_fastp.fastq
FastpFastqR2=${StrX}_${DateRun}_R2_fastp.fastq

# Prinseq parameters that can be changed
OutFormat=3   # 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL) 
MinLen=40  # Filter sequence shorter than $MinLen
MinQualMin=17  # Filter sequence with quality score mean below $MinQualMin
NsMaxN=3  # Filter sequence with more than $NsMaxN Ns
TrimQualRight=16  # Trim sequence by quality score from the 3'-end with this threshold score
TrimQualWindow=10   # The sliding window size used to calculate quality score by type
TrimQualStep=5  # Step size used to move the sliding window
TrimQualType=min  # Type of quality score calculation to use. Allowed options are min, mean, max and sum
TrimQualRule=lt  # Rule to use to compare quality score to calculated value. Allowed options are lt (less than), gt (greater than) and et (equal to)



if [[ ! -d $SLURM_SUBMIT_DIR/${StrainX}/Filtering ]]; then 
   mkdir $SLURM_SUBMIT_DIR/${StrainX}/Filtering
fi
#==================================================================================================


cd ${StrainX}/Filtering

ln -s ${FastpDir}/${FastpFastqR1}
ln -s ${FastpDir}/${FastpFastqR2}

LC_ALL=C ${PrinseqDir}/prinseq-lite.pl -fastq ${FastpFastqR1} \
					-fastq2 ${FastpFastqR2} \
					-out_format $OutFormat \
					-min_len $MinLen \
					-min_qual_mean $MinQualMin \
					-ns_max_n $NsMaxN \
					-noniupac \
					-trim_qual_right $TrimQualRight \
					-trim_qual_type $TrimQualType \
					-trim_qual_rule $TrimQualRule \
					-trim_qual_window $TrimQualWindow \
					-trim_qual_step $TrimQualStep \
					-out_good ${StrX}_${DateRun}_prinseq_good_R \
					-out_bad ${StrX}_${DateRun}_prinseq_bad_R \
					-log ${StrX}_${DateRun}_prinseq_trimming.log

# -derep 1 -derep_min 4 -lc_method entropy -lc_threshold 66 -trim_tail_right 20 -rm_header -seq_id Str07_Jun18_ -trim_qual_left 20


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
