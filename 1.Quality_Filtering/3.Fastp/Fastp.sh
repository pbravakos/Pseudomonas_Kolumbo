#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=16
# #SBATCH --mem=128000
#SBATCH --job-name="Fastp"
#SBATCH --output=Fastp_job_%j.out

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
	This script runs from the master folder of Fastp and then we change directory to the folder of each specific Strain folder.

	Attention:
	Fastp can take up to 16 threads no more. If more than 16 threads are inserted, it runs normally but uses only the 16 threads.
	
	NOTE:
	Fastp cannot produce singleton reads, meaning, it filters out either both reads or none of them. In order to be able to produce singleton reads we have to repeat 

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

# Start a job dependency to move the Sdtout file to the correct folders
${HOME}/sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}

FastpDir=${HOME}/Software/fastp-0.19.5
BBSplitDir=${HOME}/Pseudomonas/BBsplit/${StrainX}

DateRun=Nov18b
BBSplitR1=${StrainX}_${DateRun}_R1_NOphix.fastq
BBSplitR2=${StrainX}_${DateRun}_R2_NOphix.fastq

# Fastp parameters that can be changed
Adapter1="GATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # The Illumina adapter1
Adapter2="GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # The Illumina adapter2
TrimFront1=11  # Trimm that number of bases in the front of the read 1
TrimFront2=11  # Trimm that number of bases in the front of the read 2
TrimTail1=10  # Trimm that number of bases in the tali of the read 1
TrimTail2=58  # Trimm that number of bases in the tail of the read 2
CutMeanQuality=17  # When triiming by a mean, this is the mean phred quality score that will be kept in the sliding window
QualifyBaseQuality=15  # The phred quality score that determines whether a base will be characterized as qualified or not.
UnqualPercLim=40   # Percentage of reads that are allowed to be unqualified in a read passing the filters.
NBaseLimit=5   # This number is the phred score which makes the whole read to be discarded, if something less than this is found in any read.
LengthRequired=40


cd ${StrainX}

echo
echo "Fastp Version:"
${FastpDir}/fastp --version
echo

ln -s ${BBSplitDir}/${BBSplitFastq}
ln -s ${BBSplitDir}/${BBSplitR1}
ln -s ${BBSplitDir}/${BBSplitR2}


${FastpDir}/fastp --in1 ${BBSplitR1} --in2 ${BBSplitR2} --out1 ${StrX}_${DateRun}_R1_fastp.fastq --out2 ${StrX}_${DateRun}_R2_fastp.fastq --adapter_sequence ${Adapter1} --adapter_sequence_r2 ${Adapter2} --trim_front1 ${TrimFront1} --trim_tail1 ${TrimTail1} --trim_front2 ${TrimFront2} --trim_tail2 ${TrimTail2} --disable_trim_poly_g --cut_by_quality3 --cut_mean_quality ${CutMeanQuality} --qualified_quality_phred ${QualifyBaseQuality} --unqualified_percent_limit ${UnqualPercLim} --n_base_limit ${NBaseLimit} --length_required ${LengthRequired} --correction --overrepresentation_analysis --html ${StrX}_${DateRun}_Fastp.html --json ${StrX}_${DateRun}_Fastp.json --report_title "${StrX} ${DateRun} Fastp report" --thread $SLURM_NTASKS


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
