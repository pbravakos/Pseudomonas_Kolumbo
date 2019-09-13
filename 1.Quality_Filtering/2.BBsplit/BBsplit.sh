#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000
#SBATCH --job-name="BBsplit"
#SBATCH --output=BBsplit_job_%j.out
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
	This script runs from the master folder of BBsplit and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	BBSplit is a tool that bins reads by mapping to multiple references simultaneously, using BBMap.
	We use it to filter out all the PhiX sequences we have, and keep only the non phix fastq reads.
	
	For this reason we have already downloaded the phiX, from Illumina iGenomes found here: http://emea.support.illumina.com/sequencing/sequencing_software/igenome.html
	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz

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
BBMaptDir=${HOME}/Software/bbmap/
PhiXFasta=${HOME}/Pseudomonas/BBsplit/PhiX_genome.fa
RawReadsDir=${HOME}/RawReads_Nov2018b
DateRun=Nov18b
Pair1Fastq=${StrX}_R1_${DateRun}.fastq
Pair2Fastq=${StrX}_R2_${DateRun}.fastq


cd ${StrainX}

${BBMaptDir}/bbsplit.sh in=${RawReadsDir}/${Pair1Fastq}  in2=${RawReadsDir}/${Pair2Fastq} ref_PhiX=${PhiXFasta} ambiguous=toss basename=${StrainX}_${DateRun}_contaminated_%.fq outu1=${StrainX}_${DateRun}_R1_NOphix.fastq outu2=${StrainX}_${DateRun}_R2_NOphix.fastq

<< ////

	OTHER WAYS TO RUN THE PROGRAM
	
	java -Xmx3g -cp /home/pbravakos/Software/BBMap/sh/current/ align2.BBSplitter
	/home/pbravakos/Software/BBMap/sh/bbsplit.sh
	/home/pbravakos/Software/bbmap_old/bbsplit.sh


	LOCAL PC
	java -Xmx3g -cp /home/panos/Downloads/BBMap/sh/current align2.BBSplitter in=/mnt/020490155FD7A1A7/Voula_Pseudomonas/Analysis/1_S1_L001_R1_001.fastq in2=/mnt/020490155FD7A1A7/Voula_Pseudomonas/Analysis/1_S1_L001_R2_001.fastq ref_Phix=/mnt/020490155FD7A1A7/Voula_Pseudomonas/Analysis/PhiX_Illumina_RTA/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa basename=Strain07_June18_contaminated_%.fq outu1=Strain07_June18_Read1_NOphix.fq outu2=Strain07_June18_Read2_NOphix.fq

////

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


