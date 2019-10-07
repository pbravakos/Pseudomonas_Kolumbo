#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
# The default limit is 3000MB per core.
#SBATCH --job-name="BBnorm"
#SBATCH --output=BBnorm_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr


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

	This script takes as input one argument. It runs from the master folder of BBnorm and then we change directory to the directory of each specific Strain.	
	For Strain 01 the correct usage is: 
	scipt.sh Strain01

	NOTE:
	http://seqanswers.com/forums/showthread.php?t=49763 http://seqanswers.com/forums/showthread.php?t=49763&page=2 http://seqanswers.com/forums/showthread.php?t=49763&page=3
	I recommend pre-processing (adapter trimming, contaminant removal, quality-trimming or filtering) prior to normalization, because those processes all remove spurious kmers that make it harder to determine read depth, and thus improve the normalization results. 

END
}

# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[  ! $1 =~ Strain[0-9]{2} ]]; then
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi

# Check that the working directory exists.
[[ ! -d ${1} ]] && { echo "Directory ${1} cannot be not found!" >&2; exit 1; }

# Start a job dependency to move the Sdtout file to the correct folders
${HOME}/sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1

# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
BBMapDir=${HOME}/Software/bbmap

DateRun=Nov18b
PrinseqFilteredDir=${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering

PrinseqR1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
PrinseqR2=${StrX}_${DateRun}_prinseq_good_R_2.fastq

# Parameters for BBnorm that can change
TmpDir=tmp
PreFilter=f # True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable. Whether or not to use "prefilter" just depends on the amount of memory you have rather than the workflow. It basically makes BBNorm take twice as long but increases accuracy in cases where you have a very large dataset compared to memory - so, there's no penalty for using it, and it always increases accuracy, but the increase is trivial if you have a lot of memory. So if you have lots of ram or a small dataset you don't need it.
FixSpikes=f  
Passes=1  # 1 pass is the basic mode.  2 passes (default) allows greater accuracy, error detection, better contol of output depth.  The first pass will normalize to some depth higher than the ultimate desired depth, and the second pass will normalize to the target depth. This allows, in the first pass, preferential discarding of reads that are low quality. So the result from the second pass should still be a target of $Target.
Target=999999999  # Target normalization depth.  NOTE:  All depth parameters control kmer depth, not read depth.
Min=1 # This will pass only reads with a depth of at least $Min to “out”, and low-depth reads under $Min to “outt” (outtoss).
Kmer=70 # The smaller the kmer, the less reads we filter out.

# Next we will deal with special cases that need more strict filtering i.e. more reads being filtered out.
if [[ ${StrainX} = Strain18 ]]; then
    Min=2
    Kmer=105
    echo
    echo "${StrainX} normalization is running with parameters: Minimum Passes = $Min and Kmer length = $Kmer"
    echo
fi

#==================================================================================================

cd $StrainX

${BBMapDir}/bbnorm.sh in=${PrinseqFilteredDir}/${PrinseqR1} in2=${PrinseqFilteredDir}/${PrinseqR2} out=${StrX}_${DateRun}_R1_bbnorm.fastq out2=${StrX}_${DateRun}_R2_bbnorm.fastq outt=${StrX}_${DateRun}_removed_reads_bbnorm.fastq tmpdir=$TmpDir prefilter=$PreFilter fixspikes=$FixSpikes passes=$Passes target=$Target min=$Min k=$Kmer hist=${StrX}_${DateRun}_kmer_depth.hist histout=${StrX}_${DateRun}_kmer_removed_reads_depth.hist t=$SLURM_NTASKS -Xmx${SLURM_MEM_PER_NODE}m 


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

