#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
# #SBATCH --mem-per-cpu=6400
# #SBATCH --mem=128000

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

	This script takes as input one argument. It runs from the master folder of BWA and then we change directory to the directory of the Strain we are interested in.	
	For Strain01 the correct usage is: 
	./$0 Strain01

	This is BWA to retrieve genomes. We are are mapping the Prinseq filtered fastq reads to the Spades scaffolds selected by Bandage manually. 
	
	IMPORTANT
	We have already run Bandage with the results of Spades, exported the contigs of interest to a new file (named $ContigBandage -see the parameters bellow) and uploaded this file to the server.

END
}

# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
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
StrNum=${StrainX/Strain/}

ReadDir=${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering
BandageDir=${HOME}/Pseudomonas/Bandage/${StrainX}

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "05"|"08" ]]; then
     DateRun=Nov18a # This has to be changed to "Jun18" or "Nov18a", and run this script twice, one for each case.
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi
    

Paired1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
Paired2=${StrX}_${DateRun}_prinseq_good_R_2.fastq
Single1=${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq
Single2=${StrX}_${DateRun}_prinseq_good_R_2_singletons.fastq

Bam_PE_File=${StrX}_${DateRun}_PE_BandageSpades_sorted.bam
BAM_SE1_File=${StrX}_${DateRun}_SE1_BandageSpades_sorted.bam
BAM_SE2_File=${StrX}_${DateRun}_SE2_BandageSpades_sorted.bam

ContigBandage=${StrX}_genome1_Bandage_filtered.fasta
FqOnlyGenomePE1=${StrX}_${DateRun}_PE1_ONLY_genome.fastq
FqOnlyGenomePE2=${StrX}_${DateRun}_PE2_ONLY_genome.fastq
FqOnlyGenomeSE1=${StrX}_${DateRun}_SE1_ONLY_genome.fastq
FqOnlyGenomeSE2=${StrX}_${DateRun}_SE2_ONLY_genome.fastq

#--------------------------SANITY CHECKS -------------------------------------------------------------------#
if [[ ! -s ${BandageDir}/${ContigBandage} ]]; then
    echo "The filtered contigs by Bandage are missing. Please upload them to directory ${BandageDir} and name them as ${ContigBandage}." >&2
    generalInfo >&2
    exit 1
fi



echo
echo "Which bwa?"
which bwa
echo
echo
samtools --version
echo
bedtools --version
echo

cd ${StrainX}

ln -s ${BandageDir}/${ContigBandage}
echo
echo "		BWA Analysis Begins"
echo
# Create bwa index for Pilon
bwa index ${ContigBandage}

# Align paired end reads to genome contigs with bwa mem.
# Turn sam to bam BUT output only mapped reads 
# For flag expanation go to https://broadinstitute.github.io/picard/explain-flags.html
bwa mem -t $SLURM_NTASKS ${ContigBandage} ${ReadDir}/${Paired1} ${ReadDir}/${Paired2} | \
samtools view -hu -F 4 -@ $SLURM_NTASKS - | \
samtools sort -n -@ $SLURM_NTASKS - -o ${Bam_PE_File}

echo


# Repeat the alignments for the Single ends.
bwa mem -t $SLURM_NTASKS ${ContigBandage} ${ReadDir}/${Single1} | \
samtools view -hu  -F 4 -@ $SLURM_NTASKS - -o ${BAM_SE1_File}


bwa mem -t $SLURM_NTASKS ${ContigBandage} ${ReadDir}/${Single2} | \
samtools view -hu -F 4 -@ $SLURM_NTASKS - -o ${BAM_SE2_File}

echo "-------------------------------------BWA ANALYSIS COMPLETE!!!---------------------------------------------------------"

# Bedtools bamtofastq need to have sorted bam files and for this reason we have to sort the paired aligned bam files.
echo
echo
echo "			Create the fastq files"
# Get fastq reads from bam alignment

# For the mapped reads genome fastq output!!)
bedtools bamtofastq -i ${Bam_PE_File} \
			-fq ${FqOnlyGenomePE1} \
			-fq2 ${FqOnlyGenomePE2}

samtools fastq -0 ${FqOnlyGenomeSE1} ${BAM_SE1_File}
samtools fastq -0 ${FqOnlyGenomeSE2} ${BAM_SE2_File}

# Remove useless (in my case!!) intermediate files
rm -r *fasta* *.bam

echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
