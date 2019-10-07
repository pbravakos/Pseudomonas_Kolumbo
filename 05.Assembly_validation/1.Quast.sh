#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
# #SBATCH --mem=128000
#SBATCH --mem-per-cpu=6400
# Memory per node specification is in MB. It is optional.
#SBATCH --job-name="Quast"
#SBATCH --output=Quast_job_%j.out
# #SBATCH --error=Sderr_Quast_job_%j.txt
#SBATCH --mail-user=pbravakos@hcmr.gr
#SBATCH --mail-type=FAIL,END
#SBATCH --no-requeue

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
	./$0 Strain01
	This script runs from the master folder of Pilon and then we change directory to the folder of each specific Strain folder.
	
	NOTE:
	Here we use the fastq files after Prinseq even if further downstream filtering has been done (e.g. contamination removal) because we want as much information as possible to correct our contigs!!

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
QuastDir=${HOME}/Software/quast-5.0.2
ReadDir=${HOME}/Pseudomonas/Reads/${StrainX}
SpadesDir=${HOME}/Pseudomonas/Spades/${StrainX}/1sttrial
# CeleraDir=${HOME}/Pseudomonas/Celera/${StrainX}/9-terminator
PilonDir=${HOME}/Pseudomonas/Pilon/${StrainX}
RefDir=${HOME}/Pseudomonas/References/${StrainX}
ContigDir=${HOME}/Pseudomonas/Seqtk/${StrainX}

if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16"|"05"|"08" ]]; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "01"|"02"|"03"|"04"|"06"|"07" ]]; then
    DateRun=Jun18
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

Paired1=${StrX}_R1_${DateRun}.fastq
Paired2=${StrX}_R2_${DateRun}.fastq
Single1=${StrX}_SE1_${DateRun}.fastq
Single2=${StrX}_SE2_${DateRun}.fastq
MeReads=${StrX}_MeReads_${DateRun}.fastq

RefGenome=`find ${RefDir} -name '*.fasta'`
RefFeatures=`find ${RefDir} -name '*.gff3'`


#------------------------------------------------------------------------------------------------------------

# Change directory to the Strain specific folder
cd ${StrainX}

echo "			This is the QUAST analysis for ${StrainX}"

# Run Quast!!
LC_ALL=C python3 ${QuastDir}/quast.py -r ${RefGenome} --min-contig 200 -o ${StrX}_results --threads $SLURM_NTASKS --split-scaffolds --features ${RefFeatures} --labels ${StrX}_CLA,${StrX}_Spades --k-mer-stats --gene-finding --conserved-genes-finding --ambiguity-usage all --rna-finding --pe1 ${ReadDir}/${Paired1} --pe2 ${ReadDir}/${Paired2} --single ${ReadDir}/${Single1} --single ${ReadDir}/${Single2} --single ${ReadDir}/${MeReads} ${PilonDir}/${StrX}_Pilon_CLA.fasta ${SpadesDir}/scaffolds.fasta

echo "======================================================================"
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
