#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="CLA"
#SBATCH --output=CLA_job_%j.out

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
echo

generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	$0 Strain01
	This script runs from the master folder of CLA and then we change directory to the folder of each specific Strain folder.
	
	IMPORTANT!
	I have modified the program in order to run most of the sub-programs like bwa, blast and bowtie2 (in place of the default bowtie) from the path!!!

	NOTE:
	Scaffolds fasta files used here as input is the Spades scaffolds output.
	References were selected by blast search of the Barnap 16S and 23S results on NCBI nr database.

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

SpadesDir=${HOME}/Pseudomonas/Spades/${StrainX}/1sttrial
InputScaff=scaffolds.fasta

RefDir=${HOME}/Pseudomonas/References/${StrainX}

CLADir=${HOME}/Software/C-L-Authenticator_v.1.0_Linux

# We will use the Prinseq fastq files here, in order to have as much (filtered!) reads as possible.
PrinseqReadsDir=${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16"|"05"|"08" ]] ; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

PE1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
PE2=${StrX}_${DateRun}_prinseq_good_R_2.fastq

# Next we will assign the correct Reference to each of our Strains.
if [[ $StrNum = "09" ]]; then 
    RefFasta=Pseudomonas_sp_MT-1.fasta
elif [[ $StrNum = "10" ]]; then 
    RefFasta=Pseudomonas_sp_MT-1.fasta
elif [[ $StrNum = "11" ]]; then 
    RefFasta=Pseudomonas_balearica_DSM_6083.fasta
elif [[ $StrNum = "12" ]]; then 
    RefFasta=Pseudomonas_sp_R2A2.fasta
elif [[ $StrNum = "13" ]]; then 
    RefFasta=Pseudomonas_stutzeri_CCUG_29243.fasta
elif [[ $StrNum = "14" ]]; then 
    RefFasta=Pseudomonas_stutzeri_1W1-1A.fasta
elif [[ $StrNum = "16" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_MTB-1.fasta
elif [[ $StrNum = "18" ]]; then 
    RefFasta=Pseudomonas_sp_R2A2.fasta
elif [[ $StrNum = "19" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_MTB-1.fasta
elif [[ $StrNum = "20" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_strain_NCTC9433.fasta
elif [[ $StrNum = "21" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_MTB-1.fasta
elif [[ $StrNum = "22" ]]; then 
    RefFasta=Pseudomonas_balearica_DSM_6083.fasta
elif [[ $StrNum = "23" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_strain_NCTC9433.fasta
elif [[ $StrNum = "24" ]]; then 
    RefFasta=Pseudomonas_stutzeri_1W1-1A.fasta
elif [[ $StrNum = "25" ]]; then 
    RefFasta=Pseudomonas_stutzeri_CCUG_29243.fasta
elif [[ $StrNum = "05" ]]; then 
    RefFasta=Pseudomonas_aeruginosa_strain_AR_455.fasta
elif [[ $StrNum = "08" ]]; then 
    RefFasta=Pseudomonas_xanthomarina_strain_LMG_23572.fasta
else
    echo "$StrainX fasta reference file has not been set. Please check the reference list again. " >&2
    exit 1
fi


# Parameters of CLA that can be changed
# Two different libraries we constructed for the two separate Miseq runs. 
# IMPORTANT
# Strains in the same Miseq run (i.e. $DateRun) can potentially be from different libraries!
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"18"|"05" ]]; then
    InsertSize=583
elif [[ $StrNum =~ "16"|"19"|"20"|"21"|"22"|"23"|"24"|"25"|"08" ]]; then
     InsertSize=612
else
    echo "Insert Size for $StrainX has not been set. Please check the Spades output log and update the configuration on this script."  >&2
    exit 1
fi

#===========================================================================================================================================

cd $StrainX


# Create links for input files:
# 1. Two read files in fastq format generated by paired end sequencing
# 2. Contig file from any of the desired denovo assemblers
# 3. Reference genome in fasta format
ln -s ${PrinseqReadsDir}/$PE1
ln -s ${PrinseqReadsDir}/$PE2
ln -s ${SpadesDir}/${InputScaff}
ln -s ${RefDir}/${RefFasta}

# Run the program!!
LC_ALL=C ${CLADir}/C-L-Authenticator -i ${InputScaff} -r1 $PE1 -r2 $PE2 -ref ${RefFasta} -f 1 -ins $InsertSize

rm CLA-Results/j*

# Repeat the CLA run, but this time with the fastq files from the second Miseq run and input scaffolds from the first CLA run. This applies only for strains that have a second miseq run!
if [[ $StrNum =~ "05"|"08" ]]; then
     DateRun=Jun18
     InsertSize=610
     InputScaff=$SLURM_SUBMIT_DIR/${StrainX}/CLA-Results/final_Scaffolds.fa
     OutputDir=CLA-Results-2nditeration
     ln -s ${PrinseqReadsDir}/$PE1
     ln -s ${PrinseqReadsDir}/$PE2
     # Run the program again, with the new parameters!!
     LC_ALL=C ${CLADir}/C-L-Authenticator -i ${InputScaff} -r1 $PE1 -r2 $PE2 -ref ${RefFasta} -f 1 -ins $InsertSize -o $OutputDir
     rm ${OutputDir}/j*
fi


echo "==============================="
echo "PBS job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0


