#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Spades"
#SBATCH --output=Spades_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

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

	This script takes as input one argument. It runs from the master folder of Spades and then we change directory to the directory of each specific Strain.	
	For Strain 01 the correct usage is: 
	./$0 Strain01

END
}

# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ Strain[0-9]{2} ]]; then
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi

# Check that the working directory exists.
[[ ! -d ${1} ]] && { echo "Directory ${1} cannot be not found!" >&2; exit 1; }

# Start a job dependency to move the Sdtout file to the correct folders
sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1

# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
Iteration=1sttrial  # This parameter can be changed every time we want to have more than one Spades runs saved on separate directories. 

SpadesDir=${HOME}/Software/SPAdes-3.13.0-Linux/bin
OutputDir=$SLURM_SUBMIT_DIR/${StrainX}/${Iteration}


YamlFileTemplate=Reads.yaml
YamlFileStrainX=${StrainX}_reads.yaml

# Fist, create the yaml file with the correct reads, for each 
cat > $YamlFileTemplate <<"EOF"

- "left reads":
  - "${HOME}/Pseudomonas/BBmerge/${StrainX}/${StrX}_${DateRun}_Read1_unmerged.fastq"
  "merged reads":
  - "${HOME}/Pseudomonas/BBmerge/${StrainX}/${StrX}_${DateRun}_merged.fastq"
  "orientation": "fr"
  "right reads":
  - "${HOME}/Pseudomonas/BBmerge/${StrainX}/${StrX}_${DateRun}_Read2_unmerged.fastq"
  "single reads":
  - "${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering/${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq"
  - "${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering/${StrX}_${DateRun}_prinseq_good_R_2_singletons.fastq"

  "type": "paired-end"
EOF





#DateRun=Nov18a  # This parameter is needed for the bash substitution in the reads yaml file. Not to be deleted!
# Kmer selection for each Strain based on the kmergenie results.
#if [[ $StrainX = Strain09 ]]; then
#    kmers="29,49,61,85,95,105,113" 
#    echo
#    echo "Kmers for Strain09 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain10 ]]; then
#    kmers="41,51,65,73,95,105,113"
#    echo
#    echo "Kmers for Strain10 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain11 ]]; then
#    kmers="27,51,59,65,79,83,109"
#    echo
#    echo "Kmers for Strain11 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain12 ]]; then
#    kmers="37,39,43,51,59,75,93,103"
#    echo
#    echo "Kmers for Strain12 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain13 ]]; then
#    kmers="35,55,61,71,85,89,93,101,117"
#    echo
#    echo "Kmers for Strain13 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain14 ]]; then
#    kmers="29,41,51,63,65,69,85,93,99,105"
#    echo
#    echo "Kmers for Strain14 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain16 ]]; then
#    kmers="25,41,65,79,87,97,115"
#    echo
#    echo "Kmers for Strain16 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain05 ]]; then   # Strain05 for Nov18a fastq files
#    kmers="35,49,59,67,75,81,85,105"
#    echo
#    echo "Kmers for Strain05 Nov18 selected for ${Iteration}!!"
#    echo
#elif [[ $StrainX = Strain08 ]]; then   # Strain08 for Nov18a fastq files
#    kmers="21,51,67,69,77,99,107,115"
#    echo
#    echo "Kmers for Strain08 Nov18 selected for ${Iteration}!!"
#    echo
#else
#    echo "kmers for $StrainX have not been set. Please run first kmergenie, and then select the kmers to be used as input for Spades."
#    exit 1
#fi

DateRun=Nov18b  # This parameter is needed for the bash substitution in the reads yaml file. We declare it here, to indicate that the next block of code is for the Nov18b Miseq Run.
if [[ $StrainX = Strain18 ]]; then
    kmers="27,39,49,59,67,75,79,83,97,101" 
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain19 ]]; then
    kmers="25,27,29,31,35,57,67,79,105,113"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain20 ]]; then
    kmers="39,47,55,57,71,85,95,103,107,109"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain21 ]]; then
    kmers="35,45,53,61,63,67,73,81,85,101,109"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain22 ]]; then
    kmers="31,45,55,65,69,79,89,101,103,113,117"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain23 ]]; then
    kmers="35,45,55,65,67,77,85,91,99,109"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain24 ]]; then
    kmers="35,45,55,57,63,75,85,89,97,101"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain25 ]]; then   # Strain08 for Nov18a fastq files
    kmers="31,41,51,61,71,85,97,107"
    echo
    echo "Kmers for $StrainX Nov18 selected for ${Iteration}!!"
    echo
else
    echo "kmers for $StrainX have not been set. Please run first kmergenie, and then select the kmers to be used as input for Spades." >&2
    exit 1
fi



#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that the Yaml Template file has been uploaded to the server
if [[ ! -s $SLURM_SUBMIT_DIR/${YamlFileTemplate} ]]; then
    echo "Yaml file is missing!" >&2
    exit 1
fi

# Check the output folder exists and if not create it! 
if [[ ! -d $SLURM_SUBMIT_DIR/${StrainX}/${Iteration} ]]; then 
    mkdir $SLURM_SUBMIT_DIR/${StrainX}/${Iteration}
fi

#-------------------------------------------------------------------------------------------------------------

cd ${StrainX}/${Iteration}


# Create a new yaml file with bash parameter substitution.
( echo "cat <<EOF >$SLURM_SUBMIT_DIR/${StrainX}/${Iteration}/${YamlFileStrainX}";
  cat $SLURM_SUBMIT_DIR/${YamlFileTemplate};
  echo "EOF"; 
) > $SLURM_SUBMIT_DIR/${StrainX}/${Iteration}/temp.yml
. $SLURM_SUBMIT_DIR/${StrainX}/${Iteration}/temp.yml # Source the temp file, to actually create the desired yaml file!

rm temp.yml  # The temp file is no longer needed and can be deleted.



python3 ${SpadesDir}/spades.py --only-assembler --dataset ${YamlFileStrainX} --careful -k $kmers -m $SLURM_MEM_PER_NODE -t $SLURM_NTASKS -o ${OutputDir}


# NEXT STEPS:
# Open fastg output with Bandage to check for contamination, and either continue with scaffolding, or start the optional decontamination pipeline.

echo "==============================="
echo "SBATCH job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0

