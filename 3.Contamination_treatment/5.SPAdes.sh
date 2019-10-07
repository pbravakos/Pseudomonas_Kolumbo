#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Spades"
#SBATCH --output=Spades_job_%j.out


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
	$0 Strain01

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
StrNum=${StrainX/Strain/}

Iteration=1sttrial  # This parameter can be changed every time we want to have more than one Spades runs saved on separate directories. 1sttrial is reserved by the first run we did with the reads right after the filtering process.

SpadesDir=${HOME}/Software/SPAdes-3.13.0-Linux/bin
OutputDir=$SLURM_SUBMIT_DIR/${StrainX}/${Iteration}

# Next we wil define the $DateRun parameter based on the Strain number input, which is needed for the reads file creation.
if [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "05"|"08" ]]; then
    DateRun=Nov18a
    DateRun2=Jun18  # These Strains have been sequenced more than one time, and we have fastq files from each run.
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

YamlFileTemplate1=Reads_after_BWA.yaml
YamlFileTemplate2=Reads_after_BWA_from_different_runs.yaml

YamlFileStrainX=${StrainX}_reads.yaml

# Next, we will assign the $kmer parameter a value based on the kmergenie results.
if [[ $StrainX = Strain18 ]]; then
    kmers="" 
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain19 ]]; then
    kmers="27,31,35,45,55,63,73,83,93,103,111"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain25 ]]; then   
    kmers="31,41,51,63,67,77,85,91,95,105,119"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain05 ]]; then   
    kmers="29,39,49,59,69,77,81,91,101,109"
    echo
    echo "Kmers for $StrainX selected for ${Iteration}!!"
    echo
elif [[ $StrainX = Strain08 ]]; then
    kmers="29,35,45,59,65,75,85,95,105,111"
else
    echo "kmers for $StrainX have not been set. Please run first kmergenie, and then select the kmers to be used as input for Spades."
    exit 1
fi



#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that the Yaml Template file has been uploaded to the server
if [[ ! -s $SLURM_SUBMIT_DIR/${YamlFileTemplate1} ]]; then
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
  if [[ $StrNum =~ "05"|"08" ]]; then
      cat $SLURM_SUBMIT_DIR/${YamlFileTemplate2};
  else
      cat $SLURM_SUBMIT_DIR/${YamlFileTemplate1};
  fi
  echo "EOF"; 
) > $SLURM_SUBMIT_DIR/${StrainX}/${Iteration}/temp.yml
. $SLURM_SUBMIT_DIR/${StrainX}/${Iteration}/temp.yml # Source the temp file, to actually create the desired yaml file!

#rm temp.yml  # The temp file is no longer needed and can be deleted.



python3 ${SpadesDir}/spades.py --only-assembler --dataset ${YamlFileStrainX} --careful -k $kmers -m $SLURM_MEM_PER_NODE -t $SLURM_NTASKS -o ${OutputDir}


# NEXT STEPS:
# Open fastg output with Bandage to check for contamination, and either continue with the next steps of pipeline or go back to the filtering steps to reduce contamination.

echo "==============================="
echo "SBATCH job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0

