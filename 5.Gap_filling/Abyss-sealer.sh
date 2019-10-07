#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="AbyssSealer"
#SBATCH --output=AbyssSealer_job_%j.out
#SBATCH --mail-user=pbravakos@hcmr.gr

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

	This script takes as input one argument. It runs from the master folder of Abyss Sealer and then we change directory to the directory of each specific Strain.	
	For Strain 01 the correct usage is: 
	$0 Strain01

	NOTE:
	We are using the pre-filtered reads here! These are the ones right after fastp and before prinseq. This is why trim-quality is required as an option.

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

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16"|"05"|"08" ]] ; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

AbyssDir=${HOME}/Software/abyss
KrakenDir=${HOME}/Pseudomonas/Kraken/${StrainX}
ReadDir=${HOME}/Pseudomonas/Fastp/${StrainX}
Read1=${StrX}_${DateRun}_R1_fastp.fastq
Read2=${StrX}_${DateRun}_R2_fastp.fastq
ContigKraken=${StrX}_CLA_Kraken.fasta

#---------------------------------------------------------------------------------------------------------------------------------

cd ${StrainX}


${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=25 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k25.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=39 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k39.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=49 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k49.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=59 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k59.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=69 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k69.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=81 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k81.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=91 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k91.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=101 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k101.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=117 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k117.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

${AbyssDir}/Sealer/abyss-sealer -vv -k25 -k39 -k49 -k59 -k69 -k81 -k91 -k101 -k117 -o ${StrX}_CLA_Sealer -S ${KrakenDir}/${ContigKraken} -i k25.bloom -i k39.bloom -i k49.bloom -i k59.bloom -i k69.bloom -i k81.bloom -i k91.bloom -i k101.bloom -i k117.bloom --max-gap-length=1500 --max-branches=3000 --max-paths=500 --fix-errors --threads=$SLURM_NTASKS --trim-quality=16 


#rm *.bloom

# Repeat the process for the Strains that we have sequenced twice, with the second library reads!
if [[ $StrNum =~ "05"|"08" ]]; then

   rm *.bloom
   mv ${StrX}_CLA_Sealer_scaffold.fa ${StrX}_CLA_Sealer_1st_scaffold.fa

   DateRun=Jun18
   Read1=${StrX}_${DateRun}_R1_fastp.fastq
   Read2=${StrX}_${DateRun}_R2_fastp.fastq

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=25 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k25.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=39 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k39.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=49 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k49.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=59 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k59.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=69 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k69.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=81 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k81.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=91 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k91.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=101 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k101.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Bloom/abyss-bloom build -vv --kmer=117 --bloom-size=20G --threads=$SLURM_NTASKS --levels=3 k117.bloom --trim-quality=16 ${ReadDir}/${Read1} ${ReadDir}/${Read2}

   ${AbyssDir}/Sealer/abyss-sealer -vv -k25 -k39 -k49 -k59 -k69 -k81 -k91 -k101 -k117 -o ${StrX}_CLA_Sealer -S ${StrX}_CLA_Sealer_1st_scaffold.fa -i k25.bloom -i k39.bloom -i k49.bloom -i k59.bloom -i k69.bloom -i k81.bloom -i k91.bloom -i k101.bloom -i k117.bloom --max-gap-length=1500 --max-branches=3000 --max-paths=500 --fix-errors --threads=$SLURM_NTASKS --trim-quality=16 
   
fi

echo "============================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
