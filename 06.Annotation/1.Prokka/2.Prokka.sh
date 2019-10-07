#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Prokka"
#SBATCH --output=Prokka_job_%j.out

# for calculating the amount of time the job takes
begin=`date +%s`

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
		
	This script runs from the master folder of Prokka and then we change directory to the folder of each specific Strain folder.
	We assume that in the master folder there are already folders named Strain01 Strain02 etc.

	IMPORTANT!!
	tbl2asn is a dependency for Prokka and has to be updated every year! 
	Check before running Prokka that you have the latest version on the path! 

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

ProkkaBin=${HOME}/Software/prokka/bin
ReadsDir=${HOME}/Pseudomonas/Seqtk/${StrainX}

# Parameters that could be changed
Genus=Pseudomonas
Kingdom=Bacteria
Gram=neg
Gcode=11
ContigLen=200
Eval=1e-9

#-------------------------------------------------------------------------
# Print the Prokka version we are going to use in the analysis.
${ProkkaBin}/prokka --version
echo
cd $StrainX

# Create a soft link of the contigs fasta file
ln -s ${ReadsDir}/${StrX}_Pilon_CLA_Blast.fasta

LC_ALL=C ${ProkkaBin}/prokka --force \
			     --prefix ${StrX}_Prokka \
			     --outdir ${StrX}_Prokka \
			     --genus ${Genus} \
			     --kingdom ${Kingdom} \
			     --gcode ${Gcode} \
			     --gram ${Gram} \
			     --usegenus \
			     --evalue ${Eval} \
			     --mincontiglen ${ContigLen} \
			     --cpus $SLURM_NTASKS \
			     --rfam \
			     ${StrX}_Pilon_CLA_Blast.fasta


echo "============================================"
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
