#!/bin/bash
#SBATCH --partition=minibatch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name="ElConcatenero"
#SBATCH --output=ElConcatenero_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

# Some job specific info
echo "Job ID is = " "$SLURM_JOBID"
echo "SLURM cluster name = " "$SLURM_CLUSTER_NAME"
echo "SLURM partition = " "$SLURM_JOB_PARTITION"
echo "SLURM node list = " "$SLURM_JOB_NODELIST"
echo "SLURM num of nodes = " "$SLURM_JOB_NUM_NODES"
echo "SLURM number of tasks = " "$SLURM_NTASKS"
echo "SLURM memory per node = " "$SLURM_MEM_PER_NODE"
echo "SLURM memory per cpu = " "$SLURM_MEM_PER_CPU"
echo "working directory = " "$SLURM_SUBMIT_DIR"
echo "=================================================="
echo "SBATCΗ job started " "$(date)"
echo "=================================================="
echo


generalInfo () {
    cat <<END

	This script runs in the Anvio Phylogenomics ElConcatenero master directory. In order to run it, the command should be like this:
	sbatch $0

	# IMPORTANT
	We should have finished with the alignment of each separate gene before running this program. These MSAs should be in ${MSADir}.


	ElConcatenero was downloaded from https://github.com/ODiogoSilva/ElConcatenero3

END
}

# INITIAL PARAMETERS
ElConcateneroDir=${HOME}/Software/ElConcatenero3
MSADir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Final_Selected_genes
Outgroup=Cellvibrio_japonicus_Ueda107



echo
echo "		ElConcatenero Analysis begins!!!"
echo

# First, we will create a soft link of all the MSA into our working directory and then we will create a variable with all the MSA that will be used as an input in the ElConcatenero command later on.
for i in `ls -l ${MSADir} | tail -n +2 | awk '{print $9}'`
do
    ln -s ${MSADir}/${i}
    ConcInput=${i}" "${ConcInput}
done


# Finally, we will run the program creating three different outputs.

# One for use in MrBayes
${ElConcateneroDir}/PhD_Easy.py -if fasta -of nexus -o Pseudomonas_MSA_MrBayes -outgroup ${Outgroup} -in ${ConcInput}

# One for use in PhyML. This command also creates a second file with all the partitions.
${ElConcateneroDir}/PhD_Easy.py -if fasta -of phylip -o Pseudomonas_MSA_PhyML -outgroup ${Outgroup} -in ${ConcInput}

# One more phylip non interleaved format
${ElConcateneroDir}/PhD_Easy.py -if fasta -of mcmctree -o Pseudomonas_MSA_McMc -outgroup ${Outgroup} -in ${ConcInput}



echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0


