#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Prank_PhyloGen"
#SBATCH --output=Prank_PhyloGen_job_%j.out


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

	This script runs in the Prank directory as part of another script which calls this one with a simple for loop and assigns one gene name as the first argument of this script. The fasta file with the genes was created with Sed in the previous steps of this pipeline.
       
	In order to run it standallone the command should be like this:
	./$0 POG090903HI
	
	For the list of genes, please check the previous steps of this pipeline.

END
}


# INITIAL PARAMETERS
Gene=$1
SedFormatedFastaDir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Genes_fasta_headers_and_gaps_stripped

# Check that the working directory exists.
[[ ! -d $SedFormatedFastaDir ]] && { echo "Directory ${SedFormatedFastaDir} cannot be not found!" >&2; exit 1; }

# Check that the fasta file exists!
if [[ ! -e ${SedFormatedFastaDir}/Pseudomonas_${Gene}.fas ]]; then
    echo "The file Pseudomonas_${Gene}.fas cannot be found in ${SedFormatedFastaDir}. Please run the previous steps of this pipeline first!" >&2
    generalInfo >&2
    exit 1
fi

# Now run the command!
prank -d=${SedFormatedFastaDir}/Pseudomonas_${Gene}.fas \
	-o=Pseudomonas_${Gene}_prank.fasta \
	-protein \
	-iterate=200 \
	-f='fasta' \
	+F \
	-verbose


echo "======================================================"
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
