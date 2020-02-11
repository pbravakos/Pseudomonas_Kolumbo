#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Anvio_ANI"
#SBATCH --output=Anvio_ANI_job_%j.out

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

	This script calculates the Average Nucleotide Identity between all the internal genomes we intend to run for pangenomics. It runs in the Internal Prokka Anvio pangenome directory which can be found here: ${HOME}/Pseudomonas/Anvio/Pangenome/Internal_genomes/Prokka/ANI
        
        In order to run it standallone the command should be like this:
        sbatch $0

	IMPORTANT:
	We need to create the internal tsv file first and the Pangenome database in Anvio, before running this command.
	
END
}


# INITIAL PARAMETERS
AnvioDir=${HOME}/Pseudomonas/Anvio/Pangenome/Internal_genomes/Prokka
AnvioPangenomeDir=${AnvioDir}/PSEUDOMONAS_PROKKA
AnvioAniOutput=${AnvioDir}/ANI/Output



if [[ ! -e ${AnvioDir}/Internal_genomes.tsv ]]; then
    echo "File 'Internal_genomes.tsv' cannot be found on ${AnvioDir}. Please first create the file and then run this script again"
    exit 1
fi


# anvi-compute-ani --internal-genomes ${AnvioDir}/Internal_genomes.tsv --pan-db ${AnvioPangenomeDir}/Pseudomonas_Prokka_Pangenome-PAN.db  --output-dir ${AnvioAniOutput} --method ANIb --num-threads $SLURM_NTASKS --just-do-it --log-file anvi-ANI.log

anvi-compute-ani --internal-genomes ${AnvioDir}/Internal_genomes.tsv --pan-db ${AnvioPangenomeDir}/Pseudomonas_Prokka_Pangenome-PAN.db  --output-dir ${AnvioAniOutput} --method ANIm --num-threads $SLURM_NTASKS --just-do-it --log-file anvi-ANI.log


echo "==========================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
