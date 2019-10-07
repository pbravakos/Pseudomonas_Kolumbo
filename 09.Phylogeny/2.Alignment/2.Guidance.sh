#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Guidance"
#SBATCH --output=Guidance_job_%j.out


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

	This script runs in the Guidance master directory. In order to run it the command should be like this:
	./$0 POG090903HI    # POG090903HI is the gene of interest

	# IMPORTANT
	In the guidance command bellow we have as input an MSA and we want to replicate the MSA command that produced this particular MSA. In our case that was prank with the +F option.

END
}


# INITIAL PARAMETERS
Gene=$1
SeqType=AminoAcid
MSADir=${HOME}/Pseudomonas/Prank/Phylogenomics/${SeqType}
SeqDir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Genes_fasta_headers_and_gaps_stripped
GuidanceDir=${HOME}/Software/guidance.v2.02/www/Guidance
OutDir=${HOME}/Pseudomonas/Guidance/Phylogenomics/${SeqType}/Results_MSA_${Gene}
MSAfile=Pseudomonas_${Gene}_prank.fasta.best.fas
Seqfile=Pseudomonas_${Gene}.fas

# Guidance parameters that can be changed
SeqCutOFF=0.6   # Please, set this to a small number to prevent removing sequences!
ColCutOFF=0.93   # This is the important number, the lower it is set the more columns will be removed from the final MSA.
BootStraps=200

mkdir ${OutDir}

# Now run the command!
LC_ALL=C perl ${GuidanceDir}/guidance.pl --proc_num ${SLURM_NTASKS} \
					--seqFile ${SeqDir}/${Seqfile} \
					--msaFile ${MSADir}/${MSAfile} \
					--msaProgram PRANK \
					--MSA_Param \\-F \
					--outDir ${OutDir} \
					--dataset ${Gene} \
					--seqCutoff ${SeqCutOFF} \
					--colCutoff ${ColCutOFF} \
					--seqType aa \
					--bootstraps ${BootStraps}

# Rename the final output that will be used for downstream analyses!
cp ${OutDir}/${Gene}.PRANK.Without_low_SP_Col.With_Names ${OutDir}/${Gene}_Guidance_ColFiltered.fas


echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
