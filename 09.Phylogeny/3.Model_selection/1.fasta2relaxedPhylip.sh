#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name="fasta2relaxedPhylip"
#SBATCH --output=fasta2relaxedPhylip_job_%j.out


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

	This script runs in the Anvio Phylogenomics found in $SLURM_SUBMIT_DIR.
	In order to run it, we should give one argument, the gene of interest e.g.:
	sbatch $0  POG090900GE # for slurm
	or
	bash $0 
	
	NOTE:
	Here we want to create an MSA to be used as input in CodeML (part of PAML). For this reason we use the MSA output from Guidance.
 
	NOTE:
	We have cloned the Phylogenomic repo from github https://github.com/npchar/Phylogenomic. 

END
}


# INITIAL PARAMETERS
MSADir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Final_Selected_genes
PhylogenomicDir=${HOME}/Software/Phylogenomic

if [[ ! -d ${MSADir} ]]; then
    echo "Directory ${MSADir} cannot be found!"
    generalInfo
    exit 1
fi



declare -a arr=("POG09090273" "POG090903FS"
		"POG090902FY" "POG0909018I"
		"POG090901ND" "POG090900GE"
		"POG090901XJ" "POG0909006H"
		"POG090903B5" "POG090901S2"
		"POG090902JS"
		)


for i in "${arr[@]}"
do
    
    GuidanceMSA=${i}_Guidance_ColFiltered.fas
    if [[ ! -e ${MSADir}/${GuidanceMSA} ]]; then
        echo "${GuidanceMSA} cannot be found in ${MSADir}. Please upload it and run again the script!" >&2
    else
        ln -s ${MSADir}/${GuidanceMSA}
        # We will use the fasta2relaxedPhylip perl script to turn the fasta MSA to a relaxed, sequential phylip format.
        ${PhylogenomicDir}/fasta2relaxedPhylip.pl -f ${GuidanceMSA} -o ${i}_MSA_relaxed.phy > ${i}.log
        # The output from fasta2relaxedPhylip is separating headers from sequences by a tab character. According to the paml manual, we have to have at least two spaces as separators between these fields.
        sed -ir 's/\t/   /g' ${i}_MSA_relaxed.phy
        sleep 1
    fi

done

# We will use the fasta2relaxedPhylip perl script to turn the fasta MSA to a relaxed, sequential phylip format.
${PhylogenomicDir}/fasta2relaxedPhylip.pl -f ${GuidanceMSA} -o ${Gene}_MSA_relaxed.phy > ${Gene}.log

# The output from fasta2relaxedPhylip is separating headers from sequences by a tab character. According to the paml manual, we have to have at least two spaces as separators between these fields.
# PAML considers two consecutive spaces as the end of a species name, so that the species name does not have to have exactly 30 (or 10) characters. To make this rule work, you should not have two consecutive spaces within a species name (From PAML manual).
sed -ir 's/\t/   /g' ${Gene}_MSA_relaxed.phy




echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
