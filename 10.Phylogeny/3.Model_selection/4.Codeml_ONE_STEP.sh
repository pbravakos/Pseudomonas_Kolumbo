#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="CodeML"
#SBATCH --output=CodeML_job_%j.out


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

	This script runs in the CodeML directory found in $SLURM_SUBMIT_DIR.
	In order to run it, it needs one argument; the gene name we want to analyze e.g.:
	sbatch $0 POG090901ND

	It also runs as part of another script, which calls this one, each time with a selected gene name.

	NOTE:
	The control file was created based on the instructions of the PAML FAQ, the section titled "How can I estimate an amino acid substitution matrix from my own data, like mtmam.dat and wag.dat?". Here we follow the one step aproach.
	
	NOTE:
	The output of codeml is an amino acid Rate matrix, which will be used as input in the PhyML program.
	
END
}


# INITIAL PARAMETERS
Gene=$1
CodeMLBatchDir=${HOME}/Software/paml4.9i/automate-PAML-codeml
MSADir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Fasta2Phylip
MegaDir=${HOME}/Pseudomonas/MegaCC
MSAPhylip=${Gene}_MSA_relaxed.phy
InitialTree=${Gene}_ML_MegaCC.nwk
Output=${Gene}_PAML_One_Step.out



# Start a job dependency to move the Sdtout file to the correct folders
"${HOME}"/sbatch --dependency=afterany:"$SLURM_JOB_ID" Move_Slurm_output_files.sh "$SLURM_JOB_NAME" "$SLURM_JOB_ID" "$1"


# Check whether there is the directory with the gene name.
if [[ ! -d ${Gene} ]]; then
    mkdir ${Gene}
else
    rm -rf ${Gene}
    mkdir ${Gene}
fi


# Change to the working directory
cd ${Gene}


# Create the codeml.ctl file.
cat > codeml.ctl <<EOF
* Note that spaces are required on both sides of the equal sign, and blank lines or lines beginning with "*" are treated as comments.
* Options not used can be deleted from the control file.
* The order of the variables is unimportant.
seqfile = ${MSAPhylip} * sequence data file name
treefile = ${InitialTree} * tree structure file name
outfile = ${Output} * main result file name
noisy = 3 * 0,1,2,3,9: how much rubbish on the screen
verbose = 1 * 1: detailed output, 0: concise output
aaRatefile = ${HOME}/Software/paml4.9i/dat/wag.dat * wag.dat for nuclear proteins
model = 9
fix_alpha = 0
alpha = 0.5
EOF

# Create links to the files referenced in the codeml.ctl file.
ln -s ${MSADir}/${MSAPhylip}
ln -s ${MegaDir}/${InitialTree}


echo
# CodeML should be in the path!
echo "CodeML runs from this path:"
which codeml
echo

# run the program! 
codeml

# Not necessary to sleep, but to be on the safe side...
sleep 1

# AA substitution rate matrix is printed into a file named "AAratefile.dat"
# Finally, we prepare the matirx for PhyML!
head -n -4 AAratefile.dat | sed -r 's/ +/ /g;s/^ //g;1d' > ${Gene}_AAratefile.dat


echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
