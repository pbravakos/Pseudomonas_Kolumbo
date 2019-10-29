#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
# #SBATCH --mem=128000
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Modeller_model"
#SBATCH --output=Modeller_model_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

## Some job specific info
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

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a


<< ////
       
	This script runs in the Modeller directory. Check the tutorial https://salilab.org/modeller/tutorial/basic.html

	To run this script, first CHANGE the Initial Parameters according to the analysis,and then run the script: 
	sbatch $0
	
	In this second part of the pipeline, we actually create the model for our protein, based on the selected PDB file (for the selection criteria check the first script file of this pipeline).


	NOTE:
	We run the python scripts in python3, because we have set up Modeller in our system (PYTHONPΑΤΗ) to run only in python3, but they could also run in python2.

////


# INITIAL PARAMETERS
StrainX=Strain23
Gene=3096
PDB=1flg
Chain=A
StrX=${StrainX/ain/}
PDBChain=${PDB}${Chain}


anviGetGeneCalls="anvi-get-sequences-for-gene-calls"

echo
echo "			ANLYSIS BEGINS!!!!"
echo

cd ${StrX}_${Gene}


# Next, we create the build_profile.py python file.
cat > build_profile.py <<EOF
from modeller import *

# First, we create an 'environ' object to use as input to later commands
env = environ()
# We create an empty alignment 'aln'
aln = alignment(env)
# Then we create a new protein model 'mdl', into which we read the chosen chain segment of the PDB structure file
mdl = model(env, file="${PDB}", model_segment=("FIRST:${Chain}","LAST:${Chain}"))
# Next, we transfer the PDB sequence of this model to the alignment and assign it a name (align_codes)
aln.append_model(mdl, align_codes="${PDBChain}", atom_files="${PDB}.pdb")
# Then we add the sequence from pir file to the alignment
aln.append(file="${StrX}_gene_${Gene}.pir", align_codes="${StrX}_${Gene}")
# The align2d() command is executed to align the two sequences.
aln.align2d()
# Finally, the alignment is written out in two formats, PIR and PAP. The PIR format is used by MODELLER in the subsequent model building stage, while the PAP alignment format is easier to inspect visually. In the PAP format, all identical positions are marked with a "*".
aln.write(file="${StrX}_${Gene}-${PDBChain}.pir", alignment_format='PIR')
aln.write(file="${StrX}_${Gene}-${PDBChain}.pap", alignment_format='PAP')
EOF

# Next, we run the build_profile.py python file.
python3 build_profile.py > build_profile.log


# Next, we create the model-single.py python file
cat > model-single.py <<EOF
from modeller import *
#  loads in the automodel class and prepares it for use
from modeller.automodel import *
#from modeller import soap_protein_od

env = environ()
# We then create an automodel object, call it 'a', and set parameters to guide the model building procedure. alnfile names the file that contains the target-template alignment in the PIR format. knowns defines the known template structure(s) in alnfile. sequence defines the name of the target sequence in alnfile. assess_methods requests one or more assessment scores. 
a = automodel(env, alnfile="${StrX}_${Gene}-${PDBChain}.pir",
              knowns='${PDBChain}', sequence="${StrX}_${Gene}",
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
# starting_model and ending_model define the number of models that are calculated (their indices will run from 1 to 5).
a.starting_model = 1
a.ending_model = 5
# The last line calls the make method that actually calculates the models.
a.make()
EOF


# Next, we actually run the model-single.py python file.
python3 model-single.py > model-single.log


<< ////

	MANUAL STEP NEEDED!
	We check the columns of the last table from the file "model-single.log" named as "Summary of successfully produced models". 
	Model selection criteria:
		1)  The model with the lowest value of the MODELLER objective function. Column "molpdf".
		2)  The model with the lowest value of the DOPE assessment score. Column "DOPE score".
		3)  The model with the highest GA341 assessment score. Column "GA341 score". 

	Before any external evaluation of the model, one should check the log file from the modeling run for runtime errors ("model-single.log") and restraint violations!    
	
////



echo
echo "==============================="
echo "SLURM job finished " "$(date)"
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0




