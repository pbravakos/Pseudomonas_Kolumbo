#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
# #SBATCH --mem=128000
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Modeller_eval"
#SBATCH --output=Modeller_eval_job_%j.out

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

	To run this script, first CHANGE the Initial Parameters according to needs,and then run the script: 
	sbatch $0
	
	In this third part of the pipeline, we evaluate the model created previously.


	NOTE:
	We run the python scripts in python3, because we have set up Modeller in our system (PYTHONPΑΤΗ) to run only in python3, but they could also run in python2.

////


# INITIAL PARAMETERS
StrainX=Strain23
Gene=3096
PDB=1flg
Chain=A
StrX=${StrainX/ain/}
PDBSelected=${StrX}_${Gene}.B99990003.pdb
PDBChain=${PDB}${Chain}

anviGetGeneCalls="anvi-get-sequences-for-gene-calls"

echo
echo "			ANLYSIS BEGINS!!!!"
echo

cd ${StrX}_${Gene}


# Next, we create the evaluate_model.py python file.
cat > evaluate_model.py <<"EOF"
from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters
EOF


cat >> evaluate_model.py <<EOF
# read model file
# we use the complete_pdb script to read in a PDB file and prepare it for energy calculations (this automatically allows for the possibility that the PDB file has atoms in a non-standard order, or has different subsets of atoms, such as all atoms including hydrogens, while MODELLER uses only heavy atoms, or vice versa).
mdl = complete_pdb(env, '${PDBSelected}')

# We then create a selection of all atoms, since most MODELLER energy functions can operate on a subset of model atoms.
# Assess with DOPE:
s = selection(mdl)   # all atom selection
# we additionally request an energy profile, smoothed over a 15 residue window, and normalized by the number of restraints acting on each residue. This profile is written to a file, which can be used as input to a graphing program
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file="${StrX}_${Gene}.profile",
              normalize_profile=True, smoothing_window=15)
EOF

# Next, we run the evaluate_model.py python file.
python3 evaluate_model.py > evaluate_model.log



# Next, we create the evaluate_template.py python file.
cat > evaluate_template.py <<"EOF"
from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
EOF


cat >> evaluate_template.py <<EOF
# read model file
mdl = complete_pdb(env, "${PDB}.pdb", model_segment=('FIRST:A', 'LAST:A'))

s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='${PDBChain}.profile',
              normalize_profile=True, smoothing_window=15)
EOF


# Next, we run the evaluate_template.py python file.
python3 evaluate_template.py > evaluate_template.log

# IMPORTANT!!!
# For the plot_profiles.py python script to work, we need to have X11 forwarding (i.e. ssh -X). Otherwise we get an error!

cat > plot_profiles.py <<EOF
import pylab
import modeller

def r_enumerate(seq):
    """Enumerate a sequence in reverse order"""
    # Note that we don't use reversed() since Python 2.3 doesn't have it
    num = len(seq) - 1
    while num >= 0:
        yield num, seq[num]
        num -= 1

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in r_enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

e = modeller.environ()
a = modeller.alignment(e, file="${StrX}_${Gene}-${PDBChain}.pir")

template = get_profile("${PDBChain}.profile", a["${PDBChain}"])
model = get_profile("${StrX}_${Gene}.profile", a["${StrX}_${Gene}"])

# Plot the template and model profiles in the same plot for comparison:
pylab.figure(1, figsize=(10,6))
pylab.xlabel('Alignment position')
pylab.ylabel('DOPE per-residue score')
pylab.plot(model, color='red', linewidth=2, label='Model')
pylab.plot(template, color='green', linewidth=2, label='Template')
pylab.legend()
pylab.savefig('dope_profile.png', dpi=300)
EOF

# Next, we COULD run the plot_profiles.py python script, but we choose not to run it. We will run it later when we have X11 forwarding.
# python3 plot_profiles.py > plot_profiles.log



<< ////

	MANUAL STEP NEEDED!
	Connect with X11 forwarding to the server and run the following command in the working directory ${StrX}_${Gene}
	python3 plot_profiles.py > plot_profiles.log
	
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




