#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="PhyML"
#SBATCH --output=PhyML_job_%j.out

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

	This script runs from the PhyML master directory ($SLURM_SUBMIT_DIR). 
	In order to run it, the command should be like this:
	sbatch $0  

	IMPORTANT
	Don't forget BEFORE running this script to set the $NtaskMultiplier to any integer value (including zero) which will in turn determine the number of bootstraps, phyML is going to perform.
	
	IMPORTANT
	Model selection analysis by CodeML (from the PaML suite of tools) has preceded this analysis and we use the output AA rate matrices for this phylogeny. Optionally, we can also have a starting tree created by Mega, in newick format.

	IMPORTANΤ
	We will run PhyML with partitioned data and the --xml option. For this reason we must have prepared an XML file, according to the instructions found in the PhyML manual.
	
	IMPORTANT!
	If we want to run the parallel version of PhyML, we can only do it through MPI BUT the mpi version of PhyML does not support (yet, as of June 2019) the partitioned (through the xml option) bootstrap analysis.
	
        NOTE:
	The initial arrays, can be changed, according to the analysis. Hopefully, the result
	will still be a working xml file, ready for PhyML input.

	NOTE:
	We make use of the FreeRate model for the variation across sites. The initial values
	were copied from one example xml file from the PhyML manual, but we have chosen to 
	optimise these initial values (option "optimise.freerates"). For an example that uses
	the Gamma + I model (with 4 classes) instead of the FreeRate model, check the manual 
	and the relevant HERE document in the XML template.

END
}


# INITIAL VARIABLES
PhyMLStableDir=${HOME}/Software/phyml-3.3.20190321/src
PhyMLDevDir=${HOME}/Software/phyml/src
MSADir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Fasta2Phylip
XMLTemplate=Xml_Template_PhyML_Pseudomonas.sh

# Set an array with all the MSAs, we will use in this analysis.
declare -a Genes=("POG09090273" "POG090903FS"
		"POG090902FY" "POG0909018I"
		"POG090901ND" "POG090900GE"
		"POG090901XJ" "POG0909006H"
		"POG090901S2" "POG090902JS"
		)
# Next, we set an array with the Relative Rates, for 4 classes, for the FreeRate model (found in the xml file!).
declare -a RelativeRates=("0.197063" "0.750275"
                "1.951569" "5.161586"
                )
# Next, we set an array with the Frequencies (i.e. Weights), for 4 classes, for the FreeRate model (found in the xml file!).
declare -a Frequencies=("0.422481" "0.336848"
                "0.180132" "0.060539"               
               )

# Next, we set a new variable with the number of genes we are going to analyse.
NumGenes="${#Genes[@]}"

# Argument initialization for the loops in the xml template file.
ModelNum=1
Rates1=1
Rates2=1
NumGeneral=1
RatesNum1=1
RatesNum2=2
RatesNum3=3
RatesNum4=4

# Finally, we set the various levels of indentation for the xml template file (to have a nice xml output!)
level1="    "
level2=${level1}${level1}
level3=${level2}${level1}
level4=${level3}${level1}

#------------------------------BOOTSTRAP ASSIGNMENT---------------------------------------------------------------

# Set NtaskMultiplier to an integer data type. This variable will determine the number of bootstraps we will finally do.
declare -i NtaskMultiplier=0 # Can also be zero (0), for no bootstraps.

# The number of bootstrap replicates MUST be a multiple of the number of CPUs
BootStrap=$((${NtaskMultiplier} * $SLURM_NTASKS)) 

#-----------------------------------------------------------------------------------------------------------------

echo
# OpenMP is not slurm-aware. We have to define OMP_NUM_THREADS.
# Next variable setting, is for OpenMP when a job with --cpus-per-task is set. 
# Non parallel PhymL cannot run in OpenMP.
# We set OMP_NUM_THREADS to the same value as --cpus-per-task with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of --cpus-per-task, but only if --cpus-per-task is explicitly set.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}


echo
echo "			Start the PhyML analysis!!!!!!!"
echo
echo
echo "			PhyML is ready to run with ${BootStrap} bootstraps!"
echo


# Create soft links to all the MSAs and rate matrices used by the phyml xml file.
for i in "${Genes[@]}"
do
    PhylipMSA=${i}_MSA_relaxed.phy
    if [[ ! -e ${MSADir}/${PhylipMSA} ]]; then
    echo "File ${PhylipMSA} cannot be found in directory ${MSADir}. Please prepare a phylip (or nexus) MSA file and run this script again." >&2
    generalInfo >&2
    exit 1
    fi

    if [[ ! -e ${PhylipMSA} ]]; then
        ln -s ${MSADir}/${PhylipMSA}
    fi
    

    CodeMLDir=${HOME}/Pseudomonas/CodeML/${i}
    AArateMatrix=${i}_AAratefile.dat
    if [[ ! -e ${CodeMLDir}/${AArateMatrix} ]]; then
    echo "File ${AArateMatrix} cannot be found in directory ${CodeMLDir}. Please run CodeML to generate the AA rate matrix and then run this script again." >&2
    generalInfo >&2
    exit 1
    fi
    
    if [[ ! -e ${AArateMatrix} ]]; then
        ln -s ${CodeMLDir}/${AArateMatrix}
    fi
done

# Create the XML template.
cat > ${XMLTemplate} <<"EOF"
echo "<phyml output.file=\"Pseudomonas\" bootstrap=\"${BootStrap}\" branch.test=\"no\" quiet=\"yes\" memory.check=\"no\">"
echo
echo
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooo TOPOLOGY oooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- Tree topology: start with BioNJ and then SPRs -->"
echo "${level1}<topology> "
echo "${level2}$<instance id=\"T1\" init.tree=\"bionj\" optimise.tree=\"yes\" n.rand.starts=\"5\" search=\"spr\"/>"
echo "${level1}</topology>"
echo
echo
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooo SUBSTITUTION MODELS ooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo
echo "${level1}<!-- If optimise.weights=yes, then the probabilities (or weights) or each matrix in the set of matrices defined by this component, will be estimated from the data -->"
echo "${level1}<ratematrices id=\"RM1\" optimise.weights=\"yes\">"

for i in "${Genes[@]}"
do
    Model="M""$ModelNum"
    echo "${level2}<instance id=\"${Model}\" model=\"customaa\" ratematrix.file=\"${i}_AAratefile.dat\" optimise.rr=\"yes\"/>"
    ModelNum=$(($ModelNum + 1))   
done
echo "${level1}</ratematrices>"
echo
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooo VARIATION OF RATES ACROSS SITES ooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo 

for repeat in $( seq 1 ${NumGenes} )
do
    echo "${level1}<siterates id=\"SR${repeat}\">"
    echo "${level2}<!-- Freerate model of variation of rates across sites FreeRate is slightly more computationally demanding than the Γ (or Γ+I) model, it often provides a significantly better fit to the data -->" 
    echo "${level2}<!-- First, we define the relative rates --> "
    for i in "${RelativeRates[@]}"
    do
        Rate="R""$Rates1"
        echo "${level2}<instance id=\"${Rate}\" init.value=\"${i}\"/>"
        Rates1=$(($Rates1 + 1))
    done
    
    echo "${level3}<!-- 'optimise.freerates' option optimises the parameters of the FreeRate model, i.e., the relative rates and the corresponding frequencies -->"
    echo "${level3}<weights  id=\"D1\" family=\"freerates\" optimise.freerates=\"yes\">"
    echo "${level4}<!-- Next, we specify the frequency, or weight, of each class of rate -->"

    for i in "${Frequencies[@]}"
    do
        Rate="R""$Rates2"
        echo "${level4}<instance appliesto=\"${Rate}\" value=\"${i}\"/>"
        Rates2=$(($Rates2 + 1)) 
    done
    echo "${level3}</weights>"
    echo "${level1}</siterates>"
    echo
done

# Bellow is an example of how a Gamma + Invariant model can be applied:
<<END

<!-- It is possible to use two distinct models. Next we use the Γ4+I model -->
<!-- Note that one of the initial (relative) rate (init.value attribute) is set to zero. The corresponding rate class (the third in this example) will then correspond to the invariant site category. -->
  <siterates id="SRGTR">
    <weights  id="D2" family="gamma+inv" alpha=".1" optimise.alpha="yes" pinv="0.4" optimise.pinv="yes">
    </weights>
    <instance id="R1" init.value="1.0"/>
    <instance id="R2" init.value="1.0"/>
    <instance id="R3" init.value="0.0"/>
    <instance id="R4" init.value="1.0"/>
    <instance id="R5" init.value="1.0"/>
  </siterates>

END

echo  
echo 
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooo AMINO ACID FREQUENCIES oooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo
echo "${level1}<!-- If optimise.weights=yes, then the probabilities (or weights) or each vector of equilibrium frequencies in the set of vectors defined by this component, will be estimated from the data. -->"
echo "${level1}<equfreqs id=\"EF1\" optimise.weights=\"yes\">"
for repeat in $( seq 1 ${NumGenes} )
do
    echo "${level2}<instance id=\"F${repeat}\" aa.freqs=\"model\"/>"
done
echo  "${level1}</equfreqs>"
echo
echo
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- oooooooooooooooooo EDGE LENGTHS ooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo
echo  "${level1}<branchlengths id=\"BL1\" >"
for repeat in $( seq 1 ${NumGenes} )
do
    echo "${level2}<instance id=\"L${repeat}\" optimise.lens=\"yes\"/>"
done
echo "${level1}</branchlengths>"
echo
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo '<!-- ooooooooooo SET UP OF THE PARTITIONNED ANALYSIS ooooooooooooo -->'
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo
echo "${level1}<!-- Mixture model assemblage -->"
echo "${level1}<!-- The ordering of the different term in these list matters a lot since it is directly related to the elements in each class of the mixture model. -->"
echo "${level1}<!-- Branch lengths will be estimated separately for each partition element -->"
echo "${level1}<!-- Note that a given partition element can only have one branchlengths instance associated to it -->"

for Gene in "${Genes[@]}"
do
    echo "${level1}<partitionelem id=\"partition${NumGeneral}\" file.name=\"${Gene}_MSA_relaxed.phy\" data.type=\"aa\" interleaved=\"no\">"
    echo "${level2}<mixtureelem list=\"T1, T1, T1, T1\"/>"
    echo "${level2}<mixtureelem list=\"M${NumGeneral}, M${NumGeneral}, M${NumGeneral}, M${NumGeneral}\"/>"
    echo "${level2}<mixtureelem list=\"F${NumGeneral}, F${NumGeneral}, F${NumGeneral}, F${NumGeneral}\"/>"
    echo "${level2}<mixtureelem list=\"R${RatesNum1}, R${RatesNum2}, R${RatesNum3}, R${RatesNum4}\"/>"
    echo "${level2}<mixtureelem list=\"L${NumGeneral}, L${NumGeneral}, L${NumGeneral}, L${NumGeneral}\"/>"
    echo "${level1}</partitionelem>"
    echo
    NumGeneral=$(($NumGeneral + 1))
    RatesNum1=$(($RatesNum1 + 4))
    RatesNum2=$(($RatesNum2 + 4))
    RatesNum3=$(($RatesNum3 + 4))
    RatesNum4=$(($RatesNum4 + 4))
done
  
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo "<!-- ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo -->"
echo
echo "</phyml>"
EOF

# Create the actual xml file that will be used as input for PhyML.
source ${XMLTemplate} > PhyML_Pseudomonas.xml


# Run the non parallel version of PhyML. Bootstraps may also run in the non parallel version of PhyML, but it will take ages!
${PhyMLDevDir}/phyml --xml=PhyML_Pseudomonas.xml

# IMPORTANT!!!!!!!!!!!!
# Current verstion of PhyML DOES NOT support mpi for xml input!
# MPI is 'slurm-aware'. We do not need to specify the -np nor the -host, or hostfile options to mpirun or mpiexec
# Run the parallel version of PhyML
#mpirun -np $SLURM_NTASKS  ${PhyMLDevDir}/phyml-mpi --xml=PhyML_Pseudomonas.xml


# Or, alternatively run the same job with srun.
# srun --mpi=openmpi ${PhyMLDevDir}/phyml-mpi --xml=PhyML_Pseudomonas.xml
# --ntasks $SLURM_NTASKS



echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
