#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="RaXML"
#SBATCH --output=RaXML_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

# Some job specific info
echo "Job ID is = " "$SLURM_JOBID"
echo "SLURM cluster name = " "$SLURM_CLUSTER_NAME"
echo "SLURM partition = " "$SLURM_JOB_PARTITION"
echo "SLURM node list = " "$SLURM_JOB_NODELIST"
echo "SLURM num of nodes = " "$SLURM_JOB_NUM_NODES"
echo "SLURM number of tasks = " "$SLURM_NTASKS"
echo "SLURM cpus per task" = " "${SLURM_CPUS_PER_TASK}
echo "SLURM memory per node = " "$SLURM_MEM_PER_NODE"
echo "SLURM memory per cpu = " "$SLURM_MEM_PER_CPU"
echo "working directory = " "$SLURM_SUBMIT_DIR"
echo "=================================================="
echo "SBATCΗ job started " "$(date)"
echo "=================================================="
echo

generalInfo () {
    cat <<END

	This script runs from the RaXML master directory ($SLURM_SUBMIT_DIR). 
	In order to run it, the command should be like this:
	sbatch $0  

	IMPORTANT:
	We have concatenated all MSA with the software ElConcetenero, creating a phylip formated MSA and a partition data file. 

	IMPORTANT:
	We have selected the best model for each partition with Protest.

	IMPORTANT:
	We will use as a starting tree the tree from the PhyML run.

	NOTE ON EFFECTIVE USAGE OF MPI:
	The best and easiest way to find the recommended number of threads (and memory) for the data is to run RaXML with the option --parse before running the actual phylogeny.
	MPI and pthreads do exactly the same thing: Assign threads to analyse the MSA. There is no MPI parallelization for bootstraps. Bootstraps happen in a serial way!
	With fine-grained parallelization, the number of CPU cores that can be efficiently utilized is limited by the alignment 'width' (=number of sites). 
	For instance, using 20 cores on a single-gene protein alignment with 300 sites would be suboptimal, and using 100 cores would most probably result in a huge slowdown. 
	In order to prevent the waste of CPU time and energy, RAxML-NG will warn you -- or, 
	in extreme cases, even refuse to run -- if you try to assign too few alignment sites per core.
	In MPI-only mode, you should start 1 MPI process per physical CPU core (the number of threads will be set to 1 by default).
	Check for more here: https://github.com/amkozlov/raxml-ng/wiki/Parallelization

	NOTE ON MEMORY USAGE:
	
	Here is a formula to estimate RAxML memory requirements: 
	Given an alignment of n taxa and m distinct patterns the memory consumption is approximately:
	MEM(AA+GAMMA) = (n-2) * m * (80 * 8) bytes
	MEM(AA+CAT) = (n-2) * m * (20 * 8) bytes
	MEM(DNA+GAMMA) = (n-2) * m * (16 * 8) bytes
	MEM(DNA+CAT) = (n-2) * m * (4 * 8) bytes
	
END
}




# INITIAL VARIABLES
RaXMLMpiDir=/mnt/big/Tools
MSADir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Concatenate_Alignments/ElConcetenero
PhyMLDir=${HOME}/Pseudomonas/PhyML
PhyMLNwkTree=Pseudomonas_phyml_tree.txt
PhylipMSA=Pseudomonas_MSA_PhyML.phy
ParitionFile=Pseudomonas.part

# RaXML parameters that can be changed
outgroup=Cellvibrio_japonicus_Ueda107
BsCutoff=0.01
Prefix=Pseudomonas_RaxML
MaxNumBoot=1000


if [[ ! -e ${MSADir}/Pseudomonas_MSA_PhyML.phy ]]; then
    echo "MSA Pseudomonas_MSA_PhyML.phy cannot be found in directory ${MSADir}. Please upload it and run this script again."
    generalInfo
    exit 1
fi



echo
echo "				Analysis begins!!!!!!!!!!!!!!!!!!!!!!"
echo
echo

# Create a link to the MSA and the starting tree.
ln -s ${MSADir}/${PhylipMSA}
ln -s ${PhyMLDir}/${PhyMLNwkTree}

# Create the partion data file. Information is derived from both ElConcetenero (for the partitions)  and Protest (for the actual models).
cat > ${ParitionFile} <<"EOF"
LG+IO+G+FO, POG090903FS = 1-365
LG+IO+G+FO, POG090902JS = 366-1226
LG+IO+G, POG090902FY = 1227-1570
Dayhoff+IO+G+FO, POG09090273 = 1571-1999
Dayhoff+IO+G+FO, POG090901XJ = 2000-2209
JTT+IO+G+FO, POG090901S2 = 2210-2373
LG+IO+G+FO, POG090901ND = 2374-2769
LG+IO+G, POG0909018I = 2770-3200
LG+IO+G+FO, POG090900GE = 3201-3427
Dayhoff+IO+G+FO, POG0909006H = 3428-3571
EOF


# Optionally, but HIGHLY recommended is to run the following command before the main command.
# ${RaXMLMpiDir}/raxml-ng-mpi --parse --msa ${PhylipMSA} --model ${ParitionFile} --prefix ${PrefixParser}
# This command will generate a binary MSA file (${PrefixParser}.raxml.rba), which can be loaded by RAxML-NG much faster than the original FASTA alignment. 
# Furthermore, it will print estimated memory requirements and the recommended number of threads for this dataset

# mpirun -np $SLURM_NTASKS 
${RaXMLMpiDir}/raxml-ng-mpi_2 --all \
			--outgroup ${outgroup} \
			--msa ${PhylipMSA} \
			--brlen scaled \
			--model ${ParitionFile} \
			--bs-trees autoMRE{${MaxNumBoot}} \
			--bs-cutoff ${BsCutoff} \
			--msa-format PHYLIP \
			--data-type AA \
			--prefix ${Prefix} \
			--threads ${SLURM_CPUS_PER_TASK} \
			--tree ${PhyMLNwkTree} \
			--extra thread-pin


echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
