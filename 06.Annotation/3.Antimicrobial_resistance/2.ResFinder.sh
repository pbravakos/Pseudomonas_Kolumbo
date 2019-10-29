#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="ResFinder"
#SBATCH --output=Sdout_ResFinder_job_%j.txt

# for calculating the amount of time the job takes
begin=`date +%s`

# Some job specific info
echo "Job ID is = " $SLURM_JOBID
echo "SLURM cluster name = " $SLURM_CLUSTER_NAME
echo "SLURM partition = " $SLURM_JOB_PARTITION
echo "SLURM node list = " $SLURM_JOB_NODELIST
echo "SLURM num of nodes = " $SLURM_JOB_NUM_NODES
echo "SLURM number of tasks = " $SLURM_NTASKS
echo "SLURM memory per node = " $SLURM_MEM_PER_NODE
echo "SLURM memory per cpu = " $SLURM_MEM_PER_CPU
echo "working directory = " $SLURM_SUBMIT_DIR
echo "=================================================="
echo "SBATCΗ job started " `date`
echo "============================================================================"

generalInfo () {
cat <<END

	This is the ResFinder pipeline

	This script takes as input one argument. 
	For Strain 01 that would be: 
	$0 Strain01

	We run the program twice:
	One time on the nucleotide fna output from Anvio
	A second time on the nucleotide fasta contigs
	These results should be exactly the same, with the only difference that from the Anvio fna output, we can relate the gene callers id to the annotation we have done for Anvio.
	We are not going to parse these results into Anvio, because not all Strains have antimicrobial resistance genes, so this could pose problems to our Anvio pipeline (especially for pangenomics, where we need to have done exactly the same analyses for all of our strains).
	 
	ResFinder was downloaded from https://bitbucket.org/genomicepidemiology/resfinder.

END
}



# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}

ResFinderDir=${HOME}/Software/resfinder
ResFinderDBDir=${ResFinderDir}/resfinder_db
ContigDir=${HOME}/Pseudomonas/Seqtk/${StrainX}
AnvioProkkaDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}

ContigFnaFile=${StrX}_Pilon_CLA_Blast.fasta
AnvioFnaFile=${StrX}_gene_calls.fna

# Parameters for ResFinder that can be changed
ContigOutput=ContigSearch
AnvioFnaOutput=AnvioFnaSearch


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] &&  [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Error: Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi

mkdir -p ${StrainX}

# Start a job dependency to move the Sdtout file to the correct folders
sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_ID $1




#----------------------------------------------------------------------------------------------------------

cd ${StrainX}

#--------------------------SANITY CHECKS -------------------------------------------------------------------#

if [[ ! -e ${ResFinderDir}/resfinder.py ]]; then
    echo "Resfinder executable is missing! Please check that it can be found here: ${ResFinderDir}/resfinder.py" >&2
    exit 1
fi

if [[ ! -s ${AnvioProkkaDir}/${AnvioFnaFile} ]]; then
    echo "The anvio fna file is missing. Please check for ${AnvioProkkaDir}/${AnvioFnaFile}" >&2
    exit 1
fi

if [[ ! -s ${ContigDir}/${ContigFnaFile} ]]; then
    echo " The Contig file is missing. Please check for ${ContigDir}/${ContigFnaFile}" >&2
    exit 1
fi

# Please comment the following two (2) lines of code the second time you run this program (in order to save some time and still have the latest database!)
#rm -rf ${ResFinderDBDir}
#git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git ${ResFinderDBDir}

## We want to always, download the latest, database. So please, remove the previous database and install a fresh one, before running this code. Cannot do it automatically, because if multiple jobs run the same script at the same time then problems start to arise. This is why the else statement bellow must be commented, if multiple jobs are running at the same time! 
#if [[ ! -d ${ResFinderDBDir} ]]; then
#   git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git ${ResFinderDBDir}
##else
##    rm -rf ${ResFinderDBDir}
##    git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git ${ResFinderDBDir}
#fi
    

echo
echo
echo "			ResFinder Analysis Begins!"
echo

mkdir -p ${ContigOutput} ${AnvioFnaOutput}

# For the next command to run we have to have blastn on our path! Otherwise declare the full path of the blastn executable with the -b option.
${ResFinderDir}/resfinder.py -i ${AnvioProkkaDir}/${AnvioFnaFile} -o ${AnvioFnaOutput} -p ${ResFinderDBDir} 

${ResFinderDir}/resfinder.py  -i ${ContigDir}/${ContigFnaFile} -o ${ContigOutput} -p ${ResFinderDBDir}


# Next, rename the output, to something more meaningful to us.
mv ${AnvioFnaOutput}/results_tab.txt ${AnvioFnaOutput}/${StrX}_Anvio_results.tsv
mv ${ContigOutput}/results_tab.txt ${ContigOutput}/${StrX}_Contigs_results.tsv

mv ${AnvioFnaOutput}/results.txt ${AnvioFnaOutput}/${StrX}_Anvio_results.txt
mv ${ContigOutput}/results.txt ${ContigOutput}/${StrX}_Contigs_results.txt

mv ${AnvioFnaOutput}/results_table.txt ${AnvioFnaOutput}/${StrX}_Anvio_results_table.tsv
mv ${ContigOutput}/results_table.txt ${ContigOutput}/${StrX}_Contigs_results_table.tsv

mv ${AnvioFnaOutput}/Resistance_gene_seq.fsa ${AnvioFnaOutput}/${StrX}_Anvio_resistance_gene.fna
mv ${ContigOutput}/Resistance_gene_seq.fsa ${ContigOutput}/${StrX}_Contigs_resistance_gene.fna

mv ${AnvioFnaOutput}/Hit_in_genome_seq.fsa ${AnvioFnaOutput}/${StrX}_Anvio_hit_in_genome.fna
mv ${ContigOutput}/Hit_in_genome_seq.fsa ${ContigOutput}/${StrX}_Contigs_hit_in_genome.fna

# Finally, remove temporary files
rm -r ${ContigOutput}/tmp ${AnvioFnaOutput}/tmp

echo
echo "-----------------------------RESFINDER ANALYSIS COMPLETE!!!------------------------------------------------------------------------"
echo 


echo "========================================================================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
