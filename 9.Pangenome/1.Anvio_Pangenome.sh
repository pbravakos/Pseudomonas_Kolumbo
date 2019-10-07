#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
# #SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Pangenome_Internal"
#SBATCH --output=Pangenome_Internal_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

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
echo "=================================================="
echo

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a


<< ////

        This is the Anvio pangenomics pipeline. Check the tutorials:
        http://merenlab.org/2016/11/08/pangenomics-v2/ "An anvi'o workflow for microbial pangenomics" and
        http://merenlab.org/2016/06/22/anvio-tutorial-v2/ "Anvi'o User Tutorial for Metagenomic Workflow"
	http://merenlab.org/tutorials/infant-gut/#chapter-v-pangenomics 

	This script run in the Anvio Pangenomics Prokka directory. To run the script, we need to have all the needed files and dependencies (blastp, mcl) and then we simply run the command:
	./script.sh

	IMPORTANT
        We have already created the contigs DB from the Prokka output.

	NOTE:
	Useful command for pangenomes NOT used here is:  anvi-get-sequences-for-gene-clusters  

////

# INITIAL PARAMETERS
OutputGenomesDir=PSEUDOMONAS_PROKKA

DescriptionTxt=Pangenome_Prokka_description.txt
GenomesDatabase=PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db
StrainsInternal=StrainsUsed.txt

# Parameters that could be changed
MINBIT=0.5
MCLInflation=7


<< ////

	A NOTE ON MINBIT
	Use the minbit heuristic that was originally implemented in ITEP (Benedict et al, 2014) to eliminate weak matches between two amino acid sequences. 
	You see, the pangenomic workflow first identifies amino acid seqeunces that are somewhat similar by doing similarity searches, and then resolves gene clusters based on those similarities. 
	In this scenario, weak similarities can connect gene clusters that should not be connected. Although the network partitioning algorithm can recover from these weak connections, 
	it is always better to eliminate as much noise as possible at every step. So the minbit heuristic provides a mean to set a to eliminate weak matches between two amino acid sequences. 
	We learned it from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8), which is another comprehensive analysis workflow for pangenomes, and decided to use it because it makes a lot of sense. 
	Briefly, If you have two amino acid sequences, A and B, the minbit is defined as 
	BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B)). 
	So the minbit score between two sequences goes to 1.0 if they are very similar over the entire length of the ‘shorter’ amino acid sequence, 
	and goes to 0.0 if (1) they match over a very short stretch compared even to the length of the shorter amino acid sequence or (2) the match between sequence identity is low. 
	The default minbit is 0.5, but you can change it using the parameter --minbit.
	http://merenlab.org/2016/11/08/pangenomics-v2/#running-a-pangenome-analysis

	A NOTE ON INFATION
	The main handle for changing inflation. This is also the principal handle for regulating cluster granularity.
	Increasing the value inflation will increase cluster granularity. Conceivable values are from 1.1 to 10.0 or so, but the range of suitable values will certainly depend on your input graph. 
	For many graphs, 1.1 will be far too low, and for many other graphs, 8.0 will be far too high. You will have to find the right value or range of values by experimenting, 
	using your judgment, and using measurement tools such as clm dist and clm info. A good set of values to start with is 1.4, 2 and 6.
	https://micans.org/mcl/man/mclfaq.html#faq7.2

	Use the MCL algorithm (van Dongen and Abreu-Goodger, 2012) to identify clusters in amino acid sequence similarity search results. 
	We use 2 as the MCL inflation parameter by default. This parameter defines the sensitivity of the algorithm during the identification of the gene clusters. 
	More sensitivity means more clusters, but of course more clusters does not mean better inference of evolutionary relationships. 
	More information on this parameter and it’s effect on cluster granularity is here http://micans.org/mcl/man/mclfaq.html#faq7.2, but clearly, we, 
	the metagenomics people will need to talk much more about this. So far in the Meren Lab we have been using 2 if we are comparing many distantly related genomes 
	(i.e., genomes classify into different families or farther), and 10 if we are comparing very closely related genomes 
	(i.e., ‘strains’ of the same ‘species’ (based whatever definition of these terms you fancy)). 
	You can change it using the parameter --mcl-inflation. Please experiment yourself, and consider reporting!
	http://merenlab.org/2016/11/08/pangenomics-v2/#running-a-pangenome-analysis

////



#------------------------------------------------------------------------------------------------------------


echo
echo " 		Pangenomics Analysis Begins!!"
echo

# First, create the Internal genomes file with the correct headers that is needed by Anvio
printf '%s\t%s\t%s\t%s\t%s\n' 'name' 'bin_id' 'collection_id' 'profile_db_path' 'contigs_db_path' > Internal_genomes.tsv


# Second, create the file input for the upcoming while loop, one strain per line.
cat > ${StrainsInternal} <<EOF
Strain01
Strain02
Strain03
Strain04
Strain05
Strain06
Strain07
Strain08
Strain09
Strain10
Strain11
Strain12
Strain13
Strain14
Strain16
Strain18
Strain19
Strain20
Strain21
Strain22
Strain23
Strain24
Strain25
EOF


# Populate the Internal genomes file with the correct info. This file is used as input in the upcoming anvi-gen-genomes-storage command.
while read p; do
	StrainX=$p
	StrX=${StrainX/ain/}
	AnvioStrainsDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}
	CntgDB=${StrX}_Pilon_contigs.db
	MergedDB=${StrX}_Merged_Profiles/PROFILE.db
	printf "%s\t%s\t%s\t%s\t%s\n" "${StrainX}" "EVERYTHING" "DEFAULT" "${AnvioStrainsDir}/${MergedDB}" "${AnvioStrainsDir}/${CntgDB}" >> Internal_genomes.tsv
done <${StrainsInternal}
echo 
echo
echo " 		Create a Genomes Storage Database"
echo
# Third, generate a genomes storage database
# As gene caller name we have to set the exact same name we have set when we created the Contigs Databases.
anvi-gen-genomes-storage -i Internal_genomes.tsv -o ${GenomesDatabase} --gene-caller "prodigal"
echo
echo "----------------------GENOMES STORAGE DATABASE WAS CREATED--------------------------------------------"
echo
echo
echo " 			Create the Pangenome!"
echo
echo "This is the Pangenome created using all the Pseudomonas Strains from Columbo, Santorini. Contigs Databases were annotated by Prokka while the gene calls were searched with Prodigal. " > $DescriptionTxt

anvi-pan-genome --genomes-storage ${GenomesDatabase} \
                --project-name "Pseudomonas_Prokka_Pangenome" \
		--description $DescriptionTxt \
                --output-dir ${OutputGenomesDir} \
                --num-threads $SLURM_NTASKS \
                --minbit ${MINBIT} \
                --use-ncbi-blast \
		--mcl-inflation ${MCLInflation} 
echo
echo "--------------------------PANGENOME WAS CREATED!!!-----------------------------------------------"
echo
echo
echo 
echo "			Create an html summary of the Pangenome"
echo
# First we have to create a default collection
anvi-script-add-default-collection -p ${OutputGenomesDir}/Pseudomonas_Prokka_Pangenome-PAN.db

# Then create the actual html file in a directory called SUMMARY.
anvi-summarize -g PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db -p ${OutputGenomesDir}/Pseudomonas_Prokka_Pangenome-PAN.db \
							--init-gene-coverages \
							--collection-name DEFAULT \
							--report-aa-seqs-for-gene-calls \
							--output-dir SUMMARY \
							--just-do-it 

echo
echo "-------------------------------------------HTML SUMMARY WAS CREATED!-------------------------------------"
echo

 



<<////

	Next Steps:
	Download the anvio pangenome results to a GUI enabled environment, and then run the following command:
	
	anvi-display-pan -g PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db

	
	Alternatively, the above command can be invoked from the server but still displayed on the personal computer with the following commands:

	1) First, connect to the server, according to the tutorial "Running remote anvi'o interactive interfaces on your local computer"  here http://merenlab.org/2018/03/07/working-with-remote-interative/
	ssh -L 27000:localhost:27000 zorba | tee /dev/tty | python3 ~/.ssh/run_webbrowser.py

	2) After we are connected to the server, we start an interactive job on the server
	srun -N 1 --ntasks 2 -J Anvio_Interactive --mem 12800 --partition batch --time 08:00:00 --pty bash

	3) After we have started the interactive job, we start the python virtual environment. This is a separate anvio installation we have done, to enable us to run the display-pangenome command.
	# IMPORTANT: It only works through the python virtual environment, even if we have installed the anvio software on our path. 
	source ${HOME}/virtual-envs/anvio-5.3/bin/activate

	4) Go to the Anvio pangenome directory, where the analysis was completed.
	cd ${HOME}/Pseudomonas/Anvio/Pangenome/Internal_genomes/Prokka

	5) Finally, start the anvio pangenome interactive mode.
	anvi-display-pan -g PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db --port-number 27000


////





echo
## Delete files not needed anymore!
#rm $DescriptionTxt Internal_genomes.tsv ${StrainsInternal}
echo
echo "==============================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
