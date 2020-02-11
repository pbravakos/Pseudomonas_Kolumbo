#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Pangenome_binning"
#SBATCH --output=Pangenome_binning_job_%j.out

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

<< ////

        This is the Anvio pangenomics pipeline. Check the tutorial:
        http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions "Creating a quick pangenome with functions"
	
	This script runs in the Anvio Pangenomics Prokka directory. 
	It creates a new collection for the pangenome profile database, with the bins that are relevant to this work. 
	The bins of the new collection are used in the next part of this pipeline to do the functional enrichment analysis.
	To run the script:
	sbatch $0

	IMPORTANT
        To run this script we need to have created the "SUMMARY" directory in the previous step of this pipeline. 

	IMPORTANT
	We have to modify the current script before running it, in order to get the bins that are of interest to us!

////


# INITIAL PARAMETERS
SummaryDir=SUMMARY # This name has to be the same, as the output of the previous part of this pipeline.
SummaryFile=Pseudomonas_Prokka_Pangenome_gene_clusters_summary.txt
ProfileDB=PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db

# Next, we have to set the name of the collection and the respective genes. These names will be used in the next part of this pipeline.
CollectionName="Pseudomonas"
BinName1="Core"
BinName2="Stutzeri_unique"
BinName3="Aeruginosa_unique"
BinName4="Aeruginosa_unique_without_Strain05"



#------------------------------------------------------------------------------------------------------------------------------

# Check for the existence of the output of the previous part of this pipeline.
if [[ ! -d ${SummaryDir} ]]; then
    echo "${SummaryDir} cannot be found. Please run the first part of this pipeline first to create this directory!"
    exit 1
fi


# Extract the file.
if [[ ! -e ${SummaryDir}/${SummaryFile} ]]; then
    gunzip ${SummaryDir}/${SummaryFile}.gz
fi


# First, time find the core pangenome GCs i.e. the GCs that are common to all strains.
awk -v bin="${BinName1}" 'BEGIN{FS="\t"; OFS="\t"} $6==21 {print $2, bin}' ${SummaryDir}/${SummaryFile} | sort -u  > ${BinName1}_GC_binning.tsv

# First we want to get only the GCs that are related to the 13 stutzeri strains and can be found in exactly 13 genomes. The vast majority of these GCs will belong to stutzeri only but it is possible that some will be found in other strains also.
awk -v bin="${BinName2}" 'BEGIN{FS="\t"; OFS="\t"} $6==13 {print $2, bin}' ${SummaryDir}/${SummaryFile} | sort -u  > ${BinName2}_initial_binning.tsv

# Next, we want to create a file with all the GCs that do not belong to stutzeri strains.
awk 'BEGIN{FS="\t"; OFS="\t"} $4=="Strain05" || $4=="Strain07" || $4=="Strain23" || $4=="Strain19" || $4=="Strain20" || $4=="Strain06" || $4=="Strain16" || $4=="Strain21" {print $2}' ${SummaryDir}/${SummaryFile} | sort -u  > All_Aeruginosa_GCs

# Finally, we want to combine the two files created above, to get only the GCs belonging exclusively to stutzeri strains.
# For this reason we use grep with all Aeruginosa GCs as a pattern to search in the stutzeri GCs, and get the non matching GCs as output.
grep -vFf All_Aeruginosa_GCs ${BinName2}_initial_binning.tsv > ${BinName2}_GC_binning.tsv


# Next, repeat for the aeruginosa strains!
# First we want to get only the GCs that are related to the 8 aeruginosa strains and can be found in exactly 8 genomes. The vast majority of these GCs will belong to aeruginosa only but it is possible that some will be found in other strains also.
awk -v bin="${BinName3}" 'BEGIN{FS="\t"; OFS="\t"} ($4=="Strain05" || $4=="Strain07" || $4=="Strain23" || $4=="Strain19" || $4=="Strain20" || $4=="Strain06" || $4=="Strain16" || $4=="Strain21") && $6==8 {print $2, bin}' ${SummaryDir}/${SummaryFile} | sort -u  > ${BinName3}_initial_binning.tsv

# Next, we want to create a file with all the GCs that do not belong to aeruginosa strains.
awk 'BEGIN{FS="\t"; OFS="\t"} $4!="Strain05" && $4!="Strain07" && $4!="Strain23" && $4!="Strain19" && $4!="Strain20" && $4!="Strain06" && $4!="Strain16" && $4!="Strain21" {print $2}' ${SummaryDir}/${SummaryFile} | sort -u  > Non_Aeruginosa_GCs

# Finally, we want to combine the two files created above, to get only the GCs belonging exclusively to aeruginosa strains.
# For this reason we use grep with the non aeruginosa strains as a pattern to search in the aeruginosa GCs, and get the non matching GCs as output.
grep -vFf Non_Aeruginosa_GCs ${BinName3}_initial_binning.tsv > ${BinName3}_GC_binning.tsv

# Repeat the same procedure, but this time exclude one strain (new total is 7) from the aeruginosa binning.
awk -v bin="${BinName4}" 'BEGIN{FS="\t"; OFS="\t"} ($4=="Strain05" || $4=="Strain07" || $4=="Strain23" || $4=="Strain19" || $4=="Strain20" || $4=="Strain06" || $4=="Strain16" || $4=="Strain21") && $6==7 {print $2, bin}' ${SummaryDir}/${SummaryFile} | sort -u  > ${BinName4}_initial_binning.tsv

grep -vFf Non_Aeruginosa_GCs ${BinName4}_initial_binning.tsv > ${BinName4}_GC_binning.tsv



# Concatenate all the bins into one file.
cat ${BinName1}_GC_binning.tsv ${BinName2}_GC_binning.tsv ${BinName3}_GC_binning.tsv ${BinName4}_GC_binning.tsv > Binning.tsv


# Finally, we import these bins as a new collection named ${CollectionName} in the pangenome.
anvi-import-collection -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db -C ${CollectionName} Binning.tsv


awk 'BEGIN{FS="\t"; OFS="\t"} $4=="Strain14" || $4=="Strain24" || $4=="Strain01" || $4=="Strain12" || $4=="Strain03" || $4=="Strain02" || $4=="Strain18" || $4=="Strain04" || $4=="Strain11" || $4=="Strain22" || $4=="Strain08" || $4=="Strain09" || $4=="Strain10" {print $2}' ${SummaryDir}/${SummaryFile} | sort -u  > All_Stutzeri_GCs


# Print the number of GCs in each bin.
echo "==============================================="
echo
echo "The number of ${BinName1} GCs is:" `wc -l ${BinName1}_GC_binning.tsv | awk '{print $1}'`
echo "The number of ${BinName2} GCs is:" `wc -l ${BinName2}_GC_binning.tsv | awk '{print $1}'`
echo "The number of ${BinName3} GCs is:" `wc -l ${BinName3}_GC_binning.tsv | awk '{print $1}'`
echo "The number of ${BinName4} GCs is:" `wc -l ${BinName4}_GC_binning.tsv | awk '{print $1}'`
echo "The number of all Aeruginosa GCs is:" `wc -l All_Aeruginosa_GCs | awk '{print $1}'`
echo "The number of all Stuzteri GCs is: " `wc -l All_Stutzeri_GCs | awk '{print $1}'`
echo
echo "==============================================="


# Clear files not needed any more.
rm ${BinName1}_initial_binning.tsv ${BinName2}_initial_binning.tsv ${BinName3}_initial_binning.tsv

# Next steps:
# We should go to interactive mode with the command
# anvi-display-pan -g PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db
# and check visually the created bin, to be sure it is what we expected! 

echo
echo "==============================="
echo "job finished " "$(date)"
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
