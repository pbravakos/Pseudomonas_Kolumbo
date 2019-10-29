#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=2
# #SBATCH --mem=128000
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="RGI"
#SBATCH --output=RGI_job_%j.out


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
echo "=================================================="
echo

generalInfo () {
    cat <<END
	
	This is the RGI (Resistance Gene Identifier) pipeline

	This script takes as input one argument. 
	For Strain 01 that would be: 
	sbatch $0 Strain01

	The number of tasks for the Slurm should be equal to the number of srun commands we run as background processes.
	
	We have already downloaded the Card database. Then we loaded the Database to the rgi program.
	cd ${HOME}/Software/rgi-5.1.0	
	wget https://card.mcmaster.ca/download/0/broadstreet-v3.0.4.tar.gz
	tar xvf broadstreet-v3.0.4.tar.gz
	./rgi load -i card.json

	We run the program twice:
	One time on the protein result from Prokka
	A second time on the nucleotide fasta contigs

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] && [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
   generalInfo >&2
   exit 1
fi

# Start a job dependency to move the Sdtout file to the correct folders
sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1

# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}

RGIDir=${HOME}/Software/rgi-5.1.0
ProkkaDir=${HOME}/Pseudomonas/Prokka/${StrainX}/${StrX}_Prokka
AnvioProkkaDir=${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}
ContigDir=${HOME}/Pseudomonas/Seqtk/${StrainX}

# Parameters that can be changed
InputProtein=protein
InputContig=contig
AlignTool=BLAST


echo
echo
echo "			Resistance Gene Identifier Analysis Begins!"
echo

# Create the directory we are going to work on
mkdir -p ${StrainX}

# Change directory to working directory
cd ${StrainX}


srun --ntasks=1 ${RGIDir}/rgi main \
				--input_sequence ${AnvioProkkaDir}/${StrX}_InterProScan_amino-acid-sequences.fasta \
				--output_file ${StrX}_RGI_proteins \
				--input_type ${InputProtein} \
				--clean \
				--alignment_tool ${AlignTool} &

# Run again but change InputSequence and InputType
srun --ntasks=1 ${RGIDir}/rgi main \
				--input_sequence ${ContigDir}/${StrX}_Pilon_CLA_Blast.fasta \
				--output_file ${StrX}_RGI_contigs \
				--input_type ${InputContig} \
				--clean \
				--alignment_tool ${AlignTool} &

# Wait for background jobs to finish
wait

mv ${StrX}_RGI_proteins.txt ${StrX}_RGI_proteins.tsv
mv ${StrX}_RGI_contigs.txt ${StrX}_RGI_contigs.tsv


<<"HERE"

	NEXT STEPS:
	When we have run this script for all strains (and have created the respective files), we can run the following commands in the working directory:
	mkdir Results && cp Str*/St*_RGI_proteins.tsv Results/ && cd Results/
	for i in Str*; do
            sed -i -e "1d" "$i" 
            StrainX=${i/_RGI_proteins.tsv/}
            StrainX=${StrainX/Str/Strain}
            awk -v var=$StrainX '{OFS=FS="\t"; print var, $0}' $i > $i.new
            mv $i.new $i
        done
        cat *.tsv | awk 'BEGIN{OFS=FS="\t"; print "Strain_Code","Gene_Caller_ID","Best_Hit_ARO","ARO","Drug Class","Resistance Mechanism","AMR Gene Family","SNPs_in_Best_Hit_ARO","Other_SNPs"} {print $1,$2,$10,$12,$16,$17,$18,$14,$15}' > RGI_results.tsv
        rm Str*
	
HERE



echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
