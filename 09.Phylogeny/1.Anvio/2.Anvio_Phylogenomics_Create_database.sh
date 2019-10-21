#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="2nd_PhyloGen"
#SBATCH --output=2nd_PhyloGen_job_%j.out

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

#set -euo pipefail # Check for more info at https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a

generalInfo () {
    cat <<END

	This script runs in the master folder of Anvio Phylogenomics directory
        as part of another script which calls this one with a simple while loop 
        and assigns one fasta formatted genome as the first argument of this script.
        
        In order to run it standallone the command we should give as an argument one strain name for which we want to build the database, like this:
        sbatch $0 Pseudomonas_stutzeri_1W1-1A.fna

        WARNING!!
        It is advised fasta files names to have only ONE dot (.) and 
        the extension to be different than ".fasta", for example it could be ".fna" 

        This is the Anvio phylogenomics pipeline. Check the tutorials:
        http://merenlab.org/2017/06/07/phylogenomics/ "An anvi'o workflow for phylogenomics" and
        http://merenlab.org/tutorials/infant-gut/#chapter-iii-phylogenomics "Chapter III: Phylogenomics"
         
        We downloaded on 30-09-2018 the Complete Genomes Genomic DNA (Fasta) from pseudomonas.com
        wget http://www.pseudomonas.com/downloads/pseudomonas/pgd_r_17_2/Pseudomonas/complete/fna-complete.tar.gz
	
	IMPORTANT!
	We selected most of the strains based on the phylogenetic information found at the Genome Taxonomy Database Release 03-RS86 (19th August 2018) http://gtdb.ecogenomic.org/

	IMPORTANT!
	Before running this script, all fasta files should be downloaded from a source e.g. NCBI RefSeq, Pseudomonas Genome Database etc.

END
}

# Check that exactly one argument has been given.
if [[ $# -ne 1 ]] ; then
    echo "Exactly one argument should be given!" >&2
    generalInfo >&2
    exit 1
fi

echo
# Start a job dependency to move the Sdtout file to the correct folders
"${HOME}"/sbatch --dependency=afterany:"$SLURM_JOB_ID" Rename_Slurm_output_files.sh "$SLURM_JOB_NAME" "$SLURM_JOB_ID" "$1"
echo

# INITIAL PARAMETERS

# First we asign to $ContigBadFormat the name of the fasta file as it can be found in the Pseudomonas fna database
ContigBadFormat=$1
# Next we assign to $StrainX the name of the fasta file without the .fna extension!
StrainX=${ContigBadFormat%%.*} 
# Next we asign to $Contig the $StrainX name with the fasta extension 
Contig=${StrainX}.fasta
# Finally we assign to $CntgDB the $StrainX name with the .db extension.
CntgDB=${StrainX}.db

# This is the directory were we have downloaded all the fasta files.
ContigDir="${HOME}"/Pseudomonas/Pseudomonas_database/Strains_Downloaded_from_NCBI


# This is the directory where we are going to run this script!
AnvioPhyloDir="${HOME}"/Pseudomonas/Anvio/Phylogenomics 

# This is the directory where all the analysis will take place!
AnvioDBDir="${AnvioPhyloDir}"/Pseudomonas_Anvio_db

BuscoSingleCopy="${HOME}"/Software/busco_hmms_for_anvio/gammaproteobacteria_anvio


echo ' 				Phylogenomics analysis begins!!!'
echo
# Create the external genomes tab delimited file if it doesn't already exist.
if [[ ! -e Pseudomonas_external_genomes.tsv ]]; then
    printf "%s\t%s\n" "name" "contigs_db_path" > Pseudomonas_external_genomes.tsv
fi

# Here we append to the external_genomes tab delimited file that is needed in the downstream analysis
printf "%s\t%s\n" "${StrainX}" "${AnvioDBDir}"/"${CntgDB}" >> Pseudomonas_external_genomes.tsv

# Next we change directory to where the Contig Database files will be stored!
cd "${AnvioDBDir}"
echo
# Create a soft link of the non-formatted contigs fasta file
ln -s ${ContigDir}/${ContigBadFormat}

echo
echo "				REFORMAT FASTA FILE"
echo
anvi-script-reformat-fasta --simplify-names --output-file ${Contig} ${ContigBadFormat}
echo 
echo
echo '-------------------------------------FASTA FILE WAS REFORMATTED!-----------------------------------------'
echo
echo '				Create Contigs DATABASE'

# Create the Anvio description file
touch ${StrainX}_Anvio_decription.txt
echo "${StrainX} fasta file was downloaded from the online Genome database of Pseudomonas (www.pseudomonas.com) or from NCBI. Gene finding was accomplished by Prodigal." >  ${StrainX}_Anvio_decription.txt


echo
anvi-gen-contigs-database -f ${Contig} --project-name ${StrainX}_Anvio \
                                                --output-db-path ${CntgDB} \
                                                --description ${StrainX}_Anvio_decription.txt

echo
echo
echo '----------------------------------CONTIGS DATABASE COMPLETED----------------------------------------'
echo
echo '			Parse HMM hits to the contigs database'


# We can use the program busco_hmms_for_anvio (https://github.com/guyleonard/busco_hmms_for_anvio) to download 
# (we did it on 29/09/2018) and extract the gammaproteobacteria hmm busco dataset of single copy genes found here: 
# http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz
# This dataset can be used as input to the anvi-run-hmms command
anvi-run-hmms -c "${CntgDB}" --num-threads "$SLURM_NTASKS" -H "$BuscoSingleCopy"
echo
anvi-run-hmms -c "${CntgDB}" --num-threads "$SLURM_NTASKS" -I Ribosomal_RNAs
echo
anvi-run-hmms -c ${CntgDB} --num-threads "$SLURM_NTASKS" -I Campbell_et_al
echo 

echo '-------------------------------HMM HITS ANALYSIS COMPLETE-------------------------------------'

rm ${StrainX}_Anvio_decription.txt ${Contig} ${ContigBadFormat}


<<NOTE1

	NEXT STEPS
	First, check which single copy genes hmms exist for the contig databases we created
	anvi-get-sequences-for-hmm-hits --external-genomes Pseudomonas_external_genomes.tsv --list-hmm-sources
	Run commands to check which single copy genes exist in all genomes
	For example one should run both:
	anvi-get-sequences-for-hmm-hits --external-genomes Pseudomonas_external_genomes.tsv \
	                                --output-file Pseudomonas_concatenated_30S_proteins.fasta \
	                                --hmm-sources gammaproteobacteria_anvio \
	                                --gene-names POG090903HI,POG090902MV,POG090900O9,POG09090188,POG090900YW,POG090901TG,POG090901ZM,POG0909039Q,POG0909032E \
	                                --min-num-bins-gene-occurs 107 \
	                                --return-best-hit \
	                                --concatenate-genes \
	                                --align-with muscle \
	                                --get-aa-sequences
	and
	anvi-get-sequences-for-hmm-hits --external-genomes Pseudomonas_external_genomes.tsv \
	                                --output-file Pseudomonas_concatenated_30S_proteins.fasta \
	                                --hmm-sources gammaproteobacteria_anvio \
	                                --gene-names POG090903HI,POG090902MV,POG090900O9,POG09090188,POG090900YW,POG090901TG,POG090901ZM,POG0909039Q,POG0909032E \
	                                --max-num-genes-missing-from-bin 0 \
	                                --return-best-hit \
	                                --concatenate-genes \
	                                --align-with muscle
	To check (from the warnings!) which genes are not included in which genomes (or bins in Anvio's parlance). 
	Note: The total number of genome databases in the above example is 107! 

NOTE1

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
