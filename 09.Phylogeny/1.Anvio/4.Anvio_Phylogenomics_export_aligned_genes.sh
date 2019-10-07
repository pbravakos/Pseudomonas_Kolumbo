#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="4th_Anvio_Phylo"
#SBATCH --output=4th_Anvio_Phylo_job_%j.out

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

	This script is the 4th part of the Anio Phylogenomics pipeline and runs in the master folder of Anvio Phylogenomics directory. The input of this script is the "Pseudomonas_external_genomes.tsv" which is the ouput from a previous step of this pipeline. The genes to search for are also selected in the previous step of this pipeline.
        
        In order to run it standallone the command should be like this:
        sbatch $0
	
	IMPORTANT!
	In this script the number of the tasks reserved in SLURM equals the number of genes we are going to align.
	Each srun command will run in the background, in order to allow them to run concurrently

END
}


# INITIAL PARAMETERS
AnvioGenesDir=${HOME}/Pseudomonas/Anvio/Phylogenomics/Anvio_aligned_genes


# Check that the working directory exists.
[[ ! -d $AnvioGenesDir ]] && { echo "Directory ${AnvioGenesDir} cannot be not found!" >&2; exit 1; }

if [[ ! -e Pseudomonas_external_genomes.tsv ]]; then
    echo "The file 'Pseudomonas_external_genomes.tsv' cannot be found in the working directory ${AnvioGenesDir}. Please run the previous steps of this pipeline first!" >&2
    generalInfo >&2
    exit 1
fi

# First, we find the total number of Genome databases we have prepared
TotalLines=`wc -l Pseudomonas_external_genomes.tsv | cut -f1 -d' '`
# We have to subtract 1 (one) from the total lines because we don't want to take into account the header line.
GenomesTotal=`expr $(($TotalLines - 1))`

# We substitute the command anvi-get-sequences-for-hmm-hits with a parameter, but this is certainly not necessary!
anviGetSeq='anvi-get-sequences-for-hmm-hits'

# Next, we declare an array with all the gene names.
declare -a arr=("POG090903HI" "POG090902MV" 
                "POG09090004" "POG090900O9" 
                "POG09090188" "POG090901ZM" 
                "POG090900YW" "POG090901TG" 
                "POG0909039Q" "POG0909032E" 
                "POG090901TJ" "POG090901ZB"
                "POG09090095" "POG090900UZ"
                "POG0909008M" "POG090900IU"
		"POG090903GW" "POG0909008P"
		"POG09090288" "POG090901JA"
		"POG0909006H" "POG0909001N"
		"POG090900YH" "POG090903CH"
		"POG090901HE" "POG090901H1"
		"POG090900CF" "POG090900MN"
		"POG090901OY" "POG090900GC"
		"POG090903B5" "POG0909024W"
		"POG090900OV" "POG0909001O"
		"POG09090325" "POG090900ZF"
		"POG0909008D" "POG0909006C"
		"POG0909002C" "POG090901XG"
		"POG090900YI" "POG090900CU"
		"POG090900OV" "POG0909001O"
		"POG09090325" "POG090900ZF"
		"POG0909008D" "POG0909006C"
		"POG0909002C" "POG090901XG"
		"POG090900YI" "POG09090384"
		"POG0909028G" "POG090901S2"
		"POG090900RW" "POG0909000G" 
		"POG090900FV" "POG090903FS"
		"POG09090273" "POG090903FS"
		"POG090902FY" "POG0909018I"
		"POG090901ND" "POG090900GE"
		"POG090901RJ" "POG090902JS"
		"POG0909000T" "POG090903DD" 
		"POG090901L6" "POG090903LP"
		"POG090901XJ"
                 )


# Next, we run the srun command inside a for loop. For each gene of the array, we run the "anviGetSeq" command.
for i in "${arr[@]}"
do
   ${HOME}/srun --ntasks 1 --output Stdout_Aligned_${i}.out ${anviGetSeq} \
                              --external-genomes Pseudomonas_external_genomes.tsv \
                              --output-file ${AnvioGenesDir}/Pseudomonas_${i}.fasta \
                              --hmm-sources gammaproteobacteria_anvio \
                              --gene-names ${i} \
                              --min-num-bins-gene-occurs ${GenomesTotal} \
                              --return-best-hit \
                              --align-with muscle \
                              --get-aa-sequences &
done


# Wait till all the background processes are completed before moving to the next part.
wait






echo "==========================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
