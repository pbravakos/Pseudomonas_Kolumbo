#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="3rd_Anvio_Prokka"
#SBATCH --output=3rd_Anvio_Prokka_job_%j.out

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

	This is the 3rd part of the Anvio Prokka pipeline. 
	Check the tutorial:
	http://merenlab.org/2016/06/22/anvio-tutorial-v2/ "Anvi'o User Tutorial for Metagenomic Workflow"

	This script takes as input one argument. 
	For Strain 01 that would be: 
	script.sh Strain01
	This script runs from the master folder of Anvio and then we change directory to the folder of each specific Strain folder.

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
    echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
    generalInfo >&2
    exit 1
fi

# Check that the working directory exists.
[[ ! -d $1 ]] && { echo "Directory $1 cannot be not found!" >&2; exit 1; }

# Start a job dependency to move the Sdtout file to the correct folders
${HOME}/sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrNum=${StrainX/Strain/}

ReadDir=${HOME}/Pseudomonas/Reads/${StrainX}

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16"|"05"|"08" ]] ; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "01"|"02"|"03"|"04"|"06"|"07" ]]; then
    DateRun=Jun18
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

Paired1=${StrX}_R1_${DateRun}.fastq
Paired2=${StrX}_R2_${DateRun}.fastq
Single1=${StrX}_SE1_${DateRun}.fastq
Single2=${StrX}_SE2_${DateRun}.fastq
MeReads=${StrX}_MeReads_${DateRun}.fastq

Bam_PE_File=${StrX}_paired_end_Pilon_sorted.bam
BAM_SE1_File=${StrX}_SE1_Pilon_sorted.bam
BAM_SE2_File=${StrX}_SE2_Pilon_sorted.bam
BAM_Merged_File=${StrX}_merged_Pilon_sorted.bam

ProfileDirPE=${StrX}_PE_BamProfile
ProfileDirSE1=${StrX}_SE1_BamProfile
ProfileDirSE2=${StrX}_SE2_BamProfile
ProfileDirMeReads=${StrX}_MeReads_BamProfile

CntgDB=${StrX}_Pilon_contigs.db
ContigPilon=${StrX}_Pilon_CLA_Blast.fasta

# Anvio command options that can be altered
ProfileDBMinCngLength=1000 # Default 2500, do not go below 1000


cd ${StrainX}

#--------------------------SANITY CHECKS -------------------------------------------------------------------#
# Check that all the reads can be found
if [[ ! -e ${ReadDir}/$Paired1 ]] && [[ ! -e ${ReadDir}/$Paired1 ]] && [[ ! -e ${ReadDir}/$Single1 ]] && [[ ! -e ${ReadDir}/$Single2 ]] && [[ ! -e ${ReadDir}/$MeReads ]]; then
    echo "Some or all of the fastq files are missing. Please check the folder ${ReadDir}" >&2
    exit 1
fi

# Check that the contigs fasta file exists in the working folder
if [[ ! -e $ContigPilon ]]; then
    echo "The contigs fasta file is not present in the working directory." >&2
    exit 1
fi


## Start a job dependency for the next part of the Anvio pipeline to start as soon as this one has finished
#${HOME}/sbatch --dependency=afterok:$SLURM_JOB_ID 4.Anvio_Prokka_Jun18_PROFILING_AND_HMMs.sh ${StrainX}

echo
echo "				The analysis begins!!" 
echo

# ATTENTION!
# BWA needs the contig file in the same directory where it is running!!

echo '		BWA MAPPING BEGINS'
echo 

# Create bwa index for Pilon
bwa index ${ContigPilon}
# Align reads to contigs with bwa mem.
# Turn sam to bam and finally sort by coordinates
bwa mem -t $SLURM_NTASKS ${ContigPilon} ${ReadDir}/${Paired1} ${ReadDir}/${Paired2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${Bam_PE_File}
# Index the sorted bam file
samtools index ${Bam_PE_File}

bwa mem -t $SLURM_NTASKS ${ContigPilon} ${ReadDir}/${MeReads} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${BAM_Merged_File}
samtools index ${BAM_Merged_File}

bwa mem -t $SLURM_NTASKS ${ContigPilon} ${ReadDir}/${Single1} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${BAM_SE1_File}
samtools index ${BAM_SE1_File}

bwa mem -t $SLURM_NTASKS ${ContigPilon} ${ReadDir}/${Single2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${BAM_SE2_File}
samtools index ${BAM_SE2_File}


echo 
echo '-------------------------BWA MAPPING COMPLETE----------------------------------------'
echo

echo
echo '		Profile BAM files'
anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${Bam_PE_File} \
             --output-dir ${ProfileDirPE} --overwrite-output-destinations \
             --sample-name ${StrX}_PE_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs

anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${BAM_SE1_File} \
             --output-dir ${ProfileDirSE1} --overwrite-output-destinations \
             --sample-name ${StrX}_SE1_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs

anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${BAM_SE2_File} \
             --output-dir ${ProfileDirSE2} --overwrite-output-destinations \
             --sample-name ${StrX}_SE2_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs

anvi-profile --contigs-db ${CntgDB} --num-threads $SLURM_NTASKS --input-file ${BAM_Merged_File} \
             --output-dir ${ProfileDirMeReads} --overwrite-output-destinations \
             --sample-name ${StrX}_MeReads_Profile --min-contig-length ${ProfileDBMinCngLength} \
             --cluster-contigs --profile-SCVs

echo
echo '--------------------------------BAM FILE PROFILING ANALYSIS COMPLETED-----------------------'
echo
echo
# touch ${StrX}_ProfileDB_Anvio_decription.txt
echo "${StrX} Profile database created from the Pilon output scaffold fasta files. Gene finding was done with Prodigal (Prodigal V2.6.3: February, 2016) and annotation in Prokka v1.13.3, using a genus specific Database, after downloading all Pseudomonas complete genomes in genbank format from Refseq. Also genes from GenemarkS2 (http://exon.gatech.edu/GeneMark/genemarks2.cgi) gene caller have been added, only when these gene calls were not found by Prodigal. There are four different profiles. One corresponds to the paired end reads, one to the merged reads (created by BBMerge) and the last two correspond to the singleton reads as were filtered out from the original paired end reads by the program Prinseq. Scaffolds were created only from one Illumina Miseq run!" >  ${StrX}_ProfileDB_Anvio_decription.txt
echo
echo '		Merge BAM PROFILES'
echo

anvi-merge --contigs-db ${CntgDB} --output-dir ${StrX}_Merged_Profiles \
           --sample-name ${StrX}_Merged_Profiles --overwrite-output-destinations \
           --description ${StrX}_ProfileDB_Anvio_decription.txt \
           ${ProfileDirPE}/PROFILE.db ${ProfileDirSE1}/PROFILE.db \
           ${ProfileDirSE2}/PROFILE.db ${ProfileDirMeReads}/PROFILE.db

echo
echo '------------------------MERGING PROFILES COMPLETED----------------------------------------'
echo
echo "			Create a new Collection uniting all the bins into one!"

# In our case probably it is meaningful to unite all bins in one new collection and examine them all together.
anvi-script-add-default-collection -p ${StrX}_Merged_Profiles/PROFILE.db -c ${StrX}_Pilon_contigs.db

# Subsequently, in order to view the new collection use the following command:
# ONLY ON GUI ENVIRONMENTS! anvi-interactive -c {CntgDB} -p ${StrX}_Merged_Profiles/PROFILE.db --gene-mode --collection-name "DEFAULT" --bin-id EVERYTHING

echo 
echo "-------------------------------DEFAULT COLLECTION CREATED!---------------------------------"
echo

rm ${StrX}_ProfileDB_Anvio_decription.txt *.bam* *.fasta.*
echo

<< ////
	NEXT STEPS:
	Download to a local pc the following files and folders (with all their contents): $CntgDB $ProfileDirMeReads
	Start the interactive interface on a gui enabled pc by typing:
	anvi-interactive -p ${StrX}_Merged_Profiles/PROFILE.db -c ${CntgDB} --title ${StrX}_Anvio_results --taxonomic-level t_species

	The next commands are not actually needed! Everything can be done from the interactive interface of the above command!
	Instead of contigs one can also explore the genes AND their annotations:
	anvi-interactive -c{CntgDB} -p ${StrX}_Merged_Profiles/PROFILE.db --gene-mode \
                         --collection-name "CONCOCT" --bin-id Bin_1 #Bin_1 can be substituted with Bin_2 and so on
	To check the available collections and the number of bins:
	anvi-show-collections-and-bins -p {StrX}_Merged_Profiles/PROFILE.db
	anvi-interactive -p {StrX}_Merged_Profiles/PROFILE.db -c {CntgDB} --list-collections
////

echo "======================================================================="
echo "SLURM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
