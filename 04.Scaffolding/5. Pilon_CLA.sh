#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="Pilon"
#SBATCH --output=Pilon_job_%j.out


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

generalInfo () {
    cat <<END
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	./$0 Strain01
	This script runs from the master folder of Pilon and then we change directory to the folder of each specific Strain folder.
	
	NOTE:
	Here we use the fastq files after Prinseq even if further downstream filtering has been done (e.g. contamination removal) because we want as much information as possible to correct our contigs!!

END
}


# Check that an argument has been given in the correct form.
if [[ $# -ne 1 ]] || [[ ! $1 =~ ^Strain[0-9]{2}$ ]]; then
   echo "Exactly one argument should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9!" >&2
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

PilonDir="${HOME}/Software/pilon"
SealerDir="${HOME}/Pseudomonas/Sealer/${StrainX}"
PrinseqDir="${HOME}/Pseudomonas/Prinseq/${StrainX}/Filtering"
ContigSealer=${StrX}_CLA_Sealer_scaffold.fa

PESorted=${StrX}_PE_Pilon_sorted.bam
SE1Sorted=${StrX}_SE1_Pilon_sorted.bam
SE2Sorted=${StrX}_SE1_Pilon_sorted.bam

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16"|"05"|"08" ]] ; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
else
    echo "$StrainX has not been defined in the if else bash expression. Please check it to resolve the issue."  >&2
    exit 1
fi

PE1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
PE2=${StrX}_${DateRun}_prinseq_good_R_2.fastq
SE1=${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq
SE2=${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq


# Parameters to Pilon which can be changed
FlankBases=1 # Controls how much of the well-aligned reads will be used; this many bases at each end of the good reads will be ignored (default 10).
GapMarginBases=2000   # Closed gaps must be within this number of bases of true size to be closed (100000)
Kmer=77    # Kmer size used by internal assembler (default 47)
#--------------------------------------------------------------------------------------------------------------------------

cd ${StrainX}

# Create soft link to contig fasta files. This is done because contig fasta files need to be in the same folder as the bwa.sh file for the pipeline to work!
ln -s ${SealerDir}/${ContigSealer}

echo
# Print Bwa version
( bwa 3>&1 1>&2- 2>&3- ) | head -n 3
echo

echo
# Print samtools version
( samtools 3>&1 1>&2- 2>&3- ) | head -n 3
echo

# Create bwa index for Spades
bwa index ${ContigSealer}
# Align reads to contigs with bwa mem and bwa bwasw. bwasw seems not to work well with paired end reads. For this reason we use bwa mem.
#Turn sam to bam and finally sort by coordinates
bwa mem -t $SLURM_NTASKS ${ContigSealer} ${PrinseqDir}/${PE1} ${PrinseqDir}/${PE2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${PESorted}
#Index the sorted bam file
samtools index ${PESorted}

#Repeat the same pipeline for the single end reads.
#1st
bwa mem -t $SLURM_NTASKS ${ContigSealer} ${PrinseqDir}/${SE1} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${SE1Sorted}
#Index the sorted bam file
samtools index ${SE1Sorted}

#2nd
bwa mem -t $SLURM_NTASKS ${ContigSealer} ${PrinseqDir}/${SE2} |\
samtools view -hu -@ $SLURM_NTASKS - |\
samtools sort -@ $SLURM_NTASKS - -o ${SE2Sorted}
#Index the sorted bam file
samtools index ${SE2Sorted}

# Merge the single end bam files
#samtools merge {PESorted} ${SE1Sorted} ${SE2Sorted}

#Remove files that are not needed for downstream analysis 
rm *.fasta.*


java -Xmx${SLURM_MEM_PER_NODE}M -jar ${PilonDir}/pilon-1.23.jar --genome ${ContigSealer} \
								--frags ${PESorted} \
								--unpaired ${SE1Sorted} \
								--unpaired ${SE2Sorted} \
								--output ${StrX}_Pilon_ \
								--K $Kmer \
								--flank $FlankBases \
								--gapmargin $GapMarginBases \
								--fix "gaps","local","amb","breaks" \
								--iupac  \
								--threads $SLURM_NTASKS \
								--verbose \
								--debug \
								--tracks
							
# Remove files not needed anymore
rm *.bam* *.fa.*

# Rename the Pilon fasta scaffolds output to something better!
mv ${StrX}_Pilon_.fasta ${StrX}_Pilon_CLA.fasta


# Repeat the process for those strains that have been sequenced twice, using this time the reads from another library, and use the output from th
if [[ $StrNum =~ "05"|"08" ]]; then 
    DateRun=Jun18
    PE1=${StrX}_${DateRun}_prinseq_good_R_1.fastq
    PE2=${StrX}_${DateRun}_prinseq_good_R_2.fastq
    SE1=${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq
    SE2=${StrX}_${DateRun}_prinseq_good_R_1_singletons.fastq
    # Rename the output of the previous library
    mv ${StrX}_Pilon_CLA.fasta ${StrX}_Pilon_CLA_old.fasta
    # Create bwa index for Spades
    bwa index ${StrX}_Pilon_CLA_old.fasta
    # Align reads to contigs with bwa mem and bwa bwasw. bwasw seems not to work well with paired end reads. For this reason we use bwa mem.
    # Turn sam to bam and finally sort by coordinates
    bwa mem -t $SLURM_NTASKS ${StrX}_Pilon_CLA_old.fasta ${PrinseqDir}/${PE1} ${PrinseqDir}/${PE2} |\
    samtools view -hu -@ $SLURM_NTASKS - |\
    samtools sort -@ $SLURM_NTASKS - -o ${PESorted}
    # Index the sorted bam file
    samtools index ${PESorted}

    # Repeat the same pipeline for the single end reads.
    # 1st
    bwa mem -t $SLURM_NTASKS ${StrX}_Pilon_CLA_old.fasta ${PrinseqDir}/${SE1} |\
    samtools view -hu -@ $SLURM_NTASKS - |\
    samtools sort -@ $SLURM_NTASKS - -o ${SE1Sorted}
    # Index the sorted bam file
    samtools index ${SE1Sorted}
    # 2nd
    bwa mem -t $SLURM_NTASKS ${StrX}_Pilon_CLA_old.fasta ${PrinseqDir}/${SE2} |\
    samtools view -hu -@ $SLURM_NTASKS - |\
    samtools sort -@ $SLURM_NTASKS - -o ${SE2Sorted}
    # Index the sorted bam file
    samtools index ${SE2Sorted}

    #Remove files that are not needed for downstream analysis 
    rm *.fasta.*


    java -Xmx${SLURM_MEM_PER_NODE}M -jar ${PilonDir}/pilon-1.23.jar --genome ${StrX}_Pilon_CLA_old.fasta \
								--frags ${PESorted} \
								--unpaired ${SE1Sorted} \
								--unpaired ${SE2Sorted} \
								--output ${StrX}_Pilon_ \
								--K $Kmer \
								--flank $FlankBases \
								--gapmargin $GapMarginBases \
								--fix "gaps","local","amb","breaks" \
								--iupac  \
								--threads $SLURM_NTASKS \
								--verbose \
								--debug \
  								--tracks
    
    
    # Remove files not needed anymore
    rm *.bam* *.fa.*
    # Rename the Pilon fasta scaffolds output to something better! 
    mv ${StrX}_Pilon_.fasta ${StrX}_Pilon_CLA.fasta 

    # Remove excess information from fasta headers that is not needed. This will be the final output of our Pilon analysis!
    sed -i 's/pilon_pilon/pilon/g' ${StrX}_Pilon_CLA.fasta
fi

echo "==============================="
echo "SLRUM job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
