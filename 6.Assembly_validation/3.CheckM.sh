#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem-per-cpu=20000
#SBATCH --job-name="CheckM"
#SBATCH --output=CheckM_job_%j.out

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
	
	This script takes as input one argument. 
	For Strain 01 that would be: 
	$0 Strain01
	This script runs from the master folder of CheckM and then we change directory to the folder of each specific Strain folder.
	
	In order to run this bash file for different strains changes have to take place:
	a) In the initial parameters given in the start of this file like the fastq read files (and associated changes of directory structure in each command).
	b) In the taxon_set command parameters for the taxonomy analysis (i.e. change the genus or species). 
	c) Someone could also change (for better or worse!) the individual parameters of each given command.

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
sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID $1


# INITIAL PARAMETERS
StrainX=$1
StrX=${StrainX/ain/}
StrNum=${StrainX/Strain/}


CheckMDir=${HOME}/Software/CheckM-1.0.13/bin

ReadDir=${HOME}/Pseudomonas/Reads/${StrainX}
ContigDir=${HOME}/Pseudomonas/Seqtk/${StrainX}

# Next we wil define the $DateRun parameter based on the Strain number input.
if [[ $StrNum =~ "09"|"10"|"11"|"12"|"13"|"14"|"16" ]] ; then
    DateRun=Nov18a
elif [[ $StrNum =~ "18"|"19"|"20"|"21"|"22"|"23"|"24"|"25" ]]; then
    DateRun=Nov18b
elif [[ $StrNum =~ "01"|"02"|"03"|"04"|"06"|"07"|"05"|"08" ]]; then
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

Bam_PE_File=${StrX}_PE_Pilon_sorted.bam
BAM_SE1_File=${StrX}_SE1_Pilon_sorted.bam
BAM_SE2_File=${StrX}_SE2_Pilon_sorted.bam
BAM_Merged_File=${StrX}_MeReads_Pilon_sorted.bam

ContigPilon=${StrX}_Pilon_CLA_Blast.fasta
StrXPilon=${StrX}_Pilon

genus='Pseudomonas'
# species='Pseudomonas stutzeri'

BinFolder=${StrX}_bin_folder
TempFolder=${StrX}_temp
PlotFolder=${StrX}_plot_folder
TaxFolder=${StrX}_taxonomy_output
LinFolder=${StrX}_lineage_output
TreeFolder=${StrX}_output_tree

CovFile=${StrX}_Pilon.coverage
LinMarkers=${StrX}_Pilon_lineage.markers
TaxMarkers=${StrX}_Pilon_taxonomy.markers
TetraFile=${StrX}_Pilon.tetra


#-------------------------------------------------------------------------------------------------
cd ${StrainX}

mkdir ${TempFolder}
mkdir ${BinFolder}
cd ${BinFolder}
ln -s ${ContigDir}/${ContigPilon} ${StrXPilon}.fna
cd ..
# We also need the same Pilon fasta at the working directory for BWA to work!
ln -s ${ContigDir}/${ContigPilon}


echo '		CHECKM PIPELINE FOR GENOME SCAFFOLDS'
echo

echo 
echo '		BWA MAPPING BEGINS'
echo 
# Create bwa index for Pilon
bwa index ${ContigPilon}
# Align reads to contigs with bwa mem and bwa bwasw.
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
echo '!!!!!!!!!!!!!!!BWA MAPPING COMPLETE!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		GC PLOT ANALYSIS BEGINS'
echo 

# ATTENTION!! IMPORTANT!!!
# matplotilib version HAS to be less than 2.2 for this to run. Except for some plots (like for tetranuclotides) for which the version HAS to be 1.3.1!!!!!!


python2 ${CheckMDir}/checkm gc_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 8 --gc_window_size 2000 --gc_bin_width 0.01 ${BinFolder} ${PlotFolder} 50

echo 
echo '!!!!!!!!!!!!!!!ANALYSIS OF GC PLOT COMPLETE!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		Nx ANALYSIS BEGINS'
echo 

python2 ${CheckMDir}/checkm nx_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --step_size 0.01 ${BinFolder} ${PlotFolder}

echo ''
echo '!!!!!!!!!!!!Nx ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!'
echo ''
echo ''
echo '		Sequence Length ANALYSIS BEGINS'
echo ''

python2 ${CheckMDir}/checkm len_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder}

echo 
echo '!!!!!!!!!!!!Length ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		LENGTH HISTOGRAM ANALYSIS BEGINS'
echo 

LC_ALL=C python2 ${CheckMDir}/checkm len_hist --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder}

echo ''
echo '!!!!!!!!!!!!LENGTH HISTOGRAM ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo ''
echo ''
echo '		COVERAGE ANALYSIS BEGINS'
echo ''

LC_ALL=C python2 ${CheckMDir}/checkm coverage --all_reads --threads $SLURM_NTASKS ${BinFolder} ${CovFile} ${Bam_PE_File} ${BAM_SE1_File} ${BAM_SE2_File} ${BAM_Merged_File}

echo 
echo '!!!!!!!!!!!!COVERAGE ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		TETRANUCLEOTIDE ANALYSIS BEGINS'
echo 

LC_ALL=C python2 ${CheckMDir}/checkm tetra -t 1 ${BinFolder}/${StrXPilon}.fna ${TetraFile}

echo ''
echo '!!!!!!!!!!!!!!!!TETRA ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!'
echo ''
echo '		TREE ANALYSIS BEGINS'
echo ''

LC_ALL=C python2 ${CheckMDir}/checkm tree --ali --nt --threads $SLURM_NTASKS --pplacer_threads $SLURM_NTASKS --tmpdir ${TempFolder} ${BinFolder} ${TreeFolder}
echo ''
LC_ALL=C python2 ${CheckMDir}/checkm tree_qa --out_format 1 --file ${StrXPilon}_tree_format1.tab --tab_table --tmpdir ${TempFolder} ${TreeFolder}
echo ''
LC_ALL=C python2 ${CheckMDir}/checkm tree_qa --out_format 2 --file ${StrXPilon}_tree_format2.tab --tab_table --tmpdir ${TempFolder} ${TreeFolder}
echo ''
# Format 3 is nwk tree but tips have only img genome ids.
LC_ALL=C python2 ${CheckMDir}/checkm tree_qa --out_format 3 --file ${StrXPilon}_tree_format3.nwk --tmpdir ${TempFolder} ${TreeFolder}
echo ''
# Format 3 is nwk tree but tips of the tree have img genome ids and genome names
LC_ALL=C python2 ${CheckMDir}/checkm tree_qa --out_format 4 --file ${StrXPilon}_tree_format4.nwk --tmpdir ${TempFolder} ${TreeFolder}
echo ''
LC_ALL=C python2 ${CheckMDir}/checkm tree_qa --out_format 5 --file ${StrXPilon}_tree_format5.msa --tmpdir ${TempFolder} ${TreeFolder}

echo 
echo '!!!!!!!!!!!!!TREE ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		LINEAGE ANALYSIS BEGINS'

LC_ALL=C python2 ${CheckMDir}/checkm lineage_set --tmpdir ${TempFolder}  --unique 10 --multi 10 ${TreeFolder} ${LinMarkers} 

LC_ALL=C python2 ${CheckMDir}/checkm analyze --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${BinFolder} ${LinFolder}

# AAI stands for Amino Acid Identity in the checkm qa manual. Look here https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands
LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 1 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format1.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 2 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format2.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 3 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format3.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 4 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format4.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 5 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format5.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 6 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format6.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 7 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format7.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 8 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format8.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 9 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_lineage_qa_format9.fasta --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${LinMarkers} ${LinFolder}

echo 
echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!LINEAGE ANALYSIS COMPLETE!!!!!!!!!!!!!!!'
echo 
echo 
echo '		TAXONOMY ANALYSIS BEGINS'
echo 

# The list of available taxonomic-specific marker sets which can be inserted as input in the taxon_set command can be viewed with the taxon_list command.

LC_ALL=C python2 ${CheckMDir}/checkm taxon_set --tmpdir ${TempFolder} genus ${genus} ${TaxMarkers}
# A more specialized command depending on the specific species under analysis would be to specify a species name instead of a genus e.g.:
# NOT NEDED! LC_ALL=C ${CheckMDir}/checkm taxon_set --tmpdir ${TempFolder} species ${species} ${TaxMarkers}

LC_ALL=C python2 ${CheckMDir}/checkm analyze --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${BinFolder} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 1 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format1.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 2 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format2.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 3 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format3.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 4 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format4.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 5 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format5.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 6 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format6.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 7 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format7.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 8 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format8.tab --tab_table --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

LC_ALL=C python2 ${CheckMDir}/checkm qa --out_format 9 --e_value 1e-10 --length 0.7 --coverage_file ${CovFile} --file ${StrXPilon}_taxonomy_qa_format9.fasta --threads $SLURM_NTASKS --tmpdir ${TempFolder} ${TaxMarkers} ${TaxFolder}

echo 
echo '!!!!!!!!!!!!!!!!!!!TAXONOMY ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		PLOTS FOR completeness AND contamination ANALYSIS BEGINS'
echo 

LC_ALL=C python2 ${CheckMDir}/checkm bin_qa_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --row_height 0.3 ${TaxFolder} ${BinFolder} ${PlotFolder}
cd ${PlotFolder}
mv bin_qa_plot.pdf ${StrXPilon}_taxonomy_Contamination.pdf
cd ..
echo 
LC_ALL=C python2 ${CheckMDir}/checkm bin_qa_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --row_height 0.3 ${LinFolder} ${BinFolder} ${PlotFolder}
cd ${PlotFolder}
mv bin_qa_plot.pdf ${StrXPilon}_lineage_Contamination.pdf
cd ..
echo ''
echo '!!!!!!!!!!!!!!!PLOTS FOR completeness AND contamination ANALYSIS COMPLETE!!!!!!!!!!!!!'
echo 
echo 
echo '		Coding plot analysis BEGINS'
echo 

LC_ALL=C python2 ${CheckMDir}/checkm coding_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --cd_window_size 5000 --cd_bin_width 0.01 ${TaxFolder} ${BinFolder} ${PlotFolder} 50

echo 
echo '!!!!!!!!!!!!!CODING PLOT ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		PLOT OF GC and COVERAGE BEGINS'

LC_ALL=C python2 ${CheckMDir}/checkm par_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${TaxFolder} ${BinFolder} ${PlotFolder} ${CovFile}

echo 
echo '!!!!!!!!!!!!!!!!!!!PLOT OF GC and COVERAGE COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		GC, CD, and TD distribution plots ANALYSIS BEGINS'
LC_ALL=C python2 ${CheckMDir}/checkm dist_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 8 --gc_window_size 5000 --td_window_size 5000 --cd_window_size 5000 ${TaxFolder} ${BinFolder} ${PlotFolder} ${TetraFile} 50
echo '!!!!!!!!!!!!!!!!!!!!!GC, CD, and TD distribution plots ANALYSIS COMPLETE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo ''

echo 
echo '		PLOTS ANALYSIS BEGINS'
echo 
echo 
# ATTENTION!!!! NOT WORKING!!!
# Tetra pca needs matplotlib 1.3.1 to run!!! Otherwise it doesn't run at all!!!!!!!!!!!!
LC_ALL=C python2 ${CheckMDir}/checkm tetra_pca --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 6.5 ${BinFolder} ${PlotFolder} ${TetraFile}

echo 
echo 

# ATTENTION!!!! NOT WORKING!!!
# Tetra plot needs matplotlib 1.3.1 to run!!! Otherwise it doesn't run at all!!!!!!!!!!!!

LC_ALL=C python2 ${CheckMDir}/checkm tetra_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --td_window_size 1000 --td_bin_width 0.01 ${TaxFolder} ${BinFolder} ${PlotFolder} ${TetraFile} 50

echo 
echo 


# ATTENTION!!
# NOT WORKING!! Probably needs matplotlib 1.3.1 to run!!!!!!
LC_ALL=C python2 ${CheckMDir}/checkm gc_bias_plot --image_type pdf --dpi 600 --font_size 10 --width 6.5 --height 3.5 --window_size 1000 --all_reads --threads $SLURM_NTASKS ${BinFolder} ${PlotFolder} ${Bam_PE_File}

echo 
echo '!!!!!!!!!!!!!!!!PLOTS ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
echo 
echo 
echo '		MAPPED READS % ANALYSIS BEGINS'

LC_ALL=C python2 ${CheckMDir}/checkm profile --file ${StrXPilon}_mapped_reads_percentage.tab --tab_table ${CovFile} 

echo '!!!!!!!!!!!!!!!!!!!!MAPPED READS % ANALYSIS COMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!'


rm -r ${TempFolder} *.bam* *fasta.*

echo "==============================="
echo "PBS job finished " `date`
echo

# finished commands

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: `printf '%dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60))`

exit 0
