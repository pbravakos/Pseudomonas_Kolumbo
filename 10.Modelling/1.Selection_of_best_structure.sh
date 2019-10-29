#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
# #SBATCH --mem=128000
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Modeller_compare"
#SBATCH --output=Modeller_compare_job_%j.out

# for calculating the amount of time the job takes
begin=$(date +%s)

## Some job specific info
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

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a


<< ////
       
	This script runs in the Modeller directory. Check the tutorial https://salilab.org/modeller/tutorial/basic.html

	To run this script, first CHANGE the Initial Parameters according to needs,and then run the script: 
	sbatch $0


	NOTE:
	We run the python scripts in python3, because we have set up Modeller in our system (PYTHONPΑΤΗ) to run only in python3, but they could also run in python2.

////


# Here we will try to export from the anvio Databases the genes of interest and their annotations.

# INITIAL PARAMETERS
StrainX=Strain23
Gene=3096
StrX=${StrainX/ain/}

anviGetGeneCalls="anvi-get-sequences-for-gene-calls"

echo
echo "			ANLYSIS BEGINS!!!!"
echo


if [[ ! -e pdb_95.pir ]]; then
    wget https://salilab.org/modeller/downloads/pdb_95.pir.gz
    gunzip pdb_95.pir.gz
fi

if [[ ! -d ${StrX}_${Gene} ]]; then
    mkdir ${StrX}_${Gene}
fi

cd ${StrX}_${Gene}

ln -s ../pdb_95.pir

$anviGetGeneCalls -c ${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}/${StrX}_Pilon_contigs.db \
                  --gene-caller-ids ${Gene} \
                  -o ${HOME}/Pseudomonas/Modeller/${StrX}_${Gene}/${StrX}_gene_${Gene}.fasta \
                  --report-extended-deflines  \
                  --get-aa-sequences

# Next, we want to create the pir file needed as the Modeller input
# First we remove all newlines from the fasta file.
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' ${StrX}_gene_${Gene}.fasta > ${StrX}_gene_${Gene}_no_newlines.fasta
# Finally, we add a new line bellow the fasta header (second line) with stuctural information and then we add an asterisk (*) at the end of the pir file.
sed "s/>.*/>P1;${StrX}_${Gene}/;2isequence:${StrX}_${Gene}:::::::0.00: 0.00" ${StrX}_gene_${Gene}_no_newlines.fasta | sed -E 's/([A-Z]{4}$)/\1\*/' > ${StrX}_gene_${Gene}.pir


# Next, we export all the functions from the contigs database
anvi-export-functions -c ${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX}/${StrX}_Pilon_contigs.db \
                      -o ${HOME}/Pseudomonas/Modeller/${StrX}_${Gene}/${StrX}.annot
# Finally, we export only the function of the gene of interest in a new file. This file is not going to be used as an input in downstream analyses, it only useful for informational purposes.
awk -v gen="$Gene" 'BEGIN{FS="\t";OFS="\t"} gen==$1 {print}' ${StrX}.annot | sort -t$'\t' -k2 > ${StrX}_gene_${Gene}.annot


# Next, we create the build_profile.py python3 file.
cat > build_profile.py <<EOF
from modeller import *

log.verbose()
env = environ()

#-- Prepare the input files

#-- Read in the sequence database
sdb = sequence_db(env)
sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
          chains_list='ALL')

#-- Now, read in the binary database
sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
         chains_list='ALL')

#-- Read in the target sequence/alignment
aln = alignment(env)
aln.append(file='${StrX}_gene_${Gene}.pir', alignment_format='PIR', align_codes='ALL')

#-- Convert the input sequence/alignment into profile format
prf = aln.to_profile()
EOF

cat >> build_profile.py <<"EOF"

#-- Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

#-- Write out the profile in text format
prf.write(file='build_profile.prf', profile_format='TEXT')

#-- Convert the profile back to alignment format
aln = prf.to_alignment()

#-- Write out the alignment file
aln.write(file='build_profile.ali', alignment_format='PIR')
EOF


# Next, run the python script to create the profiles
python3 build_profile.py > build_profile.log


# Next, filter from the output of the previous python script only the rows with interest and print only the column with the pdb IDs. Finally, separate the chain from the pdb ID by a colon. 
awk 'NF > 10 && $12=="0.0" && $10 > 100 && $11 > 25 {print $2}' build_profile.prf | sed -E 's/(^.{4})/\1:/g' > pdb_IDs.txt

# awk 'BEGIN { ORS = " " } { print }' 

# Next, create a new file containing only the pdb IDs, without the chain information.
cut -d: -f1 pdb_IDs.txt | uniq > pdb_IDs_without_Chains.txt

# Next, we download all pdb files that are of interest, which are going to be used as input in the upcoming compare.py python script.
while read -r pdb ; do
    wget https://files.rcsb.org/download/${pdb}.pdb
done < pdb_IDs_without_Chains.txt

# Next, we create a variable, with all the pdb IDs and relevant chains, that can be used in the upcoming compare.py python script.
totalPDB=`wc -l pdb_IDs.txt | cut -d" " -f1`
num=0
pdb_compare=''
while read -r line ; do
    ((num++))
    addPar1=("`echo $line | cut -d: -f1`")
    addPar2=("`echo $line | cut -d: -f2`")
    addPar3="'"${addPar1}"'"", ""'"$addPar2"'"
    pdb_compare=$pdb_compare"("$addPar3")"
    # We don't want to have a comma at the end, so we rule this out.
    if [[ num -lt $totalPDB ]]; then
        pdb_compare+=", "
    fi
done < pdb_IDs.txt

# Next, create the compare.py python script using the variable $pdb_compare created in the previous while loop.
cat > compare.py <<EOF
from modeller import *

env = environ()
# We create an (initially empty) alignment object 'aln'
aln = alignment(env)
# use for loop to instruct MODELLER to read each of the PDB files
for (pdb, chain) in (${pdb_compare}):
    #  We use the model_segment argument to ask only for a single chain to be read from each PDB file (by default, all chains are read from the file).
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    # We use the append_model method to add the structure to the alignment
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
# At the end of the loop, all of the structures are in the alignment, but they are not ideally aligned to each other (append_model creates a simple 1:1 alignment with no gaps). Therefore, we improve this alignment by using malign to calculate a multiple sequence alignment.
aln.malign()
# malign3d command then performs an iterative least-squares superposition of the six 3D structures, using the multiple sequence alignment as its starting point
aln.malign3d()
# Compare_structures command compares the structures according to the alignment constructed by malign3d. It does not make an alignment, but it calculates the RMS and DRMS deviations between atomic positions and distances, differences between the mainchain and sidechain dihedral angles, percentage sequence identities, and several other measures
aln.compare_structures()
# Finally, the id_table command writes a file with pairwise sequence distances that can be used directly as the input to the dendrogram command (or the clustering programs in the PHYLIP package)
aln.id_table(matrix_file='family.mat')
# dendrogram calculates a clustering tree from the input matrix of pairwise distances, which helps visualizing differences among the template candidates.
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
EOF


# Finally, run the compare.py in python3.
python3 compare.py > compare.log





<< ////
       
	MANUAL STEP NEEDED!
	Next, we have to select manually one of the structures compared by the compare.py and use it for further analysis. In order to do this we select based on the dendrogram found in the compare.log. From the same tree we can select based on the crystallographic resolution (in angstrom, the lower the better). We also take into considaration the sequence identity from the build_profile.prf file we created earlier (the higher the better). Finally, we can choose based on the crystallographic R-factor (R-Value Work) found on each pdb file or the website https://www.rcsb.org/ (the lower the the better). We could also use the functional annotation information we already have as an additional indication of which which structure to select.

In the file build_profile.prf the important columns are the following:
	Column 2nd: reports the code of the PDB sequence that was compared with the target sequence. The PDB code in each line is the representative of a group of PDB sequences that share 95% or more sequence identity to each other and have less than 30 residues or 30% sequence length difference.
	Column 10th: reports the lengths of the alignment_format
	Column 11th: reports the percentage sequence identities between query and a PDB sequence normalized by the lengths of the alignment
	Column 12th: reports the e-value of the alignment.
	
////

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




