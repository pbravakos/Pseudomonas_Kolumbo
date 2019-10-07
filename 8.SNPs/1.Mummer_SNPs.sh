#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6400
#SBATCH --job-name="Mummer_SNPs"
#SBATCH --output=Mummer_SNPs_job_%j.out

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


generalInfo () {
    cat << ////
       
	This script runs in the mummer directory.

	To run this script two arguments (two Strains) must be given. The first argument will be used as the Reference strain and the second as the Query strain e.g.:
	sbatch $0 Strain16 Strain20

	IMPORTANT
	We want to find a reliable set of SNPs between two highly similar multi-FastA sequence sets ref.fasta and qry.fasta. We use mummer to accomplish it. Moreover, we want to create a file that correlates Anvio gene caller IDs, function calls and the SNPs.
        We can select which annotation sources we will extract, and probably we should check this selection before running this script.
        
        IMPORTANT
        Strains searched for SNPs should have very high similarity (ANI values).
	
////
}

# Check that an argument has been given in the correct form.
if [[ $# -ne 2 ]]; then
    echo "Exactly two arguments should be given!" >&2
    generalInfo >&2
    exit 1
fi

if [[ ! $1 =~ ^Strain[0-9]{2}$ ]] || [[ ! $2 =~ ^Strain[0-9]{2}$ ]] ; then
    echo "Exactly two arguments should be given, which has to be in the form of 'StrainXX' where X is a number from 0 up to 9'!" >&2
    generalInfo >&2
    exit 1
fi


# INITIAL PARAMETERS
StrainX1=$1
StrX1=${StrainX1/ain/}
StrainX2=$2
StrX2=${StrainX2/ain/}

# Available functional annotation sources .....................: Gene3D, Kegg Pathways + Enzymes, PfamInterPro, Prokka, TMHMM2.0, ProSitePatterns, KEGG Genes, Hamap, COG_FUNCTION, Rast FigFam, PANTHER GO-Slim Cellular, MobiDBLite, CDD, SMART, LipoP, Blast2Go, SignalP-4.1, PANTHER GO-Slim Function, KEGG Function, SFLD, PANTHER Protein, EggnogKegg, COG_CATEGORY, KEGG Brite/Pathway, PRINTS, PIRSF, EGGNOG, ProSiteProfiles, MetaCyc, Interpro, Phobius, TIGRFAM, ProDom, smCOG, SUPERFAMILY, Reactome, Interpo GoTerms, PANTHER GO-Slim Process
AnnotSource1=Prokka

AnnotSource2=PfamInterPro

AnnotSource3="Rast FigFam"

AnnotSource4="KEGG Genes"

ContigPilonDir1=${HOME}/Pseudomonas/Seqtk/${StrainX1}
ContigPilonDir2=${HOME}/Pseudomonas/Seqtk/${StrainX2}

anviGetSeq="anvi-get-sequences-for-gene-calls"


#------------------------------------------------------------------------------------------------------------------------------------
echo
echo " 					START THE ANALYSIS!!!!!"
echo
echo

# Create the directory where the analysis will take place
if [[ ! -d ${StrX1}_vs_${StrX2} ]]; then
    mkdir ${StrX1}_vs_${StrX2}
fi

# Start a job dependency to move the Sdtout file to the correct folders
sbatch --dependency=afterany:$SLURM_JOB_ID Move_Slurm_output_files.sh $SLURM_JOB_NAME $SLURM_JOB_ID ${StrX1}_vs_${StrX2}

# Change directory
cd ${StrX1}_vs_${StrX2}


if [[ ! -e ${StrX1}_Pilon_CLA_Blast.fasta ]]; then
    ln -s ${ContigPilonDir1}/${StrX1}_Pilon_CLA_Blast.fasta
fi

if [[ ! -e ${StrX2}_Pilon_CLA_Blast.fasta ]]; then
    ln -s ${ContigPilonDir2}/${StrX2}_Pilon_CLA_Blast.fasta
fi

echo
echo "Nucmer version used for this analysis is: `nucmer -V`"
echo
echo

echo
echo "			${StrX1} vs ${StrX2} calculation of SNPs begins"
echo
# First, run nucmer to find the differences between the two fasta files and store these differences in a delta file.
nucmer --prefix=${StrX1}_${StrX2} ${StrX1}_Pilon_CLA_Blast.fasta ${StrX2}_Pilon_CLA_Blast.fasta

# Then, use the delta file created previously as an input, to the next command.
# The -C option in show-snps assures that only SNPs found in uniquely aligned sequence will be reported, thus excluding SNPs contained in repeats.
show-snps -Clr ${StrX1}_${StrX2}.delta > ${StrX1}_${StrX2}.snps
 
# An alternative method which first attempts to determine the "correct" repeat copy is:
# Conflicting repeat copies will first be eliminated with delta-filter and the SNPs will be re-called in hopes of finding some that were previously masked by another repeat copy.
delta-filter -r -q ${StrX1}_${StrX2}.delta > ${StrX1}_${StrX2}.filter

show-snps -CTlr ${StrX1}_${StrX2}.filter > ${StrX1}_${StrX2}_filtered.snps

echo
echo "--------------------------------------------SNPs analysis with MUMMER was completed--------------------------------------------------------------------"
echo
echo
echo
echo

# Create a link for the Anvio contigs database (both Strains).
if [[ ! -e ${StrX1}_Pilon_contigs.db ]]; then
    ln -s ${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX1}/${StrX1}_Pilon_contigs.db
fi

if [[ ! -e ${StrX2}_Pilon_contigs.db ]]; then
    ln -s ${HOME}/Pseudomonas/Anvio/Strains/Prokka/${StrainX2}/${StrX2}_Pilon_contigs.db
fi

# Export information from the contig database for the Reference strain that is going to prove useful later on, and save to a new file.
anvi-export-table --table genes_in_contigs --fields 'gene_callers_id, contig, start, stop, direction' --output-file ${StrX1}_Anvio_gene_callers ${StrX1}_Pilon_contigs.db

# Repeat for the query strain. For the query strain the goal is to have in one file the directionality of the SNP gene and the gene caller ids of the Reference strain.
anvi-export-table --table genes_in_contigs --fields 'gene_callers_id, contig, start, stop, direction' --output-file ${StrX2}_Anvio_gene_callers ${StrX2}_Pilon_contigs.db

# Remove first line (header) of the newly created gene callers file, and rearange the columns order, so as to have the Scaffold column first.
sed '1d' ${StrX1}_Anvio_gene_callers | awk -F'[\t]' '{print $2"\t"$1"\t"$3"\t"$4"\t"$5}' > ${StrX1}_Anvio_gene_callers_inverse

# Repeat for the Query strain.
sed '1d' ${StrX2}_Anvio_gene_callers | awk -F'[\t]' '{print $2"\t"$1"\t"$3"\t"$4"\t"$5}' > ${StrX2}_Anvio_gene_callers_inverse


# Next, we will manipulate the Mummer SNPs output, in order to have a space separated file, with only informative columns i.e. gene caller IDs and SNP positions.
# Remove the first 4 lines (headers).
sed '1,4d' ${StrX1}_${StrX2}_filtered.snps > ${StrX1}_${StrX2}_fixed.snps

# Extract the data (columns) of the first (reference) strain from the snps mummer output, and create a new file with this data.
awk '{print $11"\t"$1"\t"$2"\t"$3}' ${StrX1}_${StrX2}_fixed.snps > ${StrX1}.snps

# Extract the data of the second (query) strain from the snps mummer output, and create a new file with this data.
awk '{print $12"\t"$4"\t"$3"\t"$11"\t"$1}' ${StrX1}_${StrX2}_fixed.snps > ${StrX2}.snps

{
# Next, we will try to create a new file that shows all the genes from the reference strain that have a snp, along with all the relevant information.
# IMPORTANT
# Any SNPs that cannot be found in a coding region, are excluded from the analysis!

# Fist, we check if there is a file with the same name that we are going to create in the while loop in the next step. If there is, we delete it because inside the while loop we append to this file, and things can go terribly wrong if we already have some information in this file before going to the while loop!
OuputWhileLoop1=deleteMe
if [[ -e $OuputWhileLoop1 ]]; then
    rm $OuputWhileLoop1
fi

# Next, we remove any spaces from the file that is going to be read in the while loop, because by default the read command separates columns by space. We have to be careful though because this change (substitution of spaces with underscores) may change the columns we are comparing in the awk command. Happily, in our case columns compared in the awk command DO NOT have any spaces, so remain unaffected by the change! 
sed -i 's/ /_/g' $OuputWhileLoop1

# Now the while loop which will create the file we want for the Reference strain.
while read -r line ; do
    echo $line > line1
    Scaffold=`cut -d' ' -f1 line1`
    SNP_pos=`cut -d' ' -f2 line1`
    SNP_Nucl=`cut -d' ' -f3 line1`
    Nucl2=`cut -d' ' -f4 line1`
    rm line1
    awk -v scaf="$Scaffold" -v pos="$SNP_pos" -v nucl="$SNP_Nucl" -v nuc2="$Nucl2" 'BEGIN{FS="\t";OFS="\t"} (scaf==$1 && pos>=$3 && pos<=$4) {print $2,$1,$3,$4,$5,pos,nucl,nuc2}' ${StrX1}_Anvio_gene_callers_inverse >> $OuputWhileLoop1
done < ${StrX1}.snps
}


{

# Repeat for Query strain
Ouput2WhileLoop1=deleteMe2
if [[ -e $Ouput2WhileLoop1 ]]; then
    rm $Ouput2WhileLoop1
fi

# We substitute any spaces with underscores for the file that is going to be read in with the read command.
sed -i 's/ /_/g' ${StrX2}.snps

# IMPORTANT
# We may loose some gene caller ids (rows) in the following while loop, if there is a snp in the query strain scaffold that does not correspond to a coding region of the query strain!

# Only difference in the query strain are the columns we have as output from the while loop. In reality we do not need all of them! This is why later on we will continue working with a short version of the $Ouput2WhileLoop1.
while read -r line ; do
    echo $line > line1
    Scaffold=`cut -d' ' -f1 line1`
    SNP_pos=`cut -d' ' -f2 line1`
    SNP_Nucl=`cut -d' ' -f3 line1`
    ScaffRef=`cut -d' ' -f4 line1`
    SnpPosRef=`cut -d' ' -f5 line1`
    rm line1
    awk -v scaf="$Scaffold" -v pos="$SNP_pos" -v nucl="$SNP_Nucl" -v scafRef="$ScaffRef" -v posRef="$SnpPosRef" 'BEGIN{FS="\t";OFS="\t"} (scaf==$1 && pos>=$3 && pos<=$4) {print $2,$1,$3,$4,$5,pos,nucl,scafRef,posRef}' ${StrX2}_Anvio_gene_callers_inverse >> $Ouput2WhileLoop1
done < ${StrX2}.snps


# We add a header to the file we created it previously inside the while loop, and rename it.
awk -v str1="${StrX1}" -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","Scaffold_"str1,"Start","End","Direction","SNP_position_"str1,"SNP_Reference_"str1,"SNP_Query_"str2} {print}' $OuputWhileLoop1 > ${StrX1}_Anvio_gene_callers_SNPs.tsv

# We repeat for the query strain, adding the correct headers.
awk -v str1="${StrX1}" -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id_"str2,"Scaffold_"str2,"Start_"str2,"End_"str2,"Direction_"str2,"SNP_position_"str2,"SNP_Query_"str2,"Scaffold_"str1,"SNP_position_"str1} {print}' $Ouput2WhileLoop1 > ${StrX2}_Anvio_gene_callers_SNPs.tsv


# Finally, we can safely delete the $OuputWhileLoop1 and $Ouput2WhileLoop1  files (because we renamed them and saved them in the previous step!)
rm $OuputWhileLoop1 $Ouput2WhileLoop1
}

{
# We repeat (more or less) the same procedure as above, to append the functional annotation to our file, but only for the Reference strain!

# First, we export the functions from the Anvio Database 
# In order to check the available sources in the Anvio Database we can give the following command:
# anvi-export-functions -c ${StrX1}_Pilon_contigs.db --output-file dokimi --list-annotation-sources
# The selection of sources we extract can change, in future iterations of this script.
anvi-export-functions -c ${StrX1}_Pilon_contigs.db --output-file ${StrX1}_Prokka_Pfam.func --annotation-sources "$AnnotSource1","$AnnotSource2","$AnnotSource3","$AnnotSource4"

# Again, we don't want to have any file with the same name as the output of the upcoming while loop in our directory, so we check for it.
OuputWhileLoop2=${StrX1}_versus_${StrX2}_Anvio_gene_callers_function_SNPs.tsv
if [[ -e $OuputWhileLoop2 ]]; then
    rm $OuputWhileLoop2
fi

# We substitute any spaces with underscores for the file that is going to be read with the read command.
sed -i 's/ /_/g' ${StrX1}_Anvio_gene_callers_SNPs.tsv

# IMPORTANT
# Not all genes have functional annotation. For this reason the file may be missing some or many gene caller ids! As a result in the next while loop it is possible to filter out many of our data, without having the intention to do so!

# In the next while loop we will append to "${StrX1}_Anvio_gene_callers_SNPs.tsv" the functional annotation and name it as $OuputWhileLoop2
while read -r line ; do
    echo $line > line1
    GeneID=`cut -d' ' -f1 line1`
    Scaffold=`cut -d' ' -f2 line1`
    Start=`cut -d' ' -f3 line1`
    End=`cut -d' ' -f4 line1`
    Direct=`cut -d' ' -f5 line1`
    Posit=`cut -d' ' -f6 line1`
    Nucleo=`cut -d' ' -f7 line1`
    Nucleo2=`cut -d' ' -f8 line1`
    rm line1
    awk -v GenID="$GeneID" -v scaf="$Scaffold" -v star="$Start" -v telos="$End" -v kateuth="$Direct" -v pos="$Posit" -v nuc="$Nucleo" -v nuc2="$Nucleo2" 'BEGIN{FS="\t";OFS="\t"} GenID==$1 {print $1,$2,$3,$4,scaf,star,telos,pos,kateuth,nuc,nuc2}' ${StrX1}_Prokka_Pfam.func >> $OuputWhileLoop2
done < ${StrX1}_Anvio_gene_callers_SNPs.tsv # Good part is that we also get a nice (and correct) header for $OuputWhileLoop2 from the above command!



# Next we will try to append to the $OuputWhileLoop2 file the missing (if there are any!) gene callers, that did not have any functional information, and for this reason were excluded from the analysis. 

# First we extract all the gene caller ids that were included in the final ouput.
sort -u $OuputWhileLoop2 | awk '{print $1}' > Reference_Genes_with_function_calls

# Second, we extract all the gene caller ids which should be in the final output.
sort -u ${StrX1}_Anvio_gene_callers_SNPs.tsv | awk '{print $1}' > All_Reference_set_genes

# Next, we compare the two files created earlier, and create a new one with the gene callers ids that exist in one but not in the other file.
awk 'FNR==NR {a[$0]++; next} !a[$0]' Reference_Genes_with_function_calls All_Reference_set_genes | uniq > ${StrX1}_genes_with_SNPs_removed_from_analysis

# Finally, if the file created earlier is not empty, we append to it all the missing information (except the functional annotation!).
if [[ -s ${StrX1}_genes_with_SNPs_removed_from_analysis ]]; then

    awk 'NR==FNR{c[$1]++;next};c[$1]' ${StrX1}_genes_with_SNPs_removed_from_analysis ${StrX1}_Anvio_gene_callers_SNPs.tsv > ${StrX1}_Anvio_gene_callers_SNPs_no_functions.tsv
    
    awk -v str1="${StrX1}" -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t"} {print $1,"none","none","none",$2,$3,$4,$6,$5,$7,$8 >> str1"_versus_"str2"_Anvio_gene_callers_function_SNPs.tsv"}' ${StrX1}_Anvio_gene_callers_SNPs_no_functions.tsv

fi
}


{
# Next, we create an awk file to help us find the complement nucleotides for all reverse strand SNPs in the Query strain.
cat > awkfileQuerySnps <<"EOF"
BEGIN{FS="\t";OFS="\t"}
{
    if ($5=="r" && $7=="A") 
    {
      print $1,$2,$3,$4,$5,$6,"T",$8,$9
    }
    else if ($5=="r" && $7=="T") 
    { 
      print $1,$2,$3,$4,$5,$6,"A",$8,$9
    }
    else if ($5=="r" && $7=="C") 
    { 
      print $1,$2,$3,$4,$5,$6,"G",$8,$9
    }
    else if ($5=="r" && $7=="G") 
    { 
      print $1,$2,$3,$4,$5,$6,"C",$8,$9
    }
   else
   {
     print $0
   }
}
EOF

# Finally, we will print a file where there is a column that indicates, the exact position of the SNP in each Gene. 
# IMPORTANT!
# When the gene is in reverse order ( $9 == r) then the actual SNP we expect to find in the gene is the complement of what we see i.e. if we have a "C" then we expect to find in the gene sequence a "G". This is why we use the awk file to turn all the reverse strand SNPs to their complements. The next awk help us calculate the exact SNP position in each gene (i.e. the first nucleotide of each gene is always no 1).
awk -F"\t" -FS"\t" -f awkfileQuerySnps ${StrX2}_Anvio_gene_callers_SNPs.tsv | awk -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t"} NR==1{print $1,$2,$6,"SNP_position_in_gene_"str2,$7,$8,$9} NR>=2 && $5=="f"{print $1,$2,$6,$6-$3,$7,$8,$9} NR>=2 && $5=="r"{print $1,$2,$6,$4-$6+1,$7,$8,$9}' > ${StrX2}_intermediate1
}


{

# Check if the output we append to in the upcoming while loop, exists or not. If yes then delete it.
if [[ -e ${StrX2}_intermediate1_loop ]]; then
    rm ${StrX2}_intermediate1_loop
fi

# First we make the rows unique based on the  columns with the Scaffold Reference name AND the Reference Snp position in the Scaffold.
awk 'BEGIN{FS="\t";OFS="\t"} !seen[$5,$8]++ {print}' $OuputWhileLoop2 > ${StrX1}_vs_${StrX2}_temp1

# Once again we want to make sure that the input of the read command has no spaces in the columns. We also check if this space substitution affects the equality if statement (==) inside the awk command.
sed -i 's/ /_/g' ${StrX1}_vs_${StrX2}_temp1


# IMPORTANT: 
# In the next awk command, we may loose any genes that are duplicates (because were found by different programs) e.g. have different gene caller ids and different start and/or end position in the scaffold but the same SNP position in the scaffold!
# IMPORTANT:
# In the next awk command we may also loose any Reference genes that cannot be found in the Query strain in coding regions. This happens here because it is the first time we merge the information we have from the Refernce and the Query strains together into one file.

# We want to have in one file the gene callers ID from the query strain and the gene callers Id from the Reference strain.
while read -r line ; do
    echo $line > line1
    GenIDRef=`cut -d' ' -f1 line1`
    ScaffRef=`cut -d' ' -f5 line1`
    SnpPosRef=`cut -d' ' -f8 line1`
    SnpNucRef=`cut -d' ' -f10 line1`
    rm line1
    awk -v genID="$GenIDRef" -v scafR="$ScaffRef" -v snpR="$SnpPosRef" -v snpNucR="$SnpNucRef" 'BEGIN{FS="\t";OFS="\t"} scafR==$6 && snpR==$7 {print $0,genID,snpNucR}' ${StrX2}_intermediate1 >> ${StrX2}_intermediate1_loop
done < ${StrX1}_vs_${StrX2}_temp1

}


{
# Next we create the DNA codon table, which we will use in the upcoming for loops, to translate the Query and Reference codons into amino acids.
cat > DNA_codon_table  <<"EOF"
s/ACT/Thr/g
s/ACC/Thr/g
s/ACA/Thr/g
s/ACG/Thr/g
s/GCT/Ala/g
s/GCC/Ala/g
s/GCA/Ala/g
s/GCG/Ala/g
s/TAT/Tyr/g
s/TAC/Tyr/g
s/TAA/Stop/g
s/TAG/Stop/g
s/CAT/His/g
s/CAC/His/g
s/CAA/Gln/g
s/CAG/Gln/g
s/AAT/Asn/g
s/AAC/Asn/g
s/AAA/Lys/g
s/AAG/Lys/g
s/GAT/Asp/g
s/GAC/Asp/g
s/GAA/Glu/g
s/GAG/Glu/g
s/TGT/Cys/g
s/TGC/Cys/g
s/TGA/Stop/g
s/TGG/Trp/g
s/CGT/Arg/g
s/CGC/Arg/g
s/CGA/Arg/g
s/CGG/Arg/g
s/AGA/Arg/g
s/AGG/Arg/g
s/TCT/Ser/g
s/TCC/Ser/g
s/TCA/Ser/g
s/TCG/Ser/g
s/AGT/Ser/g
s/AGC/Ser/g
s/GGT/Gly/g
s/GGC/Gly/g
s/GGA/Gly/g
s/GGG/Gly/g
s/TTT/Phe/g
s/TTC/Phe/g
s/TTA/Leu/g
s/TTG/Leu/g
s/CTT/Leu/g
s/CTC/Leu/g
s/CTA/Leu/g
s/CTG/Leu/g
s/ATT/Ile/g
s/ATC/Ile/g
s/ATA/Ile/g
s/ATG/Met/g
s/GTT/Val/g
s/GTC/Val/g
s/GTA/Val/g
s/GTG/Val/g
s/CCT/Pro/g
s/CCC/Pro/g
s/CCA/Pro/g
s/CCG/Pro/g
EOF
}

{
# Next, we prepare for the upcoming for loop.
numGenCalls=`sort -nu  ${StrX2}_intermediate1_loop | sed '1d' | awk -F"\t" '{print $1}' | wc -l` # Get the number of genes
sort -nu  ${StrX2}_intermediate1_loop | sed '1d' | awk -F"\t" '{print $1}' | awk 'BEGIN { ORS = " " } { print }'> ${StrX2}_GeneCalls.txt # We save in one line the actuall genes, separated by space.

# check before going to the for loop.
if [[ -e ${StrX2}_intermediate2_loop ]]; then 
    rm ${StrX2}_intermediate2_loop 
fi

# Next we will try to find all the Query codons that have SNPs, and combine them in one file with gene caller IDs, Snp posistion in gene and Snp position in the scaffold. 

for i in $(seq 1 $numGenCalls); do
    GeneName=`cut -d' ' -f${i} ${StrX2}_GeneCalls.txt`
    $anviGetSeq -c ${StrX2}_Pilon_contigs.db --gene-caller-ids $GeneName -o ${GeneName}_gene_call.fasta
    awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' ${GeneName}_gene_call.fasta | sed '1d' > ${GeneName}_gene  
    wc -L ${GeneName}_gene > wcMaxLW
    read geneLen filename < wcMaxLW
    genLen=$geneLen
#    if [[ -e Codon_start_num ]]; then 
#       rm Codon_start_num
#    fi
    for e in $(seq 1 3 $geneLen); do echo $e>> Codon_start_num; done
    awk '{print $1,$1+2}' Codon_start_num > Codon_intervals
    SameGeneSNPsNum=`awk '{print $1}' ${StrX2}_intermediate1_loop | grep -c "$GeneName"`
    if [[ $SameGeneSNPsNum -gt 1 ]]; then
        for u in $(seq 1 $SameGeneSNPsNum); do
            SnpPos=`awk -v genName=$GeneName 'genName==$1 {print}' ${StrX2}_intermediate1_loop | awk -v nu=$u 'NR==nu{print $4}'`
            awk -v GeneLen="$geneLen" -v snpP="$SnpPos" 'GeneLen>=$2 && snpP>=$1 && snpP<=$2 {print}' Codon_intervals > ${GeneName}_${SnpPos}_codon_range
            StartPos=`awk '{print $1}' ${GeneName}_${SnpPos}_codon_range`
            EndPos=`awk '{print $2}' ${GeneName}_${SnpPos}_codon_range`
            cut -c ${StartPos}-${EndPos} ${GeneName}_gene > ${GeneName}_${SnpPos}_codon_seq
            echo ${GeneName} > geneName.txt
            echo ${SnpPos} > SnpPos.txt
            sed -f DNA_codon_table ${GeneName}_${SnpPos}_codon_seq > Ref_aminoAcids
            awk -v genName=$GeneName 'genName==$1 {print}' ${StrX2}_intermediate1_loop | awk -v nu=$u 'NR==nu{print $3}' > ${GeneName}_${SnpPos}_SnpPosScaff
            paste geneName.txt SnpPos.txt Ref_aminoAcids ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_SnpPosScaff >> ${StrX2}_intermediate2_loop         
            rm ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_codon_range ${GeneName}_${SnpPos}_SnpPosScaff
        done
    else
        SnpPos=`awk -v geneN="$GeneName" 'geneN==$1 {print $4}' ${StrX2}_intermediate1_loop` 
        awk -v GeneLen="$geneLen" -v snpP="$SnpPos" 'GeneLen>=$2 && snpP>=$1 && snpP<=$2 {print}' Codon_intervals > ${GeneName}_${SnpPos}_codon_range
        StartPos=`awk '{print $1}' ${GeneName}_${SnpPos}_codon_range`
        EndPos=`awk '{print $2}' ${GeneName}_${SnpPos}_codon_range`
        cut -c ${StartPos}-${EndPos} ${GeneName}_gene > ${GeneName}_${SnpPos}_codon_seq
        echo ${GeneName} > geneName.txt
        echo ${SnpPos} > SnpPos.txt
        sed -f DNA_codon_table ${GeneName}_${SnpPos}_codon_seq > Ref_aminoAcids
        awk -v geneN="$GeneName" 'geneN==$1 {print $3}' ${StrX2}_intermediate1_loop > ${GeneName}_${SnpPos}_SnpPosScaff
        paste geneName.txt SnpPos.txt Ref_aminoAcids ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_SnpPosScaff >> ${StrX2}_intermediate2_loop        
        rm ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_codon_range ${GeneName}_${SnpPos}_SnpPosScaff
    fi
     
    rm Codon_start_num # Do not comment this line! NEVER!!!!
    rm ${GeneName}_gene_call.fasta ${GeneName}_gene wcMaxLW Codon_intervals geneName.txt SnpPos.txt
done

# Next, we add a header to the $OutputForLoop we just created. The last column is the snp position in the scaffold
awk -v str1="${StrX1}" -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id_"str2,"SNP_position_in_gene_"str2,"AA_Query_"str2,"Codon_Query_"str2,"SNP_position_"str2} {print}' ${StrX2}_intermediate2_loop > ${StrX2}_intermediate2_loop_head
}



{
# Next, prepare for the upcoming while loop.
if [[ -e ${StrX2}_intermediate3_loop ]]; then
    rm ${StrX2}_intermediate3_loop
fi

# Again, we substitute any spaces in the input file to the read command. We take notice that this substitution does not affect the if equality statement (==) in the awk command.
sed -i 's/ /_/g' ${StrX2}_intermediate2_loop_head

# We want to add the amino acid and codon information for the query strain, we got earlier to our file.
while read -r line ; do
    echo $line > line1
    GenIDQuer=`cut -d' ' -f1 line1`
    SnpGeneQuer=`cut -d' ' -f2 line1`
    AAQuer=`cut -d' ' -f3 line1`
    CodonQuer=`cut -d' ' -f4 line1`
    SnpScaffQuer=`cut -d' ' -f5 line1`
    rm line1
    awk -v genID="$GenIDQuer" -v snpGQ="$SnpGeneQuer" -v AAQ="$AAQuer" -v cdnQ="$CodonQuer" -v snpScQ="$SnpScaffQuer" 'BEGIN{FS="\t";OFS="\t"} genID==$1 && snpScQ==$3 {print $0,AAQ,cdnQ}' ${StrX2}_intermediate1_loop >> ${StrX2}_intermediate3_loop
done < ${StrX2}_intermediate2_loop_head


# Next, we will try to append to the output of the previous while loop, the SNps that exist in coding regions in the Reference strains, but do not exist in the Query strains and have been removed from the analysis.
awk '{print $8}' ${StrX2}_intermediate3_loop  | sort -u > ${StrX1}_inter1

awk '{print $1}' ${StrX1}_vs_${StrX2}_temp1 | sort -u > ${StrX1}_temp1

awk 'FNR==NR {a[$0]++; next} !a[$0]' ${StrX1}_inter1 ${StrX1}_temp1 | uniq > Reference_genes_with_SNPs_removed_from_analysis

# Finally, if the file created earlier is not empty, we append to it all the missing information of the Reference strains (except the info concerning the query strains!).
if [[ -s Reference_genes_with_SNPs_removed_from_analysis ]]; then

    awk 'NR==FNR{c[$1]++;next};c[$1]' Reference_genes_with_SNPs_removed_from_analysis ${StrX1}_vs_${StrX2}_temp1  > ${StrX1}_Anvio_gene_SNPs_NOT_in_Query_strain.tsv
    
    awk -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t"} {print "0","none","0","0","non_coding",$5,$8,$1,$10,"none","none" >> str2"_intermediate3_loop"}' ${StrX1}_Anvio_gene_SNPs_NOT_in_Query_strain.tsv

fi
}


{
# Prepare for the upcoming while loop.
if [[ -e ${StrX1}_${StrX2}_combined ]]; then
    rm ${StrX1}_${StrX2}_combined
fi

# Make the rows unique, based on the gene calls IDs, Scaffold AND SNP nucleotide of the Refernce strain.
awk 'BEGIN{FS="\t";OFS="\t"} !seen[$7,$8,$9]++ {print}' ${StrX2}_intermediate3_loop > ${StrX2}_intermediate3_loop_seen

sed -i 's/ /_/g' ${StrX2}_intermediate3_loop_seen

# Here, we want to have a (huge) file which combines most of the information of Reference and Query strains.
while read -r line ; do
    echo $line > line1
    GenIDQuer=`cut -d' ' -f1 line1`
    SnpScaffQuer=`cut -d' ' -f3 line1`
    SnpGeneQuer=`cut -d' ' -f4 line1`
    SnpNucQuer=`cut -d' ' -f5 line1`
    SnpScafRef=`cut -d' ' -f7 line1`
    GenIDRef=`cut -d' ' -f8 line1`
    SnpNucRef=`cut -d' ' -f9 line1`
    AAQuer=`cut -d' ' -f10 line1`
    CodonQuer=`cut -d' ' -f11 line1`
    
    rm line1
    awk -v genIDQ="$GenIDQuer" -v snpScQ="$SnpScaffQuer" -v snpPosQ="$SnpNucQuer" -v snpGQ="$SnpGeneQuer" -v snpScR="$SnpScafRef" -v genIDR="$GenIDRef" -v snpPosR="$SnpNucRef" -v AAQ="$AAQuer" -v cdnQ="$CodonQuer" 'BEGIN{FS="\t";OFS="\t"} genIDR==$1 && snpScR==$8 && snpPosR==$10 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,genIDQ,snpScQ,snpGQ,snpPosQ,AAQ,cdnQ}' $OuputWhileLoop2 >> ${StrX1}_${StrX2}_combined
done < ${StrX2}_intermediate3_loop_seen
}



{
# Next, we create an awk file, which will help us turn all the SNPs found in the reverse strand of the Reference strain to their complements.
cat > awkfileSnpReference <<"EOF"
BEGIN{FS="\t";OFS="\t"}
{
    if ($9=="r" && $10=="A") 
    {
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,"T",$14,$15,$16,$11,$13
    }
    else if ($9=="r" && $10=="T") 
    { 
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,"A",$14,$15,$16,$11,$13
    }
    else if ($9=="r" && $10=="C") 
    { 
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,"G",$14,$15,$16,$11,$13
    }
    else if ($9=="r" && $10=="G") 
    { 
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,"C",$14,$15,$16,$11,$13
    }
   else
   {
     print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,$15,$16,$11,$13
   }
}
EOF

# Next, we use the awk file to find the complements in the reverse strand, and then calculate the Snp position with start (1st position) the start of each gene that the SNP is part of.
awk -F"\t" -FS"\t" -f awkfileSnpReference ${StrX1}_${StrX2}_combined  | awk 'BEGIN{FS="\t";OFS="\t"} NR==1{print $1,$2,$3,$4,"SNP_position_in_gene",$10,$11,$12,$13,$14,$15} NR>=2 && $9=="f"{print $1,$2,$3,$4,$8-$6,$10,$11,$12,$13,$14,$15} NR>=2 && $9=="r"{print $1,$2,$3,$4,$7-$8+1,$10,$11,$12,$13,$14,$15}' > ${StrX1}_${StrX2}_before_final
}



{
# Next, we prepare for the upcoming for loop.
numGenCalls=`sort -nu  ${StrX1}_${StrX2}_before_final | sed '1d' | awk -F"\t" '{print $1}' | wc -l` # Find the number of genes.
sort -nu  ${StrX1}_${StrX2}_before_final | sed '1d' | awk -F"\t" '{print $1}' | awk 'BEGIN { ORS = " " } { print }'> GeneCalls.txt # Sort all genes in a line, separated by a space.
awk -F"\t" '!seen[$1,$5]++ {print $1,$5}' ${StrX1}_${StrX2}_before_final | sed 's/\.//g' > SNPs_unique_gene_calls.tsv  # Remove duplicates and export only gene callers IDs and SNP position in gene for the Refernce strains only.

# Check for the existance of the upcoming for loop output file. If it exists remove it!
if [[ -e ${StrX1}_intermediate_for_loop ]]; then
    rm ${StrX1}_intermediate_for_loop
fi

for i in $(seq 1 $numGenCalls); do
    GeneName=`cut -d' ' -f${i} GeneCalls.txt`
    $anviGetSeq -c ${StrX1}_Pilon_contigs.db --gene-caller-ids $GeneName -o ${GeneName}_gene_call.fasta
    awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' ${GeneName}_gene_call.fasta | sed '1d' > ${GeneName}_gene
    wc -L ${GeneName}_gene > wcMaxLW
    read geneLen filename < wcMaxLW
    genLen=$geneLen
    for e in $(seq 1 3 $geneLen); do echo $e>> Codon_start_num; done
    awk '{print $1,$1+2}' Codon_start_num > Codon_intervals
    SameGeneSNPsNum=`awk '{print $1}' SNPs_unique_gene_calls.tsv | grep -c "$GeneName"`
    
    if [[ $SameGeneSNPsNum -gt 1 ]]; then
        for u in $(seq 1 $SameGeneSNPsNum); do
            SnpPos=`awk -v genName=$GeneName 'genName==$1 {print}' SNPs_unique_gene_calls.tsv | awk -v nu=$u 'NR==nu{print $2}'`
            awk -v GeneLen="$geneLen" -v snpP="$SnpPos" 'GeneLen>=$2 && snpP>=$1 && snpP<=$2 {print}' Codon_intervals > ${GeneName}_${SnpPos}_codon_range
            StartPos=`awk '{print $1}' ${GeneName}_${SnpPos}_codon_range`
            EndPos=`awk '{print $2}' ${GeneName}_${SnpPos}_codon_range`
            cut -c ${StartPos}-${EndPos} ${GeneName}_gene > ${GeneName}_${SnpPos}_codon_seq
            echo ${GeneName} > geneName.txt
            echo ${SnpPos} > SnpPos.txt
            sed -f DNA_codon_table ${GeneName}_${SnpPos}_codon_seq > Ref_aminoAcids
            paste geneName.txt SnpPos.txt Ref_aminoAcids ${GeneName}_${SnpPos}_codon_seq >> ${StrX1}_intermediate_for_loop       
            rm ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_codon_range
        done
    else
        SnpPos=`awk -v geneN="$GeneName" 'geneN==$1 {print $2}' SNPs_unique_gene_calls.tsv` 
        awk -v GeneLen="$geneLen" -v snpP="$SnpPos" 'GeneLen>=$2 && snpP>=$1 && snpP<=$2 {print}' Codon_intervals > ${GeneName}_${SnpPos}_codon_range
        StartPos=`awk '{print $1}' ${GeneName}_${SnpPos}_codon_range`
        EndPos=`awk '{print $2}' ${GeneName}_${SnpPos}_codon_range`
        cut -c ${StartPos}-${EndPos} ${GeneName}_gene > ${GeneName}_${SnpPos}_codon_seq
        echo ${GeneName} > geneName.txt
        echo ${SnpPos} > SnpPos.txt
        sed -f DNA_codon_table ${GeneName}_${SnpPos}_codon_seq > Ref_aminoAcids
        paste geneName.txt SnpPos.txt Ref_aminoAcids ${GeneName}_${SnpPos}_codon_seq >> ${StrX1}_intermediate_for_loop         
        rm ${GeneName}_${SnpPos}_codon_seq ${GeneName}_${SnpPos}_codon_range
    fi
     
    rm Codon_start_num # Do not comment this line! NEVER!!!!
    rm ${GeneName}_gene_call.fasta ${GeneName}_gene wcMaxLW Codon_intervals geneName.txt SnpPos.txt 

done


# Next, we add a header to the $OutputForLoop we just created.
awk -v str1="${StrX1}" -v str2="${StrX2}" 'BEGIN{FS="\t";OFS="\t";print "gene_callers_id","SNP_position_in_gene","AA_Reference_"str1,"Codon_Reference_"str1} {print}' ${StrX1}_intermediate_for_loop > ${StrX1}_intermediate_for_loop_head
}


{
# Again, delete the file before going to the while loop if it exists.
OuputWhileLoop5=${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv 
if [[ -e ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv ]]; then
    rm ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv
fi

# Remove any spaces from the input of the read command.
sed -i 's/ /_/g' ${StrX1}_intermediate_for_loop_head

# Next, we create the final file that contains all the information we wanted.
while read -r line ; do
    echo $line > line1
    GenIDRef=`cut -d' ' -f1 line1`
    SnpPos=`cut -d' ' -f2 line1`
    AARef=`cut -d' ' -f3 line1`
    CdnRef=`cut -d' ' -f4 line1`
    rm line1
    awk -v genID="$GenIDRef" -v snpP="$SnpPos" -v AArefe="$AARef"  -v cdnR="$CdnRef" 'BEGIN{FS="\t";OFS="\t"} genID==$1 && snpP==$5 {print $1,$2,$3,$4,$5,$6,AArefe,cdnR,$7,$8,$9,$10,$11}' ${StrX1}_${StrX2}_before_final >> ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv
done < ${StrX1}_intermediate_for_loop_head
}

# We want to replace ASCII encodings with characters, because in certain annotations like Prokka and Rast Figfam, we have such a problem.
sed -i 's/%2C/,/g;s/%26/\&/g;s/%3B/;/g'	 ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv

{
# Next, we will try to export all the genes that have a SNP, in a new fasta file.
# First, we get the all the gene calls, that have SNPs, in one line, separated by commas.
GeneCalls=`sort -nu  ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv | sed '1d' | awk -F"\t" '{print $1}' | awk 'BEGIN { ORS = "," } { print }' | sed 's/,$//g'`

# Finally, we create a fasta file with the nucleotide sequences
$anviGetSeq -c ${StrX1}_Pilon_contigs.db --gene-caller-ids $GeneCalls -o ${StrX1}_versus_${StrX2}_gene_calls_with_SNPs.fasta
}



#--------------------------------------------------------------------------------------------------------------
# NOTE:
# There are three ways a SNP can be removed from the analysis:
# 1) If SNP cannot be found in a coding region in the Reference strain, but it can be found in a coding region in the Query strain
# 2) If SNP cannot be found in a coding region in the Query strain, but can be found in a coding region in the Reference strain
# 3) If SNP cannot be found in a coding region in the Reference strain, but the snp position in the scaffold is the exact same with other genes.
# 4) If SNP cannot be found in a coding region either in Reference or Query strains. 

# We do not mind about the 2nd and 4th option - on the opposite, it is exactly what we want to achieve- because here we examine the snps found in coding regions for the Reference strain. 
# We have taken care of the 1st option, in this script. 
# Only the 3rd option remains unresolved, but it is not a big issue because, these genes are usually duplicates 

# Next we will find which gene caller IDs (corresponding to the Reference strain) have been (if any) excluded from the analysis. We will not do any further actions on them, but it is good to know...
sort -u ${StrX1}_versus_${StrX2}_SNPs_AA_final.tsv | awk '{print $1}' > Final_Reference_genes

sort -u ${StrX1}_Anvio_gene_callers_SNPs.tsv | awk '{print $1}' > All_Reference_set_genes

awk 'FNR==NR {a[$0]++; next} !a[$0]' Final_Reference_genes All_Reference_set_genes | uniq > ${StrX1}_gene_caller_ids_with_SNPs_removed_from_analysis

# Remove files no longer needed.
rm ${StrX1}_Anvio_gene_callers ${StrX2}_Anvio_gene_callers ${StrX1}_Anvio_gene_callers_inverse ${StrX2}_Anvio_gene_callers_inverse ${StrX1}_${StrX2}_fixed.snps ${StrX1}.snps ${StrX2}.snps ${StrX1}_Anvio_gene_callers_SNPs.tsv ${StrX2}_Anvio_gene_callers_SNPs.tsv ${StrX1}_Prokka_Pfam.func ${StrX1}_Pilon_CLA_Blast.fasta ${StrX2}_Pilon_CLA_Blast.fasta ${StrX1}_Pilon_contigs.db ${StrX2}_Pilon_contigs.db ${StrX1}_${StrX2}.delta ${StrX1}_${StrX2}.filter DNA_codon_table GeneCalls.txt SNPs_unique_gene_calls.tsv Final_Reference_genes All_Reference_set_genes ${StrX1}_versus_${StrX2}_Anvio_gene_callers_function_SNPs.tsv ${StrX1}_inter1 ${StrX1}_temp1 Reference_genes_with_SNPs_removed_from_analysis Reference_Genes_with_function_calls ${StrX1}_genes_with_SNPs_removed_from_analysis ${StrX1}_Anvio_gene_callers_SNPs_no_functions.tsv awkfileQuerySnps ${StrX2}_intermediate1 ${StrX1}_vs_${StrX2}_temp1 ${StrX2}_intermediate1_loop ${StrX2}_GeneCalls.txt ${StrX2}_intermediate2_loop ${StrX2}_intermediate2_loop_head ${StrX1}_Anvio_gene_SNPs_NOT_in_Query_strain.tsv ${StrX2}_intermediate3_loop ${StrX2}_intermediate3_loop_seen ${StrX1}_${StrX2}_combined awkfileSnpReference ${StrX1}_${StrX2}_before_final Ref_aminoAcids ${StrX1}_intermediate_for_loop ${StrX1}_intermediate_for_loop_head 



# 2626 does is not part of a coding region for Strain11 while it is in Strain04 (the reference) 


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
