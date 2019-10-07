#!/bin/bash


# for calculating the amount of time the job takes
begin=$(date +%s)

echo
echo "=================================================="
echo "job started " "$(date)"
echo "=================================================="
echo

# set -euo pipefail # Check https://dev.to/pencillr/newbie-at-bash-scripting-heres-some-advice-j1a


<< ////

        This is the Anvio pangenomics pipeline. Check the tutorial:
        http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions "Creating a quick pangenome with functions"
	
	This script run in the Anvio Pangenomics Prokka directory. To run the script, we need to have all the needed files (e.g. fix_functional_occurrence_table.py) and dependencies (blastp, mcl) and then we simply run the command:
	sbatch $0

	IMPORTANT
        We run this script after we have created the pangenome!

	IMPORTANT
	We have to prepare a text file with all the strain names e.g. Strain01, that we want to analyze here, one strain per line, with a second column which will categorize the strains into groups. The selection of the groups could be based on various criteria e.g. phylogenetic, ANI. This file is ${CategoricalTxt}, so please check the creation of this file at the start of this script!
	
	IMPORTANT!
	We have set the bins for the core pangenomes of interest in the interactive mode in order to be able to extract the gene clusters.

////


# INITIAL PARAMETERS
WorkDir=/mnt/283983cd-07b2-44de-b352-e9f8ccf9668e/Analysis/8.Anvio/2.Pangenome/1.Pangenome_Prokka_internal_genomes_ONLY
OutputPanDir=${WorkDir}/PSEUDOMONAS_PROKKA
PanDatabase=Pseudomonas_Prokka_Pangenome-PAN.db
GenDatabase=${WorkDir}/PSEUDOMONAS_COLUMBO_PROKKA-GENOMES.db
CategoricalTxt=Pseudomonas_internal_layer_data.txt
Category=Species  # This should match one of the column names in the $CategoricalTxt text file
# Next define the 2 groups that we are going to compare here.
Group1=stutzeri
Group2=aeruginosa

#AnnotationSource="Rast FigFam"  # This should be one of the annotation sources already available in the Pangenome Genome Storage. If you want to check the full list of Annotation Sources run this command: 
# anvi-get-enriched-functions-per-pan-group -g ${GenDatabase} -p ${OutputPanDir}/${PanDatabase} --list-annotation-sources -o dokimi.txt --category Species
# The ouput should be something like this:
# Available functional annotation sources .....................: Gene3D, Kegg Pathways + Enzymes, PfamInterPro, Prokka, TMHMM2.0, ProSitePatterns, KEGG Genes, Hamap, COG_FUNCTION, Rast FigFam, PANTHER GO-Slim Cellular, MobiDBLite, CDD, SMART, LipoP, Blast2Go, SignalP-4.1, PANTHER GO-Slim Function, KEGG Function, SFLD, PANTHER Protein, EggnogKegg, COG_CATEGORY, KEGG Brite/Pathway, PRINTS, PIRSF, EGGNOG, ProSiteProfiles, MetaCyc, Interpro, Phobius, TIGRFAM, ProDom, smCOG, SUPERFAMILY, Reactome, Interpo GoTerms, PANTHER GO-Slim Process


#-------------------------------------------------------------------------------------------------------------------------------------
if [[ ! -e "fix_functional_occurence_table.py" ]]; then
    echo "Please download the file fix_functional_occurrence_table.py (along with many other files!) from https://figshare.com/articles/31_Prochlorococcus_Isolate_Genomes/6318833 and make it available in the working directory!"
    echo "Please note that the file fix_functional_occurrence_table.py may not be used by this script. So if it is a hassle to download the file, just check first if you actually need it!" >&2 
    exit 1
fi

    
#-------------------------------------------------------------------------------------------------------------------------------------


echo
echo
echo "		Analysis begins		"
echo
echo

# First we have to arrange the gene clusters into bins, in the interactive mode (with the "search gene clusters using filters" and "append splits to selected bin" from the search tab in the interactive mode).
# The collection we created in the interactive mode was named as "default".
# We can check for the available collections with the following command:
# anvi-export-collection -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db --list-collections 
# We have created the following 4 bins:
# "aeruginosa_core_genome", "aeruginosa_all_except_Strain05", "Core_Genome" and "stutzeri_core_genome"
# Now, we can extract the binned clusters from the "default" collection with the following command:
anvi-export-collection -p PSEUDOMONAS_PROKKA/Pseudomonas_Prokka_Pangenome-PAN.db -C default -O GC_bins
# The above command outputs two files:
# 1. "GC_bins-info.txt" which contains all the available bins for the particular collection (and the corresponding color for the interactive mode).
# 2. "GC_bins.txt" which contains all the gene clusters grouped into each bin.

# Next, we want to create a file that will be used as input for the grep command, to search for all the Gene Clusters that we are interested in.
awk 'BEGIN{FS="\t";OFS="\t"} $2=="Core_Genome" {print $1}' GC_bins.txt | tr '\n' '\|' > Grep_Core_GCs_pattern.txt

awk 'BEGIN{FS="\t";OFS="\t"} $2=="stutzeri_core_genome" {print $1}' GC_bins.txt | tr '\n' '\|' > Grep_Core_${Group1}_GCs_pattern.txt

awk 'BEGIN{FS="\t";OFS="\t"} $2=="aeruginosa_core_genome" {print $1}' GC_bins.txt | tr '\n' '\|' > Grep_Core_${Group2}_GCs_pattern.txt

awk 'BEGIN{FS="\t";OFS="\t"} $2=="aeruginosa_all_except_Strain05" {print $1}' GC_bins.txt | tr '\n' '\|' > Grep_Core_${Group2}_no_Strain05_GCs_pattern.txt

# Next, create the appropriate directories to store the results.
mkdir -p CORE Core_${Group1}_UNIQUE Core_${Group1} Core_${Group2}_UNIQUE Core_${Group2} Core_${Group2}_UNIQUE_without_Strain05 Core_${Group2}_without_Strain05


# Next declare an array with all the annotation sources we have in the Anvio database.
declare -a sour=("Gene3D" "Kegg Pathways + Enzymes" 
		"PfamInterPro" "Prokka" 
		"TMHMM2.0" "ProSitePatterns"
		"KEGG Genes" "Hamap"
		"COG_FUNCTION" "Rast FigFam"
		"PANTHER GO-Slim Cellular" "MobiDBLite"
		"CDD" "SMART"
		"LipoP" "Blast2Go"
		"SignalP-4.1" "PANTHER GO-Slim Function"
		"KEGG Function" "SFLD"
		"PANTHER Protein" "EggnogKegg"
		"COG_CATEGORY" "KEGG Brite/Pathway"
		"PRINTS" "PIRSF"
		"EGGNOG" "ProSiteProfiles"
		"MetaCyc" "Interpro"
		"Phobius" "TIGRFAM"
		"ProDom" "smCOG"
		"SUPERFAMILY" "Reactome"
		"Interpo GoTerms" "PANTHER GO-Slim Process"
		)

# Start a for loop, that will iterate through all the annotation sources that we declared in the array.
# IMPORTANT
# The for loop exits at the end of the current script!
for source in "${sour[@]}"; do
AnnotationSource=$source
# Replace any blank spaces and other non regular characters with underscores for the Annotation Source (Shouldn't have introduced them in the first place! My bad!) 
AnnotName=${AnnotationSource// /_}
AnnotName=${AnnotName//\+/}
AnnotName=${AnnotName//\//_}
OutDir=${Group1}_vs_${Group2}_functional_${AnnotName}

#AnnotationSource=Prokka
#AnnotName=${AnnotationSource// /_}
#OutDir=${Group1}_vs_${Group2}_functional_${AnnotName}


[[ -d ${OutDir} ]] || mkdir ${OutDir}

# Move to the directory we just created!
cd ${OutDir}

# Copy, from the working directory, all the needed files, in order to work in the current directory.
cp ${WorkDir}/fix_functional_occurence_table.py .
cp ${WorkDir}/Grep_Core_GCs_pattern.txt .
cp ${WorkDir}/Grep_Core_${Group1}_GCs_pattern.txt .
cp ${WorkDir}/Grep_Core_${Group2}_GCs_pattern.txt .
cp ${WorkDir}/Grep_Core_${Group2}_no_Strain05_GCs_pattern.txt .

# First, we have to create the tab delimited file which will group our strains, according to our needs.
TAB=$'\t' 
cat > $CategoricalTxt <<EOF
isolate${TAB}${Category}
Strain01${TAB}${Group1}
Strain02${TAB}${Group1}
Strain03${TAB}${Group1}
Strain04${TAB}${Group1}
Strain05${TAB}${Group2}
Strain06${TAB}${Group2}
Strain07${TAB}${Group2}
Strain08${TAB}${Group1}
Strain09${TAB}${Group1}
Strain10${TAB}${Group1}
Strain11${TAB}${Group1}
Strain12${TAB}${Group1}
Strain13${TAB}${Group1}
Strain14${TAB}${Group1}
Strain16${TAB}${Group2}
Strain18${TAB}${Group1}
Strain19${TAB}${Group2}
Strain20${TAB}${Group2}
Strain21${TAB}${Group2}
Strain22${TAB}${Group1}
Strain23${TAB}${Group2}
Strain24${TAB}${Group1}
Strain25${TAB}${Group1}
EOF



# Import manually created layers to the pangenome database.
anvi-import-misc-data ${CategoricalTxt} \
					-p ${OutputPanDir}/${PanDatabase} \
					--target-data-table layers \
					--just-do-it


# Create two tab delimited tables:
# The first has the information about which gene clusters form the core of the pangenome.
# The second one shows the occurrence (presence/absence) of each ${AnnotationSource} annotation for each ${Category}
anvi-get-enriched-functions-per-pan-group -p ${OutputPanDir}/${PanDatabase} \
					-g ${GenDatabase} \
					--category $Category \
					--annotation-source "${AnnotationSource}" \
					-o PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt \
					--functional-occurrence-table-output PSEUDOMONAS_${AnnotName}_functions_occurrence.txt

# Next, we want to only keep alphanumeric characters, and replace any sequence of non-alphanumeric characters by a single "_" from the annotations of the functions occurence file.
# We have to insert a TAB variable because sed is not recognizing \t as a tab. Check also here https://stackoverflow.com/questions/2610115/why-is-sed-not-recognizing-t-as-a-tab
sed "s/[^[:alnum:]${TAB}_]/_/g" PSEUDOMONAS_${AnnotName}_functions_occurrence.txt | \
	tr -s \_ _ | \
	sed "s/^${TAB}/name${TAB}/" \
	> PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed_a_little.txt


# Print any duplicate functions to a file 
cut -f 1 PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed_a_little.txt | sort | uniq -d > duplicate_annotations.txt

# Check if the file is empty. If it is not empty then run a python program to remove all the duplicate entries. If the file is empty then no duplicates exist, and we simply rename the file.
if [[ -s duplicate_annotations.txt ]]; then 
./fix_functional_occurence_table.py PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed_a_little.txt PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt
else 
    mv PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed_a_little.txt PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt
fi

# rm duplicate_annotations.txt

# Next, we turn the functions occurence file created earlier into a newick tree file having all the ${AnnotName} at the tips of the tree.
anvi-matrix-to-newick PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt \
		      -o PSEUDOMONAS_${AnnotName}_functions_tree.nwk

# We repeat the same once more but this time we take into account the columns i.e. the ${Category}, not the rows (by taking the traspose matrix).
anvi-matrix-to-newick PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt \
		      -o PSEUDOMONAS_${AnnotName}_functions_layers_tree.nwk \
		      --transpose


# Next, we create a manual mode profile database.
anvi-interactive -p PSEUDOMONAS_${AnnotName}_functions_manual_profile.db \
		 	--tree PSEUDOMONAS_${AnnotName}_functions_tree.nwk \
		 	--view-data PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt \
		 	--manual \
		 	--dry-run

# Next, we create a new file with three columns and column names. One of the columns is the functions layer tree we created in a previous step.
echo -e "item_name\tdata_type\tdata_value" > PSEUDOMONAS_${AnnotName}_functions_layers_order.txt
echo -e "PSEUDOMONAS_functions_tree\tnewick\t`cat PSEUDOMONAS_${AnnotName}_functions_layers_tree.nwk`" \
			     >> PSEUDOMONAS_${AnnotName}_functions_layers_order.txt

# Next, we import to the functions database as a layer order, the file we prepared in the previous step.
anvi-import-misc-data PSEUDOMONAS_${AnnotName}_functions_layers_order.txt \
		      -p PSEUDOMONAS_${AnnotName}_functions_manual_profile.db \
		      -t layer_orders \
		      --just-do-it

# Next, we EXPORT from the (initial) pangenome database the layers.
anvi-export-misc-data -p ${OutputPanDir}/${PanDatabase} \
		      -t layers \
		      -o PSEUDOMONAS-layer-additional-data.txt

# Finally, we IMPORT the layers data, we exported in the previous step, into the functions database.
anvi-import-misc-data PSEUDOMONAS-layer-additional-data.txt \
		      -p PSEUDOMONAS_${AnnotName}_functions_manual_profile.db \
		      -t layers


#--------------------------------------------------------------------------------------------------------------------------
# Next, we want to export some useful information about the functions in the pangenome.
# The problem we face is that the output file "PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt" although it has a constant number of columns, the ordering of these columns differs each time we run the program (for a reason I do not understand)

# Before starting the analysis, we want to replace ASCII encodings with characters, because in certain annotations like Prokka and Rast Figfam, we have such a problem.
sed -i 's/%2C/,/g;s/%26/\&/g;s/%3B/;/g' PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt


# First, we want to get the number of columns, that the main output file has.
numCol=`awk 'BEGIN{FS="\t";OFS="\t"};{print NF}' PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | head -n 1`

# Next, we declare an array with all the expected headers of the output file "PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt". 
# The ordering (index) of the array is very important because we will use this ordering to export later on with the awk command the correct columns. 
declare -a arr=("category" "core"
		"${AnnotName}" "core_in_group"
		"function_accession" "gene_clusters_ids"
		"occurrence_in_group" "occurrence_outside_of_group"
		)
# Headers not included in the above array:
# "portion_occurrence_in_group" "portion_occurrence_outside_of_group" "p_value" "corrected_p_value" "enrichment_score"

# Next, we will try

counter=1 # This is a counter that will count the number of columns in the output file "PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt"
for i in $(seq 1 $numCol); do
    # Get the headers
    ColHeader=$(head -n 1 PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -f"$i" -d$'\t')
    # Remove any gaps or non standard characters in headers. (These should be exactly the same as the substitutions we did for $AnnotName in the beggining of this script)!
    ColHeader=${ColHeader// /_}
    ColHeader=${ColHeader//\+/}
    ColHeader=${ColHeader//\//_}
    for z in "${arr[@]}"; do
	if [[ "$ColHeader" == "$z" ]]; then
	    # declare -x Column${counter}="$ColHeader"
	    awk -v coun="$counter" 'BEGIN{FS="\t";OFS="\t"} {print $coun}' PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt > "${ColHeader}"
	    break
	fi
    done
    counter=$((counter + 1))
done

#------------------------------------------------------------------------------------------------------------------------

# 
LineNum=`wc -l PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -f1 -d" "`
for i in $( seq 2 $LineNum ); do 
    echo "0 $i" 
done > All_lines

for i in $( seq 2 $LineNum ); do 
    echo "- $i" 
done > All_lines_GCs

# We first try to extract the number of core gene clusters that exist in each line.
grep -E -o -n -f Grep_Core_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -d : -f 1 | uniq -c | sed -E 's/([[:digit:]]) ([[:digit:]])/\1\t\2/g;s/ //g;s/\t/ /g' > Core_GCs.lines
cat Core_GCs.lines All_lines | sort -n -u -k2 | awk 'BEGIN{print "Num_Core_GCs"} {print $1}' > Num_Core_GCs
# Next, we try to extract the core gene clusters that exist in each line.
grep -E -n -o -f Grep_Core_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | awk -F':' 'NF>1{a[$1] = a[$1]","$2};END{for(i in a)print a[i],i}' | sed s/^,//g | sort -n -k2 > Core_GCs.list
cat Core_GCs.list All_lines_GCs | sort -n -u -k2 | awk 'BEGIN{print "Core_GCs"} {print $1}' > Core_GCs


# Next, we try to extract the number of ${Group1} gene clusters that exist in each line.
grep -E -o -n -f Grep_Core_${Group1}_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -d : -f 1 | uniq -c | sed -E 's/([[:digit:]]) ([[:digit:]])/\1\t\2/g;s/ //g;s/\t/ /g' > Core_${Group1}_GCs.lines
cat Core_${Group1}_GCs.lines All_lines | sort -n -u -k2 | awk -v group="$Group1" 'BEGIN{print "Num_"group"_GCs"} {print $1}' > Num_Core_${Group1}_GCs
# Next, we try to extract the ${Group1} gene clusters that exist in each line.
grep -E -n -o -f Grep_Core_${Group1}_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | awk -F':' 'NF>1{a[$1] = a[$1]","$2};END{for(i in a)print a[i],i}' | sed s/^,//g | sort -n -k2 > ${Group1}_GCs.list
cat ${Group1}_GCs.list All_lines_GCs | sort -n -u -k2 | awk -v group="$Group1" 'BEGIN{print group"_GCs"} {print $1}' > Core_${Group1}_GCs


# Next, we try to extract the number of ${Group2} gene clusters that exist in each line.
grep -E -o -n -f Grep_Core_${Group2}_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -d : -f 1 | uniq -c | sed -E 's/([[:digit:]]) ([[:digit:]])/\1\t\2/g;s/ //g;s/\t/ /g' > Core_${Group2}_GCs.lines
cat Core_${Group2}_GCs.lines All_lines | sort -n -u -k2 | awk -v group="$Group2" 'BEGIN{print "Num_"group"_GCs"} {print $1}' > Num_Core_${Group2}_GCs
# Next, we try to extract the ${Group2} gene clusters that exist in each line.
grep -E -n -o -f Grep_Core_${Group2}_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | awk -F':' 'NF>1{a[$1] = a[$1]","$2};END{for(i in a)print a[i],i}' | sed s/^,//g | sort -n -k2 > ${Group2}_GCs.list
cat ${Group2}_GCs.list All_lines_GCs | sort -n -u -k2 | awk -v group="$Group2" 'BEGIN{print group"_GCs"} {print $1}' > Core_${Group2}_GCs


# Next, we try to extract the number of no_Strain05_${Group2} gene clusters that exist in each line.
grep -E -o -n -f Grep_Core_${Group2}_no_Strain05_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | cut -d : -f 1 | uniq -c | sed -E 's/([[:digit:]]) ([[:digit:]])/\1\t\2/g;s/ //g;s/\t/ /g' > Core_${Group2}_no_Strain05_GCs.lines
cat Core_${Group2}_no_Strain05_GCs.lines All_lines | sort -n -u -k2 | awk -v group="$Group2" 'BEGIN{print "Num_"group"_no_Strain05_GCs"} {print $1}' > Num_Core_${Group2}_no_Strain05_GCs

# Next, we try to extract the no_Strain05_${Group2} gene clusters that exist in each line.
grep -E -n -o -f Grep_Core_${Group2}_no_Strain05_GCs_pattern.txt PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt | awk -F':' 'NF>1{a[$1] = a[$1]","$2};END{for(i in a)print a[i],i}' | sed s/^,//g | sort -n -k2 > ${Group2}_no_Strain05_GCs.list
cat ${Group2}_no_Strain05_GCs.list All_lines_GCs | sort -n -u -k2 | awk -v group="$Group2" 'BEGIN{print group"_no_Strain05_GCs"} {print $1}' > Core_${Group2}_no_Strain05_GCs


# Next, we will paste together all the columns that were extracted in separate files
paste category core_in_group core occurrence_in_group occurrence_outside_of_group gene_clusters_ids function_accession ${AnnotName} Num_Core_GCs Num_Core_${Group1}_GCs Num_Core_${Group2}_GCs Num_Core_${Group2}_no_Strain05_GCs Core_GCs Core_${Group1}_GCs Core_${Group2}_GCs Core_${Group2}_no_Strain05_GCs > ${Group1}_${Group2}_${AnnotName}_final.tsv


# Remove the files we do not need any more.
rm category core_in_group core occurrence_in_group occurrence_outside_of_group gene_clusters_ids function_accession ${AnnotName} Num_Core_GCs Num_Core_${Group1}_GCs Num_Core_${Group2}_GCs Num_Core_${Group2}_no_Strain05_GCs Core_GCs Core_${Group1}_GCs Core_${Group2}_GCs Core_${Group2}_no_Strain05_GCs All_lines Core_GCs.lines Core_${Group1}_GCs.lines Core_${Group2}_GCs.lines Core_${Group2}_no_Strain05_GCs.lines Core_GCs.list ${Group1}_GCs.list ${Group2}_GCs.list ${Group2}_no_Strain05_GCs.list


# To get the core genome (genes found in all strains of the genus irrespective of the specific species). We also have to make the output unique because we have two groups and each row is duplicated. 
awk -v annot="${AnnotName}" 'BEGIN{FS="\t";OFS="\t";print "source", "accession","function","Core_Gene_Clusters","Num_Core_GCs"} $3=="True" {print annot,$7,$8,$13,$9}' ${Group1}_${Group2}_${AnnotName}_final.tsv | awk -F"\t" '!_[$2]++' > Functions_core_${AnnotName}.tsv

# Next, copy the file to a directory, where all results will be kept.
cp Functions_core_${AnnotName}.tsv ../CORE

# First we want to get the functions that exist in ALL strains that belong to the first group and only to the first group.
awk -v annot="${AnnotName}" -v group="$Group1" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function","Core_Gene_Clusters","Num_"group"_GCs"}  group==$1 && $2=="True" && $5=="0.00" {print $1,annot,$7,$8,$14,$10}' ${Group1}_${Group2}_${AnnotName}_final.tsv > Functions_${Group1}_UNIQUE_${AnnotName}.tsv

cp Functions_${Group1}_UNIQUE_${AnnotName}.tsv ../Core_${Group1}_UNIQUE

# Next, we get the functions that exist in all strains of the first group.
awk -v annot="${AnnotName}" -v group="$Group1" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function","Core_Gene_Clusters","Num_"group"_GCs"}  group==$1 && $2=="True" {print $1,annot,$7,$8,$14,$10}' ${Group1}_${Group2}_${AnnotName}_final.tsv > Functions_${Group1}_CORE_${AnnotName}.tsv

cp Functions_${Group1}_CORE_${AnnotName}.tsv ../Core_${Group1}


# Next, we get the functions that exist in all strains of the second group and only to the second group.
awk -v annot="${AnnotName}" -v group="$Group2" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function",group"_Gene_Clusters","Num_"group"_GCs"} group==$1 && $2=="True" && $5=="0.00" {print $1,annot,$7,$8,$15,$11}' ${Group1}_${Group2}_${AnnotName}_final.tsv  > Functions_${Group2}_UNIQUE_${AnnotName}.tsv

cp Functions_${Group2}_UNIQUE_${AnnotName}.tsv ../Core_${Group2}_UNIQUE


# Next, we get the functions that exist in all strains of the second group.
awk -v annot="${AnnotName}" -v group="$Group2" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function",group"_Gene_Clusters","Num_"group"_GCs"} group==$1 && $2=="True" {print $1,annot,$7,$8,$15,$11}' ${Group1}_${Group2}_${AnnotName}_final.tsv  > Functions${Group2}_CORE_${AnnotName}.tsv

cp Functions${Group2}_CORE_${AnnotName}.tsv ../Core_${Group2}


# Finally, a special case for this analysis and only for this analysis. We want to get all the ${Group2} strains which do not include one strain i.e. Strain05, and are unique to ${Group2}.
awk -v annot="${AnnotName}" -v group="$Group2" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function",group"_no_Strain05_GCs","Num_"group"_no_Strain05_GCs"} group==$1 && $4=="7.00" && $5=="0.00" {print $1,annot,$7,$8,$16,$12}' ${Group1}_${Group2}_${AnnotName}_final.tsv  > Functions_${Group2}_UNIQUE_without_Strain05_${AnnotName}.tsv

cp Functions_${Group2}_UNIQUE_without_Strain05_${AnnotName}.tsv ../Core_${Group2}_UNIQUE_without_Strain05

# Finally, a special case for this analysis and only for this analysis. We want to get all the aeruginosa strains which do not include one strain i.e. Strain05, because the latter one is missing genes.
awk -v annot="${AnnotName}" -v group="$Group2" 'BEGIN{FS="\t";OFS="\t";print "group","source", "accession","function",group"_no_Strain05_GCs","Num_"group"_no_Strain05_GCs"} group==$1 && $4=="7.00" {print $1,annot,$7,$8,$16,$12}' ${Group1}_${Group2}_${AnnotName}_final.tsv  > Functions_${Group2}_CORE_without_Strain05_${AnnotName}.tsv

cp Functions_${Group2}_CORE_without_Strain05_${AnnotName}.tsv ../Core_${Group2}_without_Strain05


cd ${WorkDir}

done

rm GC_bins* 

#<< ////

#Next Steps:
#	We can visualise the functions database with the following command:
#	anvi-interactive -p PSEUDOMONAS_${AnnotName}_functions_manual_profile.db \
#                 -t PSEUDOMONAS_${AnnotName}_functions_tree.nwk \
#                 --view-data PSEUDOMONAS_${AnnotName}_functions_occurrence-fixed.txt \
#                 --title "Pseudomonas Pan - ${AnnotName} functional occurrence" \
#                 --manual

#	
#	# Optionally, we can export the core genome gene clusters in a new file.
#	# The next example ouputs the core genome for the aeruginosa strains. "aeruginosa" is one of the categories of the ${Category} column in the file ${CategoricalTxt}
#	grep aeruginosa PSEUDOMONAS_${AnnotName}_enriched-function_${Category}.txt |\
#	  	awk  -F $'\t' '$6 == "True" { print $12 }' |\
#	    	tr ',' '\n' |\
#	      	sed 's/ //g' > aeruginosa_core_${AnnotName}_functions_gene_clusters.txt


#	# if we have a bin collection and we want to check what percentage of it belongs to the core genome, then we can do the following:
#	for gc in `grep Bin_1 Pseudomonas_Prokka_Pangenome_gene_clusters_summary.txt`
#	do
#    	grep $gc aeruginosa_core_${AnnotName}_functions_gene_clusters.txt
#	done > Bin_1_included_in_functional_core.txt

#////



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
