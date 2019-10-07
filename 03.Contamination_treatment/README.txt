This pipeline, which is optional, can be selected when we have more than one genomes in our Spades assembly, when they are inspected in Bandage.
In that case, and if we actually want to separate one of the genomes, in order to continue the downstream analysis only with the sequences from the selelected genome, we follow the steps indicated bellow:

	1. Bandage: Manually separate genome of interest and download scaffolds to a fasta file.
	2. BWA: Map the fastq reads to the fasta file downloaded from Bandage in the previous step, to actually retrieve the fastq reads that constitute our genome of interest.
	3. BBMerge: Merge paired end fastq reads filtered by BWA in the previous step, to get long single end reads.
	4. Kmergenie: Search the selected fastq reads for the best kmer, for the assembly
	5. Spades: Assemble the filtered fastq reads with Spades, using as kmers the kmers selected by Kmergenie in the previous step.
	6. (0PTIONAL) Barnap: Search for ribosomal RNAs in the assembled scaffolds, created by Spades in the previous step.
