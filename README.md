# AMAS
Project: All-Mapping tool using Adaptive Seeds for high-throughput sequencing data

Author: Hieu Tran, Research Fellow at University of Waterloo, Canada. Email: hieutran1985@gmail.com

Publication: AMAS optimizing the partition and filtration of adaptive seeds to speed up read mapping. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2015.

#

AMAS is designed to exhaustively search for all possible mapping locations of high-throughput sequencing (HTS) reads in a reference sequence, with up to a given edit distance (including substitutions and indels).

AMAS features major improvements of the mapping, partition and filtration of adaptive seeds to speed up the all-mapping process while preserving sensitivity. In addition to all-mapping, AMAS also supports best-mapping and -k mapping modes.

Multi-threading is supported for the mapping step. We will implement multi-threading for the indexing step very soon. The memory footprint is a bit high. For the human genome, ~19GB of memory is required. This shall also be improved in future development.

AMAS is implemented in C++ based on the source code of and the Masai and SeqAn library. Our contributions include the index, the partition and filtration of adaptive seeds. For seed extension, we borrowed Masai's implementation of the Myers' bit-vector dynamic programming algorithm. We also borrowed I/O components from Masai for handling the input reference sequence, HTS reads, and the SAM output.

# How to use Binaries & Source Code:

(I) Pre-compiled Binaries

The provided binaries were compiled and tested on a 64-bit Linux (Red Hat) platform with GCC 4.1.2.
	
	1) Download and unzip the file "amas_v02.zip".
	
	2) Use your own data or download the sample files "ce10.fasta.gzip" and "ce10_10k.fastq.gzip".
	
	3) Unzip and copy the files "ce10.fasta" and "ce10_10k.fastq" to the folder "amas".
	
	4) Inside the folder "amas":
	
		4.1) Type "bin/amas_indexer ce10.fasta" to build the tree index from the genome.
		This indexing step should take about 5 minutes for the sample genome.
		
		4.3) Type "bin/amas_mapper ce10.fasta ce10_10k.fastq" to perform the read mapping.
		This mapping step should take about 10 seconds for the sample files.
		
	5) For large genomes, e.g. human genome, if you encounter errors when building the index, 
		you may need to change the system's temporary folder by modifying the environment variable TMPDIR.
		This is because the suffix array construction algorithm of SeqAn 
		requires external space to store temporary data (default folder is /tmp). 
		

	6) See below for more details and options.
	
(II) Compile from Source Code

To compile from the source code:

	1) Inside folder "amas": type "make".
	
	2) The compiled binaries will be stored in folder "release".
	
=========================
More Details and Options
=========================

AMAS consists of two tools: "amas_indexer" to index the reference sequence 
and "amas_mapper" to perform the read mapping.

Use option "-h/--help" to get more info and options of the tools.

For "amas_indexer", only the genome file needs to be provided.

(I) Options for "amas_mapper":
	
	-e, --errors
			Maximum number of errors allowed (mismatches & gaps). 
			In range [0..10]. Default: 5.
			
	-fse, --full-sen-error NUM
			Require full sensitivity for matches with up to "fse" errors. 
			In range [0..5]. Default: 2.
			
	-best, --best
			Best-mapping mode: find the best alignment for each read.
			
	-k, --k-mapping
			k-mapping mode: find up to k alignments for each reads. 
			In range [1..1000000]. Default: 1000000.
			
	-raw, --raw
			Output the alignments in the raw format for pair-end matching with "masai_output_pe".
			
	-p, --threads NUM
			Number of threads. In range [1..16]. Default: 1.

(II) Pair-end Reads Mapping

The current version of AMAS does not directly perform pair-end reads mapping. 
It is done with the help of "masai_output_pe" from MASAI.

For pair-end reads, each read of a pair is first mapped to the reference sequence with the output in the RAW format.

	Example: "bin/amas_mapper hg19.fasta hg19_100k_1.fastq -raw"
			 "bin/amas_mapper hg19.fasta hg19_100k_2.fastq -raw"

Then the alignments of two reads of the same pair are paired off by using the tool "masai_output_pe" from MASAI.

	Example: "bin/masai_output_pe hg19.fasta hg19_100k_1.fastq hg19_100k_2.fastq 
			hg19_100k_1.fastq.amas.raw hg19_100k_2.fastq.amas.raw 
			-o hg19_100k.fastq.amas.sam"

Use "-h/--help" to get more useful options of "masai_output_pe", such as the fragment length, etc.
