// =====================================================================
// amas_mapper: an All-Mapping tool using Adaptive Seeds (AMAS)
// =====================================================================






#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>

#include "ktab.h"
#include "parser_options.h"
#include "sequence_index.h"
#include "writer.h"
#include "mapper.h"






#define ReadBlockSize 1000000 // process reads in 1-M blocks






using namespace seqan;






template <typename TFormat>
int runMapper(Options & options)
{
	// to record the total running time
	double totalCPU = cpuTime();
	double totalSYS = sysTime();











	typedef Genome<void> TGenome;

    typedef ReadsConfig<> TReadsConfig;

	typedef Reads<void, TReadsConfig> TReads;
	
    typedef ReadsLoader<void, TReadsConfig> TReadsLoader;

	typedef Writer<TGenome, TReads, TFormat, EditDistance> TWriter;

	typedef ReadMapperConfig<> TMapperConfig;

	typedef Mapper<TGenome, TReads, TWriter, TMapperConfig> TMapper;











// =====================================================================
// Reference genome and its container
// (see "sequence_index.h" for their data structures)
// =====================================================================

    TFragmentStore store;
    
	TGenome genome(store);





    
    
    
    
 
    
// =====================================================================
// Load genome
// =====================================================================

// /*

	struct rusage ru;
	double timeCPU = cpuTime();
	double timeSYS = sysTime();

	// load contigs to object "store"  & record their lengths
	std::cout << "Loading genome..." << std::endl;
    if (!loadGenome(genome, options.genomeFile)) {
		std::cerr << "Failed to load the genome!" << std::endl;
		std::cerr << "Please make sure that a proper FASTA file is provided!" << std::endl;
		return 1;
	}
	else {
		std::cout << "Genome loaded!" <<std::endl;
		std::cout << "Total number of contigs: " << genome.numContigs << std::endl;
		std::cout << "Total length of the genome: " << genome.totalLength << " nucleotides" << std::endl;
	}
	
    timeCPU = cpuTime() - timeCPU;
    timeSYS = sysTime() - timeSYS;
    std::cout << "Loading time (CPU): " << timeCPU << " seconds" << std::endl;
    std::cout << "Loading time (wall-clock): " << timeSYS << " seconds" << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

    std::cout << std::endl;

// */
            





	
	
	


// =====================================================================
// Load the adaptive-seed tree
// (see "ktab.h" for the data structures of the tree & its nodes)
// =====================================================================

// /*
	
	timeCPU = cpuTime();
	timeSYS = sysTime();
	
	// frequency threshold of adaptive seeds
	unsigned F;
	
	// the array of IDs to core (10-mer) nodes of the adaptive-seed tree
	std::vector<unsigned> coreKTab(KOUNT,0);
	
	// the table containing all nodes of the adaptive-seed tree
	std::vector<KTab1> tableKTab;
	tableKTab.push_back(KTab1()); // add an empty KTab at index 0; all null pointers will point to this entry

	// The default folder & name to save/load the tree: "tree/genomeFile.tre"
    CharString treeFile = "tree/";
    append(treeFile, options.genomeFile);
    append(treeFile, ".tre");
	std::cout << "Loading the adaptive-seed tree from the default folder & name: " << treeFile << std::endl;

    std::fstream tree_input(toCString(treeFile), std::ios::binary | std::ios::in);
    // load the adaptive-seed tree ("F", "coreKTab", "tableKTab") from the binary file
	loadTable(F, coreKTab, tableKTab, tree_input);

	std::cout << "Tree loaded!" << std::endl;

	timeCPU = cpuTime() - timeCPU;
	timeSYS = sysTime() - timeSYS;
	std::cout << "Loading time (CPU): " << timeCPU << " seconds" << std::endl;
	std::cout << "Loading time (wall-clock): " << timeSYS << " seconds" << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

    std::cout << std::endl;

// */	










	
// ==============================================================================
// Configure the "reads"
// (see "sequence_index.h" for the data structure & the configuration of "Reads")
// ==============================================================================

    TReads reads(store); // the read sequences are actually stored in the container "store"
    
    
    
    
    






// =====================================================================
// Prepare the component "readsLoader"
// (open reads' file, compute file size, estimate reads' statistics...)
// (see "sequence_index.h")
// and then reserve the memory for the "reads"
// =====================================================================

    TReadsLoader readsLoader(reads);
    
    if (!open(readsLoader, options.readsFile))
    {
        std::cerr << "Error while opening reads file" << std::endl;
        return 1;
    }

	if (reads._countEstimate < ReadBlockSize)
    {
		reserve(reads, reads._countEstimate);
	}
    else
    {
        reserve(reads, ReadBlockSize);
	}





	
	
	
	
	
	
// =====================================================================
// Configure the component "writer"
// (see "writer.h")
// =====================================================================

	TWriter writer(store, genome, reads, options.noDump);

	
	
	
	
	





// =====================================================================
// Configure the component "mapper"
// (see "mapper.h")
// =====================================================================

	TMapper mapper(store, genome, reads, writer, options.noVerify);
	










// =====================================================================
// Prepare the output file
// =====================================================================

    CharString outputFile = options.readsFile;
    if (options.raw)
	{
	    append(outputFile, ".amas.raw");
	}
	else
	{
	    append(outputFile, ".amas.sam");
	}

    if (!open(writer, outputFile))
    {
        std::cerr << "Error while opening output file" << std::endl;
        return 1;
    }

    // always output CIGAR string.
    //writeAlignments(writer, true); //nhtran: changed! This function nows performs the writing task.

	
	
	
	
	





// =====================================================================
// LOAD & MAP reads in 1-million blocks ("ReadBlockSize" = 1 million)
// =====================================================================

	unsigned readsCount = 0;
    unsigned hitsCount = 0;
    unsigned matchesCount = 0;
    unsigned readsMatched = 0;

    std::cout << "LOADING & MAPPING..." << std::endl;

	std::cout << "Number of processors online on your device" << "\t" << omp_get_num_procs() << std::endl;;
	omp_set_num_threads(options.threads);
	#pragma omp parallel
	{
		#pragma omp master
		{
		    std::cout << "Number of threads allocated" << "\t" << omp_get_num_threads() << std::endl;
		}
	}
	
	timeCPU = cpuTime();
	timeSYS = sysTime();

    while (!atEnd(readsLoader))
    {
		// LOAD LOAD LOAD LOAD LOAD LOAD LOAD LOAD LOAD LOAD LOAD LOAD
		// load reads' sequences, names, and generate their reverse complements
        if (!load(readsLoader, ReadBlockSize))
        {
            std::cerr << "Error while loading reads" << std::endl;
            return 1;
        }
		readsCount += reads.readsCount;
		
		// clear old and reserve new space for writing the alignments
		//setReads(writer); // nhtran: replaced by functions "load" & "clear" of object "reads"

		// MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP MAP
        mapReads(mapper, F, coreKTab, tableKTab, options);
		
	    // record hits, alignments, & reads mapped
		for (TReadSeqStoreSize read_ID = 0; read_ID < reads.readsCount; ++read_ID)
	    {
			hitsCount += reads.hitsCount.at(read_ID) + reads.hitsCount.at(read_ID + reads.readsCount);
			matchesCount += reads.matchesCount.at(read_ID) + reads.matchesCount.at(read_ID + reads.readsCount);
			if (reads.matchesCount.at(read_ID) + reads.matchesCount.at(read_ID + reads.readsCount) > 0)
				readsMatched ++;
		}
		
		// WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE 
		//writeAlignments(writer);
		
        // clear mapped reads, free space for the new block
        clear(reads);
	}
	
	std::cout << "Reads mapped!" << std::endl;

	timeCPU = cpuTime() - timeCPU;
	timeSYS = sysTime() - timeSYS;
	std::cout << "Mapping time (CPU): " << timeCPU << std::endl;
	std::cout << "Mapping time (wall-clock): " << timeSYS << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

	std::cout << std::endl;











// =====================================================================
// MAPPING ends; output the summary information
// =====================================================================

	// close the input reads' file and the output SAM file
    close(writer);
    close(readsLoader);
    
    
    
    
    

	// output the summary information
	std::cout << "MAPPING SUMMARY" << std::endl;
	//std::cout << "Total seed count: " << mapper._seeder.seedsCount << std::endl;
	//std::cout << "Sum of seed lengths: " << mapper._seeder.seedsLength << std::endl;
	std::cout << "Number of reads: " << readsCount << std::endl;
	std::cout << "Number of reads mapped: " << readsMatched << std::endl;
	std::cout << "Number of hits for extension: " << hitsCount << std::endl;
	std::cout << "Number of alignments found: " << matchesCount << std::endl;
	
	totalCPU = cpuTime() - totalCPU;
	totalSYS = sysTime() - totalSYS;
    std::cout << "Total running time (CPU): " << totalCPU << std::endl;
    std::cout << "Total running time (wall-clock): " << totalSYS << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

	std::cout << std::endl;
	
	return 0;
}   







int configurationOutputFormat(Options & options)
{
	if (options.raw)
	{
		return runMapper<Raw>(options);
	}
	else
	{
		return runMapper<Sam>(options);
	}
}






int main(int argc, char const **argv)
{
// =====================================================================
// Parse arguments from the command line to "options"
// (see "parser_options.h" for class "Options")
// =====================================================================

	Options options;
	
	ArgumentParser::ParseResult parseResult = parseMapper(options, argc, argv);

    if (parseResult != ArgumentParser::PARSE_OK)
        return parseResult == seqan::ArgumentParser::PARSE_ERROR;

	//std::cout << options.genomeFile << std::endl;
	//std::cout << options.readsFile << std::endl;
	//std::cout << options.errors << std::endl;
	//std::cout << options.fse << std::endl;
	//std::cout << options.msl << std::endl;
	
	return configurationOutputFormat(options);
}
