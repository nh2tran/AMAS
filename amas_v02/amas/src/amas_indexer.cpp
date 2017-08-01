// ========================================================================================
// amas_indexer: Build the adaptive-seed tree from the reference genome
// ========================================================================================






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

#include "parser_options.h"
#include "sequence_index.h"
#include "ktab.h"






using namespace seqan;






int main(int argc, char const **argv)
{
	// to record the total running time
	double totalCPU = cpuTime();
	double totalSYS = sysTime();

// =====================================================================
// Parse arguments from the command line to "options"
// (see "parser_options.h" for class "Options")
// =====================================================================

	Options options;
	
	ArgumentParser::ParseResult parseResult = parseHasher(options, argc, argv);

    if (parseResult != ArgumentParser::PARSE_OK)
        return parseResult == seqan::ArgumentParser::PARSE_ERROR;

	//std::cout << options.genomeFile << std::endl;
	//std::cout << options.F << std::endl;
		
	
	
	
	
	





// =====================================================================
// Reference genome and its container
// (see "sequence_index.h" for their data structures)
// =====================================================================

    TFragmentStore store;
    
	typedef Genome<void> TGenome;
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
// Reverse the contigs & build the FM index
// (see "sequence_index.h" for the data structure of the FM index)
// =====================================================================

// /*

	std::cout << "Building FM index..." << std::endl;

	timeCPU = cpuTime();
	timeSYS = sysTime();

    // FM index works on the reverse of the reference genome
    reverse(genome.contigs);
        
    TGenomeFM fmIndex(genome.contigs); // the index is not created yet
	TIter iter(fmIndex); // the index is actually created after this line
    
	std::cout << "FM index built!" << std::endl;

    timeCPU = cpuTime() - timeCPU;
    timeSYS = sysTime() - timeSYS;
    std::cout << "Building time (CPU): " << timeCPU << " seconds" << std::endl;
    std::cout << "Building time (wall-clock): " << timeSYS << " seconds" << std::endl;

   	getrusage(RUSAGE_SELF, &ru);        
	std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
	//std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
	//std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

	std::cout << std::endl;

// */
	

    
    
    
    
    
    
    
    
    
// =====================================================================
// DISABLED !!!
// Load the FM index
// (see "sequence_index.h" for the data structure of the FM index)
// =====================================================================

/*	

	timeCPU = cpuTime();
	timeSYS = sysTime();
	
	// The default folder & name to save/load the FM index: "index/genomeFile.index"
    CharString indexFile = "index/";
    append(indexFile, options.genomeFile);
    append(indexFile, ".index");
	std::cout << "Loading FM index from the default folder & name: " << indexFile << std::endl;
    TGenomeFM fmIndex(genome.contigs);
    if (!open(fmIndex, toCString(indexFile))) {
		std::cout << "Failed to open the FM index" << std::endl;
		return 1; 
	}
	else {
		std::cout << "FM index loaded!" <<std::endl;
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

*/	




     






// =====================================================================
// Build the adaptive-seed tree to store the locations of adaptive seeds
// (see "ktab.h" for the data structures of the tree & its nodes)
// =====================================================================

// /*	

	timeCPU = cpuTime();
	timeSYS = sysTime();
	
	// frequency threshold of adaptive seeds: user input or auto-selected
	//unsigned F = options.F
	unsigned F = genome.totalLength/KOUNT/10; // "KOUNT"=1048576 (number of 10-mers, defined in "ktab.h")
	if (F < 1)
	{
		F = 1;
	}
	std::cout << "Frequency threshold F of adaptive seeds: " << F << std::endl;
	
	std::cout << "Building the adaptive-seed tree ..." << std::endl;

	// the array of IDs to the core (10-mer) nodes of the adaptive-seed tree
	std::vector<unsigned> coreKTab(KOUNT,0);
	
	// the table containing all nodes of the adaptive-seed tree
	std::vector<KTab1> tableKTab;
	tableKTab.push_back(KTab1()); // add an empty KTab at index 0; all null pointers will point to this entry

	// would it be faster & less mem-allocating if we know the number of nodes ("numKTab") beforehand???
	tableKTab.reserve(50000000);
	
	unsigned numKTab = 0;
	unsigned maxDepth = 0;
	unsigned numLeaf = 0;
	
	// build the adaptive-seed tree (i.e. create its nodes, LINK them, & store them in "tableKTab"
	fillTable(fmIndex, F, coreKTab, tableKTab, numKTab, maxDepth, numLeaf);
	
	std::cout << "Tree built!" << std::endl;
	std::cout << "Number of nodes: " << numKTab << std::endl;
	std::cout << "Max branch length: " << maxDepth << std::endl;
	std::cout << "Number of leaf nodes: " << numLeaf << std::endl;
	
	timeCPU = cpuTime() - timeCPU;
	timeSYS = sysTime() - timeSYS;
	std::cout << "Building time (CPU): " << timeCPU << " seconds" << std::endl;
	std::cout << "Building time (wall-clock): " << timeSYS << " seconds" << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;
   
    std::cout << std::endl;

// */
	
	
	
	
	
	
	
	
	
	
	
// ==========================================================================
// Save the adaptive-seed tree to a local folder in binary format
// ==========================================================================

// /*	

	timeCPU = cpuTime();
	timeSYS = sysTime();
	
	// The default folder & name to save/load the tree: "tree/genomeFile.tre"
    CharString treeFile = "tree/";
    append(treeFile, options.genomeFile);
    append(treeFile, ".tre");
	std::cout << "Saving the adaptive-seed tree to the default folder & name: " << treeFile << std::endl;

    std::fstream tree_output(toCString(treeFile), std::ios::binary | std::ios::out);
    // save the adaptive-seed tree ("F", "numKTab", "coreKTab", "tableKTab") to the binary file
	saveTable(F, numKTab, coreKTab, tableKTab, tree_output);

	std::cout << "Tree saved!" << std::endl;
	
	timeCPU = cpuTime() - timeCPU;
	timeSYS = sysTime() - timeSYS;
	std::cout << "Saving time (CPU): " << timeCPU << " seconds" << std::endl;
	std::cout << "Saving time (wall-clock): " << timeSYS << " seconds" << std::endl;

	getrusage(RUSAGE_SELF, &ru);        
    std::cout << "Max RSS: " << ru.ru_maxrss << " KB" << std::endl;
    //std::cout << "Major page faults: " << ru.ru_majflt << std::endl;
    //std::cout << "Minor page faults: " << ru.ru_minflt << std::endl;

    std::cout << std::endl;
    
// */	

	totalCPU = cpuTime() - totalCPU;
	totalSYS = sysTime() - totalSYS;
	std::cout << "Total running time (CPU): " << totalCPU << " seconds" << std::endl;
	std::cout << "Total running time (wall-clock): " << totalSYS << " seconds" << std::endl;
    std::cout << std::endl;
    
    return 0;
}
