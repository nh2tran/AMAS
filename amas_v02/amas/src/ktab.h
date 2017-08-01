// =====================================================================
// nhTran's notes:
// 
// This file contains the data structures of the adaptive-seed tree
// and its nodes, including functions to build, save, and load the tree.
// =====================================================================





















#ifndef NHTRAN_KTAB_H_
#define NHTRAN_KTAB_H_






#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <typeinfo>
#include <exception> 

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "sequence_index.h"






#define K 10 // length of K-mers
#define KOUNT 1048576 // count of K-mers: 4^K






using namespace seqan;





















// ---------------------------------------------------------------------
// Class "KTab1": structure of a node in the adaptive-seed tree
//
// For each node, we store its locations if its "freq" < "F" (leaf nodes);
// otherwise, we store the links to its child nodes.
// ---------------------------------------------------------------------

struct KTab1
{
		// frequency in the genome
		unsigned freq;
		
		// locations in the genome
		unsigned char * contigs_ID;
		unsigned * contigs_Begin;

		// links to child nodes (this is how the tree is formed)
		unsigned * childs_ID;
		
		KTab1(): freq(0), contigs_ID(0), contigs_Begin(0), childs_ID(0) 
		{}
};





















// ---------------------------------------------------------------------
// Function "fillKmers": recursively generate all (4^K) K-mers
//
// They corrrespond to the core nodes of the adaptive-seed tree.
// ---------------------------------------------------------------------

void fillKmers(StringSet<DnaString> & kmers, unsigned & kmer_ID, DnaString curString, unsigned char curLength) 
{// 1st call (empty,0,"",0)
    if(curLength == K) 
    {
        appendValue(kmers,curString);
        kmer_ID ++;
    }
    else 
    {
		DnaString nextString = curString;
        for(unsigned char letter = 0; letter < 4; ++letter) {
			nextString = curString;
			appendValue(nextString, Dna(letter));
            fillKmers(kmers, kmer_ID, nextString, curLength+1);
		}
    }
}






// ---------------------------------------------------------------------
// Function "fillKTab": fill in the content of a node
//
// We store its locations if its "freq" < "F" (leaf nodes);
// otherwise, we store the links to its child nodes.
// ---------------------------------------------------------------------

void fillKTab(TGenomeFM & fmIndex, unsigned const & F, std::vector<KTab1> & tableKTab, unsigned & numKTab, TIter curIter, unsigned curDepth, unsigned & maxDepth, unsigned & numLeaf)
// 1st call: (fmIndex, F, tableKTab, 2, curIter, 1, 0, 0)
{
	// to record the maximum branch length
	if (maxDepth < curDepth) 
	{
		maxDepth = curDepth;
	}
	
	unsigned ktab_ID = numKTab;
	tableKTab.at(ktab_ID).freq = countOccurrences(curIter);

	if (tableKTab.at(ktab_ID).freq > F) // fill childs if exist
	{
		tableKTab.at(ktab_ID).childs_ID = new unsigned[4];
		
		TIter nextIter(fmIndex);
		for (unsigned char letter = 0; letter < 4; ++letter) 
		{
			nextIter = curIter;
			if (goDown(nextIter, Dna(letter))) // only fill reachable childs
			{
				numKTab++;
				tableKTab.at(ktab_ID).childs_ID[letter] = numKTab;
				tableKTab.push_back(KTab1());
				fillKTab(fmIndex, F, tableKTab, numKTab, nextIter, curDepth+1, maxDepth, numLeaf);
			}
			else  // unreachable >> null pointers point to index 0
			{
				tableKTab.at(ktab_ID).childs_ID[letter] = 0;
			}
		}
	}
	else // find and store locations of the leaf nodes
	{
		numLeaf ++;
		
		tableKTab.at(ktab_ID).contigs_ID = new unsigned char[tableKTab.at(ktab_ID).freq];
		tableKTab.at(ktab_ID).contigs_Begin = new unsigned[tableKTab.at(ktab_ID).freq];
	
		Pair<TContigStoreSize, TContigSeqSize> temp_pair;
		for (unsigned occ = 0; occ < tableKTab.at(ktab_ID).freq; ++occ) 
		{
			temp_pair = toSuffixPosition(fmIndex, getOccurrences(curIter)[occ], K+curDepth);
			tableKTab.at(ktab_ID).contigs_ID[occ] = getValueI1(temp_pair);
			tableKTab.at(ktab_ID).contigs_Begin[occ] = getValueI2(temp_pair);
		}
	}
}






// ---------------------------------------------------------------------
// Function "fillTable": fill in "tableKTab" and its nodes
// ---------------------------------------------------------------------

void fillTable(TGenomeFM & fmIndex, unsigned const & F, std::vector<unsigned> & coreKTab, std::vector<KTab1> & tableKTab, unsigned & numKTab, unsigned & maxDepth, unsigned & numLeaf) 
{
// 1st call: (fmIndex, F, coreKTab, tableKTab, 0, 0, 0)
	
	// generate the set of (K=10)-mers, corresponding to the core nodes in "tableKTab"
	StringSet<DnaString> kmers;
	reserve(kmers, KOUNT);
	unsigned kmer_ID = 0;
	fillKmers(kmers, kmer_ID, "", 0);
	
	
	
	
	
	
	// fill in "tableKTab"
	TIter curIter(fmIndex); // the iterator to traverse along the FM index
	unsigned ktab_ID = numKTab;
	
	for (kmer_ID = 0; kmer_ID < KOUNT; kmer_ID++)
	{
		numKTab ++;
		ktab_ID = numKTab;
		coreKTab.at(kmer_ID) = ktab_ID;

		tableKTab.push_back(KTab1());
		
		goRoot(curIter);
		if (goDown(curIter, kmers[kmer_ID])) // "freq" > 0
		{ 
			tableKTab.at(ktab_ID).freq = countOccurrences(curIter);
			
			if (tableKTab.at(ktab_ID).freq > F) // fill childs if exist
			{
				tableKTab.at(ktab_ID).childs_ID =  new unsigned[4];
				
				TIter nextIter(fmIndex);
				for (unsigned char letter = 0; letter < 4; letter++) 
				{
					nextIter = curIter;
					if (goDown(nextIter, Dna(letter))) // only fill reachable childs
					{
						numKTab ++;
						tableKTab.at(ktab_ID).childs_ID[letter] = numKTab;
						tableKTab.push_back(KTab1());
						fillKTab(fmIndex, F, tableKTab, numKTab, nextIter, 1, maxDepth, numLeaf);
					}
					else // unreachable >> null pointers point to index 0
					{
						tableKTab.at(ktab_ID).childs_ID[letter] = 0;
					}
				}
			}
			else // find and store locations of the leaf nodes
			{
				numLeaf ++;
				
				tableKTab.at(ktab_ID).contigs_ID = new unsigned char[tableKTab.at(ktab_ID).freq];
				tableKTab.at(ktab_ID).contigs_Begin = new unsigned[tableKTab.at(ktab_ID).freq];

				Pair<TContigStoreSize, TContigSeqSize> temp_pair;
				for (unsigned occ = 0; occ < tableKTab.at(ktab_ID).freq; occ++) 
				{
					temp_pair = toSuffixPosition(fmIndex, getOccurrences(curIter)[occ], K);
					tableKTab.at(ktab_ID).contigs_ID[occ] = getValueI1(temp_pair);
					tableKTab.at(ktab_ID).contigs_Begin[occ] = getValueI2(temp_pair);
				}
			}

		}
		
		// progress indicator; KOUNT = 10485,76)
		if ((kmer_ID % 104850) == 0)
		{
			std::cout << (unsigned) kmer_ID / 10485 << "\%\t" << std::flush;
		}
	}

	std::cout << std::endl;
}






// ---------------------------------------------------------------------
// Function "saveTable": save the adaptive-seed tree to a binary file
// ---------------------------------------------------------------------

void saveTable(unsigned const & F, unsigned const & numKTab, std::vector<unsigned> & coreKTab, std::vector<KTab1> & tableKTab, std::fstream & tree_output) 
{
	unsigned char size_char = sizeof(unsigned char);
	unsigned char size_int = sizeof(unsigned);
	
	// save "F"
    seqan::streamWriteBlock(tree_output, (char*)&(F), size_int);

	// save "numKTab"
    seqan::streamWriteBlock(tree_output, (char*)&(numKTab), size_int);
	
	// save "coreKTab"
	for (unsigned kmer_ID = 0; kmer_ID < KOUNT; ++kmer_ID)
	{
	    seqan::streamWriteBlock(tree_output, (char*)&(coreKTab.at(kmer_ID)), size_int);
	}

	// save "tableKTab"
	for (unsigned ktab_ID = 1; ktab_ID <= numKTab; ++ktab_ID)
	{
	    seqan::streamWriteBlock(tree_output, (char*)&(tableKTab.at(ktab_ID).freq), size_int);
		
		if (tableKTab.at(ktab_ID).freq > F) // child nodes
		{
			seqan::streamWriteBlock(tree_output, (char*)tableKTab.at(ktab_ID).childs_ID, 4*size_int);
		}
		else if (tableKTab.at(ktab_ID).freq > 0) // leaf nodes
		{
		    seqan::streamWriteBlock(tree_output, (char*)tableKTab.at(ktab_ID).contigs_ID, tableKTab.at(ktab_ID).freq * size_char);
		    seqan::streamWriteBlock(tree_output, (char*)tableKTab.at(ktab_ID).contigs_Begin, tableKTab.at(ktab_ID).freq * size_int);
		}
	}
}






// ---------------------------------------------------------------------
// Function "loadTable": load the adaptive-seed tree
// ---------------------------------------------------------------------

void loadTable(unsigned & F, std::vector<unsigned> & coreKTab, std::vector<KTab1> & tableKTab, std::fstream & tree_input) 
{
	unsigned char size_char = sizeof(unsigned char);
	unsigned char size_int = sizeof(unsigned);
	
	// load "F"
	seqan::streamReadBlock((char*)&(F), tree_input, size_int);

	// load "numKTab"
	unsigned numKTab = 0;
	seqan::streamReadBlock((char*)&(numKTab), tree_input, size_int);
	
	// load "coreKTab"
	for (unsigned kmer_ID = 0; kmer_ID < KOUNT; ++kmer_ID)
	{
		seqan::streamReadBlock((char*)&(coreKTab.at(kmer_ID)), tree_input, size_int);
	}
	
	// load "tableKTab"
	tableKTab.reserve(numKTab+1);
	for (unsigned ktab_ID = 1; ktab_ID <= numKTab; ++ktab_ID)
	{
		tableKTab.push_back(KTab1());

		seqan::streamReadBlock((char*)&(tableKTab.at(ktab_ID).freq), tree_input, size_int);
		
		if (tableKTab.at(ktab_ID).freq > F) // child nodes
		{
			tableKTab.at(ktab_ID).childs_ID = new unsigned[4];
			seqan::streamReadBlock((char*)tableKTab.at(ktab_ID).childs_ID, tree_input, 4*size_int);
		}
		else if (tableKTab.at(ktab_ID).freq > 0) // leaf nodes
		{
			tableKTab.at(ktab_ID).contigs_ID = new unsigned char[tableKTab.at(ktab_ID).freq];
			tableKTab.at(ktab_ID).contigs_Begin = new unsigned[tableKTab.at(ktab_ID).freq];
			
			seqan::streamReadBlock((char*)tableKTab.at(ktab_ID).contigs_ID, tree_input, tableKTab.at(ktab_ID).freq * size_char);
			seqan::streamReadBlock((char*)tableKTab.at(ktab_ID).contigs_Begin, tree_input, tableKTab.at(ktab_ID).freq * size_int);
		}
		
	}
}






// ---------------------------------------------------------------------
// Function "kmer2index": convert a K-mer to its integer ID in "coreKTab"
// ---------------------------------------------------------------------
int kmer2index(Dna5String kmer_String) 
{
	char base4_String[K];
	for (unsigned i = 0; i < K; ++i) {
		if (ordValue(kmer_String[i]) > 3) // not A,C,G,T, but letter "N" spotted, invalid kmer
		{
			return -1;
		}
		else
		{
			base4_String[i] = '0' + ordValue(kmer_String[i]);
		}
	}
	base4_String[K] = '\0';
	char * pEnd;
	return strtol(base4_String, & pEnd,4);
}






#endif
