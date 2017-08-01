// =====================================================================
// nhTran's notes:
//
// This file contains functions to perform the partition of adaptive seeds,
// followed by the filtration using last seeds & extra seeds.
// =====================================================================






















#ifndef NHTRAN_SEEDER_H_
#define NHTRAN_SEEDER_H_






#include <vector>       // std::vector
#include <cmath>        // std::abs

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "sequence_index.h"
#include "ktab.h"
#include "parser_options.h"







using namespace seqan;





















// ---------------------------------------------------------------------
// Class "Seeder"
// ---------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TExtender>
struct Seeder
{
    TFragmentStore 	& store;
    TGenome			& genome;
    TReads 			& reads;
    TExtender 		& extender;

    Seeder(TFragmentStore & store, TGenome & genome, TReads & reads, TExtender & extender) :
        store(store),
        genome(genome),
        reads(reads),
        extender(extender) 
	{}
};






// ---------------------------------------------------------------------
// Class "Candidate": to store candidate locations using a binary search tree
//
// The binary search tree allows to quickly sort the locations 
// and identify those reported by multi seeds.
// ---------------------------------------------------------------------

struct Candidate
{
		unsigned char contig_ID;
		unsigned contig_Begin; // location of the seed in the contig
		int contig_Loc; // location of the read in the contig
		unsigned char seed; // seed's ID in the read
		
		unsigned hitsCount; // number of seeds reporting this  candidate

		unsigned left; // links to other candidates in the tree
		unsigned right;
		
		Candidate(): 
			contig_ID(0), 
			contig_Begin(0), 
			contig_Loc(0), 
			seed(0),
			hitsCount(0), 
			left(0), 
			right(0)
		{}

};





















// ------------------------------------------------------------------------------------------
// Function "insertCandidate": recursively insert a candidate location to the tree "tableCan"
// ------------------------------------------------------------------------------------------

void insertCandidate(std::vector<Candidate> & tableCan, 
					unsigned curCan_ID, 
					Candidate const & newCan, 
					unsigned maxErrors)
{
	if (std::abs(newCan.contig_Loc - tableCan.at(curCan_ID).contig_Loc) <= maxErrors)
	{
		// duplicate hit
		tableCan.at(curCan_ID).hitsCount++;
	}
	else if (newCan.contig_Loc < tableCan.at(curCan_ID).contig_Loc)
	{
		if (tableCan.at(curCan_ID).left == 0)
		{
			// empty!!! Add the new candidate here
			tableCan.push_back(newCan);
			tableCan.at(curCan_ID).left = tableCan.size() - 1;
		}
		else
		{
			// extend left
			insertCandidate(tableCan, tableCan.at(curCan_ID).left, newCan, maxErrors);
		}
	}
	else 
	{
		if (tableCan.at(curCan_ID).right == 0)
		{
			// empty!!! Add the new candidate here
			tableCan.push_back(newCan);
			tableCan.at(curCan_ID).right = tableCan.size() - 1;
		}
		else
		{
			// extend right
			insertCandidate(tableCan, tableCan.at(curCan_ID).right, newCan, maxErrors);
		}
	}
}






// ---------------------------------------------------------------------
// Function "searchKTab": recursively search "tableKTab" for locations 
// of a seed and insert those locations to "tableCan"
// ---------------------------------------------------------------------

void searchKTab(std::vector<KTab1> const & tableKTab, unsigned ktab_ID, unsigned F, 
				std::vector<Candidate> & tableCan, 
				TReadSeqSize seed_Begin, 
				unsigned char seed, 
				unsigned maxErrors) 
{
	if (tableKTab.at(ktab_ID).freq > F) // not reach a leaf node yet
	{
		for (unsigned char letter=0; letter < 4; ++letter)
		{
			searchKTab(tableKTab, tableKTab.at(ktab_ID).childs_ID[letter], F, tableCan, seed_Begin, seed, maxErrors);
		}
	}
	else // reach a leaf node, extract its locations
	{
		for (unsigned hit = 0; hit < tableKTab.at(ktab_ID).freq; ++hit)
		{
			Candidate newCan;
			newCan.contig_ID = tableKTab.at(ktab_ID).contigs_ID[hit];
			newCan.contig_Begin = tableKTab.at(ktab_ID).contigs_Begin[hit];
			newCan.contig_Loc = newCan.contig_Begin - seed_Begin;
			newCan.seed = seed;
			newCan.hitsCount = 1;
			
			unsigned newCan_initID = newCan.contig_ID + 1;
			insertCandidate(tableCan, newCan_initID, newCan, maxErrors);
		}
	}
}






// ---------------------------------------------------------------------
// Function "seedsPartition": perform the partition of adaptive seeds
// followed by the filtration using last seeds & extra seeds.
// ---------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TExtender>
void seedsPartition(Seeder<TGenome, TReads, TExtender> const & seeder,
			TReadSeqStoreSize read_ID,
			unsigned F,
			std::vector<unsigned> const & coreKTab, 
			std::vector<KTab1> const & tableKTab,
			Options const & options)
{
	unsigned maxErrors = options.errors;
	unsigned fse = options.fse;
	
	TReadSeq & read = seeder.store.readSeqStore[read_ID];
	TReadSeqSize read_Length = length(read);

	// contraint on the maximum seed length
	// we want full sensitivity for "fse" errors >> "fse"+1 seeds required >> max_seed_length = read_Length/("fse"+1)
	TReadSeqSize max_seed_length = read_Length/(fse+1);
	
	// store info of seeds
	std::vector<unsigned> seeds_KTab; // ID in "tableKTab"
	std::vector<TReadSeqSize> seeds_Begin; // location in the read
	std::vector<TReadSeqSize> seeds_Length;
	
	unsigned char numSeeds = 0;
	unsigned numHits = 0;
	TReadSeqSize nextBase = 0;
	TReadSeqSize curBase = nextBase;
	Dna5String kmer_String;
	int kmer_ID = 0; // this ID could be -1 for invalid kmers with letters "N"
	unsigned ktab_ID = 0;
	





//======================================================================
// ADAPTIVE SEED PARTITION
//======================================================================

	while ((nextBase < read_Length-K+1) && (numSeeds < maxErrors + 2)) // pigeonhole principle + 1 extra seed
	{
		curBase = nextBase;
		nextBase += K;
		
		kmer_String = infix(read, curBase, nextBase);
		kmer_ID = kmer2index(kmer_String);
		if (kmer_ID >= 0) // valid kmer
		{
			ktab_ID = coreKTab.at(kmer_ID);
			
			// traverse "tableKTab" with "F" & "max_seed_length" constraints
			while ((tableKTab.at(ktab_ID).freq > F) && (nextBase < read_Length) && (nextBase < (numSeeds+1)*max_seed_length))
			{
				if (ordValue(read[nextBase]) > 3) // not A, C, G, T, but letter "N" spotted
				{
					ktab_ID = 0;
				}
				else
				{
					ktab_ID = tableKTab.at(ktab_ID).childs_ID[ordValue(read[nextBase])];
				}
				nextBase++;
			}
		}
		else
		{
			ktab_ID = 0;
		}
		
		// 3 possible situations: tableKTab.at(ktab_ID).freq == 0, <= F, > F)
		// add ALL seeds, including the last seed, even if its "freq" > "F"
		numSeeds++;
		numHits += tableKTab.at(ktab_ID).freq;
		seeds_KTab.push_back(ktab_ID);
		seeds_Begin.push_back(curBase);
		seeds_Length.push_back(nextBase-curBase);
	}
	
	// only for testing
	//seeder.reads.seedsCount.at(read_ID) = numSeeds;
	
	
	
	
	

//======================================================================
// LAST-SEED filteration
//======================================================================

	unsigned lastSeed_Length = seeds_Length.at(numSeeds-1);
	unsigned lastSeed_Count = tableKTab.at(seeds_KTab.at(numSeeds-1)).freq;
	
	// only for testing
	//seeder.reads.lastSeed_Length.at(read_ID) = lastSeed_Length;
	//seeder.reads.lastSeed_Count.at(read_ID) = lastSeed_Count;
	
	unsigned char numSeeds_filtered = numSeeds;
	// filteration criteria
	if (((double)lastSeed_Count/numHits > 0.95)
		&& (lastSeed_Count > seeder.genome.totalLength/KOUNT) // greater than average frequency of 10-mers
		//&& (lastSeed_Count > 3000) // HUMAN: 3 x 10^9
		//&& (lastSeed_Count > 100) // WORM: 100 x 10^6
		//&& (lastSeed_Count > 169) // FLY: 169 x 10^6
		&& (numSeeds > (fse+1)))
	{
		numSeeds_filtered--;
	}






//======================================================================
// EXTRA-SEED filteration
//======================================================================

	std::vector<Candidate> tableCan;
	
	 // empty element at index 0
	tableCan.push_back(Candidate());
	
	// the next "numContigs" elements are the ROOTS of the binary search trees
	unsigned char numContigs = seeder.genome.numContigs;
	for (unsigned char id = 1; id <= numContigs; id++)
	{
		tableCan.push_back(Candidate());
		tableCan.at(id).contig_Loc = seeder.genome.contigs_Size[id-1] + maxErrors + 10; // assign the largest location to the root nodes
	}
	
	// fill up the tree "tableCan"
	for (unsigned char seed = 0; seed < numSeeds_filtered; seed++)
	{
		ktab_ID = seeds_KTab.at(seed);
		
		// find all hits of the seed & add them to the trees
		searchKTab(tableKTab, ktab_ID, F, tableCan, seeds_Begin.at(seed), seed, maxErrors);
	}

	





//======================================================================
// SEED EXTENSION
//======================================================================

	unsigned char temp = 1; // min hit count to extend a candidate location
	if (numSeeds_filtered > maxErrors)
	{
		temp = numSeeds_filtered - maxErrors;
	}
	
	unsigned numCans = tableCan.size();
	
	if (options.best)
	{
		unsigned hitsCount_max = tableCan.at(0).hitsCount;
		unsigned hitsCount_can = 0;
		for (unsigned can_ID = 1; can_ID < numCans; can_ID++)
		{
			if (tableCan.at(can_ID).hitsCount > hitsCount_max)
			{
				hitsCount_max = tableCan.at(can_ID).hitsCount;
				hitsCount_can = can_ID;
			}
		}
		seeder.reads.hitsCount.at(read_ID) ++;
		onSeedHit(seeder.extender,
						read_ID,
						tableCan.at(hitsCount_can).contig_ID,
						tableCan.at(hitsCount_can).contig_Begin,
						seeds_Begin.at(tableCan.at(hitsCount_can).seed),
						seeds_Length.at(tableCan.at(hitsCount_can).seed),
						maxErrors,
						tableCan.at(hitsCount_can).seed);
		
	}
	else
	{
		for (unsigned can_ID = 0; can_ID < numCans; can_ID++)
		{
			if (seeder.reads.matchesCount.at(read_ID) >= options.k)
				break;
			
			if (tableCan.at(can_ID).hitsCount >= temp)
			{
				seeder.reads.hitsCount.at(read_ID) ++;
				onSeedHit(seeder.extender,
								read_ID,
								tableCan.at(can_ID).contig_ID,
								tableCan.at(can_ID).contig_Begin,
								seeds_Begin.at(tableCan.at(can_ID).seed),
								seeds_Length.at(tableCan.at(can_ID).seed),
								maxErrors,
								tableCan.at(can_ID).seed);
			}
		}
	}
}

#endif  // #ifndef NHTRAN_SEEDER_H_
