// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Mapper class.
// ==========================================================================






// ==========================================================================
// nhTran's notes:
//
// To get compatible with the components "writer" and "extender" of MASAI,
// we reused classes "Mapper" and "ReadMapperConfig".
//
// The class definitions below were copied from the file "mapper.h" 
// in MASAI package (see above copyright).
// 
// We modified "ReadMapperConfig" to consider only EditDistance & All-mapping.
//
// A new function "mapReads" was provided 
// since AMAS uses a different seeding technique from MASAI.
// ==========================================================================






















#ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
#define SEQAN_EXTRAS_MASAI_MAPPER_H_






#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

#include "ktab.h"

#include "tags.h"
#include "sequence_index.h"
#include "seeder.h"
#include "extender.h"
#include "manager.h"
#include "parser_options.h"






using namespace seqan;





















// ----------------------------------------------------------------------------
// Class ReadMapperConfig
// ----------------------------------------------------------------------------

template <typename TDistance_       = EditDistance,
          /*typename TStrategy_       = AnyBest,*/
          typename TStrategy_       = All
          /*,typename TBacktracking_   = MultipleBacktracking*/>
struct ReadMapperConfig
{
	typedef TDistance_          TDistance;
	typedef TStrategy_          TStrategy;
    /*typedef TBacktracking_      TBacktracking;*/
};






// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TWriter, typename TMapperConfig>
struct Mapper
{
    //typedef MatchManager<TWriter, typename TMapperConfig::TStrategy> 				TManager;
    typedef Extender<TGenome, TReads, TWriter /*TManager*/, typename TMapperConfig::TDistance> 	TExtender;
    typedef Seeder<TGenome, TReads, /*TManager,*/ TExtender> 						TSeeder;

    TFragmentStore 	& store;
    TGenome 		& genome;
    TReads 			& reads;
    
    //TManager 		manager;
    TExtender 		extender;
    TSeeder 		seeder;
    
    Mapper(TFragmentStore & store, TGenome & genome, TReads & reads, TWriter & writer, bool disableExtender = false) :
        store(store),
        genome(genome),
        reads(reads),
        //manager(writer, reads.readsCount),
        extender(store, genome, reads, writer /*manager*/, disableExtender),
        seeder(store, genome, reads, extender)
	{}
};





















// ----------------------------------------------------------------------------
// Function mapReads()                                                 [Mapper]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TWriter, typename TMapperConfig>
void mapReads(Mapper<TGenome, TReads, TWriter, TMapperConfig> const & mapper,
			unsigned F,
			std::vector<unsigned> const & coreKTab, 
			std::vector<KTab1> const & tableKTab,
			Options const & options)
{
	//unsigned readsCountPercent = mapper.reads.readsCount / 100;
	// Seed Partition & Extension
	//std::cout << "testMap" << std::endl;
	#pragma omp parallel
	{
		#pragma omp for schedule(static, 1)
		for (int read_ID = 0; read_ID < mapper.reads.readsCount; ++read_ID)
	    {
			//std::cout << "read_ID" << "\t" << read_ID << std::endl;
			TReadSeqStoreSize read_revcom_ID = read_ID + mapper.reads.readsCount;
			seedsPartition(mapper.seeder, read_ID, F, coreKTab, tableKTab, options); // read
			seedsPartition(mapper.seeder, read_revcom_ID, F, coreKTab, tableKTab, options); // & its reverse complement
			
/*		// progress indicator
			if ((read_ID % (10*readsCountPercent)) == 0)
			{
				std::cout << (unsigned) read_ID / readsCountPercent << "\%\t" << std::flush;
			}
*/
		}
	}
}






#endif  // #ifndef SEQAN_EXTRAS_MASAI_MAPPER_H_
