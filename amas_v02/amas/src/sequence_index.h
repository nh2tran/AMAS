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
// This file contains basic type definitions.
// ==========================================================================






// ==========================================================================
// nhTran's notes:
//
// In AMAS we reuse MASAI's well-developed data structures of 
// the genome sequence, its FM index, and the read sequences
// (including functions for loading the genome file & the reads file).
//
// The code below is part of the two files "store.h" and "index.h" 
// in MASAI package (see above copyright).
// 
// We modified class "Genome" to include the number of contigs,
// their sizes, and the genome's size.
//
// We modified class "Reads" to include the number of seeds, the number
// of hits, the number of matches, the last seed's count and length 
// for each read.
// ==========================================================================





















#ifndef SEQAN_EXTRAS_MASAI_SEQUENCE_INDEX_H_
#define SEQAN_EXTRAS_MASAI_SEQUENCE_INDEX_H_






#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "store/store_io.h"
#include "tags.h"






using namespace seqan;





















struct FragStoreConfig
{
    typedef String<Dna5Q>			TReadSeq;
    typedef String<Dna5>            TContigSeq; 

    typedef double                  TMean;
    typedef double                  TStd;
    typedef signed char             TMappingQuality;

    typedef void                    TReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TReadSeqStoreSpec;
    typedef Alloc<>                 TReadNameSpec;
    typedef Owner<ConcatDirect<> >  TReadNameStoreSpec;
    typedef void                    TMatePairStoreElementSpec;
    typedef void                    TLibraryStoreElementSpec;
    typedef void                    TContigStoreElementSpec;
    typedef void                    TContigFileSpec;
    typedef void                    TAlignedReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TAlignedReadTagStoreSpec;
    typedef void                    TAnnotationStoreElementSpec;
};






typedef StringSet<FragStoreConfig::TContigSeq, Dependent<> > TContigs;

namespace seqan {

    template <>
    struct Size<FragStoreConfig::TContigSeq>
    {
        typedef unsigned int            Type;
    };

    // NOTE(esiragusa): Genome can be at most 2^32 bp in total
    template <>
    struct StringSetLimits<TContigs>
    {
        typedef String<unsigned int>    Type;
    };
}

namespace seqan {
template <>
struct SAValue<TContigs>
{
    typedef Pair<unsigned char, unsigned int, Pack> Type;
};
}






typedef Index<TContigs, FMIndex<WT<>, CompressText> > TGenomeFM;
typedef Iterator<TGenomeFM, TopDown<> >::Type TIter;






// ----------------------------------------------------------------------------
// Fragment Store
// ----------------------------------------------------------------------------

typedef FragmentStore<void, FragStoreConfig> TFragmentStore;






// ----------------------------------------------------------------------------
// Fragment Store Contig Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TContigStore                    TContigStore;
typedef Size<TContigStore>::Type                        TContigStoreSize;
typedef Value<TContigStore>::Type                       TContigStoreElement;
typedef TFragmentStore::TContigSeq                      TContigSeq;
typedef Size<TContigSeq>::Type                          TContigSeqSize;
typedef Segment<TContigSeq, InfixSegment>               TContigInfix;






// ----------------------------------------------------------------------------
// Fragment Store Reads Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TReadStore                      TReadStore;
typedef Value<TReadStore>::Type                         TReadStoreElement;
typedef TFragmentStore::TReadNameStore                  TReadNameStore;
typedef TFragmentStore::TReadSeqStore                   TReadSeqStore;
typedef Size<TReadSeqStore>::Type                       TReadSeqStoreSize;
typedef Value<TReadSeqStore>::Type const                TReadSeq;
typedef Size<TReadSeq>::Type                            TReadSeqSize;






// ----------------------------------------------------------------------------
// Fragment Store Mapped Reads Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TAlignedReadStore               TAlignedReadStore;
typedef Value<TAlignedReadStore>::Type                  TAlignedReadStoreElement;
typedef TFragmentStore::TAlignQualityStore              TAlignQualityStore;
typedef Value<TAlignQualityStore>::Type                 TAlignQualityStoreElement;
typedef TFragmentStore::TAlignedReadTagStore            TAlignedReadTagStore;
typedef Value<TAlignedReadTagStore>::Type               TAlignedReadTagStoreElement;

typedef unsigned long 									TWord;
const unsigned 											TWordLen = BitsPerValue<TWord>::VALUE;

//typedef std::vector<BamAlignmentRecord> 				TAlignmentRecords;





// ----------------------------------------------------------------------------
// Class Genome
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Genome
{
    TFragmentStore          & _store;
    TContigs                contigs;
    String<TContigSeqSize>	contigs_Size;
    TContigStoreSize		numContigs;
    unsigned long 			totalLength;

    Genome(TFragmentStore & store) :
        _store(store),
        numContigs(0),
        totalLength(0)
    {}
};






// ---------------------------------------------------------------------
// Class ReadsConfig
// ---------------------------------------------------------------------

template <typename TUseReadStore_       = True,
          typename TUseReadNameStore_   = True,
          typename TForward_            = True,
          typename TReverse_            = True>
struct ReadsConfig
{
    typedef TUseReadStore_      TUseReadStore;
    typedef TUseReadNameStore_  TUseReadNameStore;
    typedef TForward_           TForward;
    typedef TReverse_           TReverse;
};






// ----------------------------------------------------------------------------
// Class Reads
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct Reads
{
    Holder<TFragmentStore>          _store;
    unsigned                        _avgSeqLengthEstimate;
    unsigned                        _avgNameLengthEstimate;
    unsigned                        _countEstimate;
    unsigned                        readsCount;
    
    std::vector<unsigned> 		matchesCount;
    //std::vector<unsigned> 		seedsCount;
    std::vector<unsigned> 		hitsCount;

    //std::vector<unsigned> 		lastSeed_Length;
    //std::vector<unsigned> 		lastSeed_Count;
    
    String<TWord> 				_readAligned;
    
    //std::vector<TAlignmentRecords> readsAlignments;
    
    Reads() :
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
        readsCount(0)
    {}

    Reads(TFragmentStore & store) :
        _store(store),
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
        readsCount(0)
    {}
};






// ----------------------------------------------------------------------------
// Class ReadsLoader
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct ReadsLoader
{
    // TODO(esiragusa): Use MMap FileReader. Support compressed streams.
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;
    typedef Reads<TSpec, TConfig>                   TReads;

    TStream                                 _file;
    unsigned long                           _fileSize;
    AutoSeqStreamFormat                     _fileFormat;
    std::SEQAN_AUTO_PTR_NAME<TRecordReader> _reader;
    Holder<TReads>                          reads;

    ReadsLoader(TReads & reads) :
        _fileSize(0),
        reads(reads)
    {}
};





















// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function loadGenome()                                               [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool loadGenome(Genome<TSpec> & genome, TString const & genomeFile)
{
    if (!loadContigs(genome._store, genomeFile)) {
		return false;
	}
	
	genome.numContigs = length(genome._store.contigStore);
	
    reserve(genome.contigs, genome.numContigs);
    reserve(genome.contigs_Size, genome.numContigs);
    for (TContigStoreSize contig_ID = 0; contig_ID < genome.numContigs; ++contig_ID) {
		shrinkToFit(genome._store.contigStore[contig_ID].seq);
        appendValue(genome.contigs, genome._store.contigStore[contig_ID].seq);
        appendValue(genome.contigs_Size, length(genome.contigs[contig_ID]));
        //std::cout << "Sequence " << contig_ID << " has " << genome.contigs_Size[contig_ID] << " nucleotides" << std::endl;
        genome.totalLength += genome.contigs_Size[contig_ID];
	}

    return true;
}






// ----------------------------------------------------------------------------
// Function open()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
bool open(ReadsLoader<TSpec, TConfig> & readsLoader, TString const & readsFile)
{
    typedef ReadsLoader<TSpec, TConfig>             TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    readsLoader._file.open(toCString(readsFile), std::ios::binary | std::ios::in);

    if (!readsLoader._file.is_open())
        return false;

    readsLoader._file.seekg(0, std::ios::end);
    readsLoader._fileSize = readsLoader._file.tellg();
    readsLoader._file.seekg(0, std::ios::beg);

    // Initialize record reader.
    readsLoader._reader.reset(new TRecordReader(readsLoader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(readsLoader._reader), readsLoader._fileFormat))
        return false;

    _estimateReadsStatistics(readsLoader);

    // Reinitialize record reader.
    readsLoader._file.seekg(0, std::ios::beg);
    readsLoader._reader.reset(new TRecordReader(readsLoader._file));

    return true;
}






// ----------------------------------------------------------------------------
// Function close()                                               [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
bool close(ReadsLoader<TSpec, TConfig> & readsLoader)
{
    readsLoader._file.close();

    return !readsLoader._file.is_open();
}






// ----------------------------------------------------------------------------
// Function _estimateReadsStatistics()                            [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _estimateReadsStatistics(ReadsLoader<TSpec, TConfig> & readsLoader)
{
    CharString                  seqName;
    FragStoreConfig::TReadSeq   seq;

    // Read first record.
    if (readRecord(seqName, seq, *(readsLoader._reader), readsLoader._fileFormat) != 0)
        return;

    value(readsLoader.reads)._avgSeqLengthEstimate = length(seq);
    value(readsLoader.reads)._avgNameLengthEstimate = length(seqName);

    // Estimate record size.
    unsigned long recordSize;
    switch (readsLoader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        recordSize = _estimateRecordSize(readsLoader, Fasta());
        break;
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        recordSize = _estimateRecordSize(readsLoader, Fastq());
        break;
    default:
        recordSize = 0;
        break;
    }

    // Estimate number of reads in file.
    if (recordSize > 0)
        value(readsLoader.reads)._countEstimate = readsLoader._fileSize / recordSize;
    else
        value(readsLoader.reads)._countEstimate = 0;
}






// ----------------------------------------------------------------------------
// Function _estimateRecordSize()                                 [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & readsLoader, Fastq const & /* tag */)
{
    // 6 stands for: @, +, and four \n.
    return value(readsLoader.reads)._avgNameLengthEstimate + 2 * value(readsLoader.reads)._avgSeqLengthEstimate + 6;
}

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & readsLoader, Fasta const & /* tag */)
{
    // 3 stands for: >, and two \n.
    return value(readsLoader.reads)._avgNameLengthEstimate + value(readsLoader.reads)._avgSeqLengthEstimate + 3;
}






// ----------------------------------------------------------------------------
// Function load()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
bool load(ReadsLoader<TSpec, TConfig> & readsLoader, TSize count)
{
    switch (readsLoader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        return load(readsLoader, count, Fasta());
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        return load(readsLoader, count, Fastq());
    default:
        return false;
    }
}

template <typename TSpec, typename TConfig, typename TSize, typename TFormat>
bool load(ReadsLoader<TSpec, TConfig> & readsLoader, TSize count, TFormat const & /* tag */)
{
    CharString                  seqName;
    FragStoreConfig::TReadSeq   seq;

    for (; count > 0 && !atEnd(readsLoader); count--)
    {
        // NOTE(esiragusa): AutoFormat seems to thrash memory allocation.
        //if (readRecord(seqName, seq, *(loader._reader), loader._fileFormat) != 0)

        if (readRecord(seqName, seq, *(readsLoader._reader), TFormat()) != 0)
            return false;

        appendSeq(value(readsLoader.reads), seq);
        appendName(value(readsLoader.reads), seqName, typename TConfig::TUseReadNameStore());
        appendId(value(readsLoader.reads), TReadStoreElement::INVALID_ID, typename TConfig::TUseReadStore());
    }

    value(readsLoader.reads).readsCount = length(getSeqs(value(readsLoader.reads)));
    value(readsLoader.reads).matchesCount.resize(2*value(readsLoader.reads).readsCount,0);
    //value(readsLoader.reads).seedsCount.resize(2*value(readsLoader.reads).readsCount,0);
    value(readsLoader.reads).hitsCount.resize(2*value(readsLoader.reads).readsCount,0);
    //value(readsLoader.reads).lastSeed_Length.resize(2*value(readsLoader.reads).readsCount,0);
    //value(readsLoader.reads).lastSeed_Count.resize(2*value(readsLoader.reads).readsCount,0);

    resize(value(readsLoader.reads)._readAligned, (value(readsLoader.reads).readsCount + TWordLen - 1) / TWordLen, 0);
    
    //value(readsLoader.reads).readsAlignments.resize(value(readsLoader.reads).readsCount);

    _loadReverseComplement(readsLoader);

    return true;
}






// ----------------------------------------------------------------------------
// Function atEnd()                                               [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline bool atEnd(ReadsLoader<TSpec, TConfig> & readsLoader)
{
    return atEnd(*(readsLoader._reader));
}






// ----------------------------------------------------------------------------
// Function _loadReverseComplement()                              [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _loadReverseComplement(ReadsLoader<TSpec, TConfig> & readsLoader)
{
    for (TReadSeqStoreSize readId = 0; readId < value(readsLoader.reads).readsCount; ++readId)
    {
        TReadSeq & read = getSeqs(value(readsLoader.reads))[readId];
        appendSeq(value(readsLoader.reads), read);
        reverseComplement(back(getSeqs(value(readsLoader.reads))));
    }
}






// ----------------------------------------------------------------------------
// Function reserve()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(Reads<TSpec, TConfig> & reads, TSize count)
{
    reserve(getSeqs(reads).concat, 2 * count * reads._avgSeqLengthEstimate, Exact());
    reserve(getSeqs(reads), 2 * count, Exact());

    reserveIds(reads, count, typename TConfig::TUseReadStore());

    reserveNames(reads, count, reads._avgNameLengthEstimate, typename TConfig::TUseReadNameStore());
}






// ----------------------------------------------------------------------------
// Function reserveIds()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & /* reads */, TSize /* space */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & reads, TSize count, True const & /* tag */)
{
    reserve(getIds(reads), count, Exact());
}






// ----------------------------------------------------------------------------
// Function reserveNames()                                              [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & /* reads */, TSize /* count */, TLength /* length */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & reads, TSize count, TLength length, True const & /* tag */)
{
    reserve(getNames(reads).concat, count * length, Exact());
    reserve(getNames(reads), count, Exact());
}






// ----------------------------------------------------------------------------
// Function appendSeq()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeq>
inline void appendSeq(Reads<TSpec, TConfig> & reads, TReadSeq const & seq)
{
    appendValue(getSeqs(reads), seq, Generous());
}






// ----------------------------------------------------------------------------
// Function appendName()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & /* reads */, TReadName const & /* seqName */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & reads, TReadName const & seqName, True const & /* tag */)
{
    appendValue(getNames(reads), seqName, Generous());
}






// ----------------------------------------------------------------------------
// Function appendId()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & /* reads */, TReadId const & /* matePairId */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & reads, TReadId const & matePairId, True const & /* tag */)
{
	typename Value<TFragmentStore::TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(getIds(reads), r, Generous());
}






// ----------------------------------------------------------------------------
// Function getSeqs()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadSeqStore & getSeqs(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readSeqStore;
}






// ----------------------------------------------------------------------------
// Function getNames()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadNameStore & getNames(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readNameStore;
}






// ----------------------------------------------------------------------------
// Function getIds()                                                    [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadStore & getIds(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readStore;
}






// ----------------------------------------------------------------------------
// Function clear()                                                     [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Reads<TSpec, TConfig> & reads)
{
    clearReads(value(reads._store));

    reads.matchesCount.clear();
    //reads.seedsCount.clear();
    reads.hitsCount.clear();
    //reads.lastSeed_Length.clear();
    //reads.lastSeed_Count.clear();
    
    clear(reads._readAligned);
    
/*    for (TReadSeqStoreSize read_ID = 0; read_ID < reads.readsCount; ++read_ID)
    {
		reads.readsAlignments.at(read_ID).clear();
	}
	reads.readsAlignments.clear();
*/	
}






#endif
