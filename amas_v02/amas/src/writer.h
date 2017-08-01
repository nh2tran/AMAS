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
// This file contains the Writer class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_WRITER_H_
#define SEQAN_EXTRAS_MASAI_WRITER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "tags.h"
#include "sequence_index.h"
#include "matches.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Writer
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec = void>
struct Writer {};






template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
struct Writer<TGenome, TReads, Sam, TDistance, TSpec>
{
    typedef BamIOContext<TFragmentStore::TContigNameStore>  TBamIOContext;
    //typedef unsigned long                                   TWord; //nhtran: moved to "sequence_index.h"
    typedef Stream<FileStream<File<> > >         TStream;

    TFragmentStore          & store;
    TGenome                 & genome;
    TReads                  & reads;
    TStream                 _stream;
    TBamIOContext           _context;
    bool                    disabled;
    bool                    _writeCigar;
    //String<TWord>           _readAligned; //nhtran: moved to object "reads"
    //const unsigned          _wordLen; //nhtran: moved to "sequence_index.h"

    Writer(TFragmentStore & store, TGenome & genome, TReads & reads, bool disabled = false) :
        store(store),
        genome(genome),
        reads(reads),
        _context(store.contigNameStore, store.contigNameStoreCache),
        disabled(disabled),
        _writeCigar(true)
        //_wordLen(BitsPerValue<TWord>::VALUE) //nhtran: moved to "sequence_index.h"
    {}
};






template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
struct Writer<TGenome, TReads, Raw, TDistance, TSpec>
{
    typedef Match<>								TMatch;
    typedef Stream<FileStream<File<>, TMatch> >	TStream;

    TFragmentStore          & store;
    TGenome                 & genome;
    TReads                  & reads;
    TStream                 _stream;
    bool                    disabled;

    Writer(TFragmentStore & store, TGenome & genome, TReads & reads, bool disabled = false) :
        store(store),
        genome(genome),
        reads(reads),
        disabled(disabled)
    {}
};






// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function onMatch()                                             [Writer<Sam>]
// ----------------------------------------------------------------------------

// Single-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Raw, TDistance, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
	if (writer.disabled)
        return;

    // Fill record.
    Match<> match;
    fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Write record.
    #pragma omp critical
    {
	    streamWriteChar(writer._stream, match);
	}		    
    return;
}






template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
	if (writer.disabled)
        return;

    typedef Align<TFragmentStore::TReadSeq, ArrayGaps>  TAlign;

    TAlignedReadStoreElement        alignedRead;
    TAlignedReadStoreElement        alignedMate;
    TAlignQualityStoreElement       alignQuality;
    TAlignedReadTagStoreElement     alignedTags;

    // Fill aligned read.
    _fillAlignedRead(alignedRead, alignQuality,
                     contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Check for secondary alignment.
    bool secondary = _checkSecondary(writer, alignedRead);

    // Align read only if cigar output is enabled.
    TAlign align;
    if (writer._writeCigar)
        _alignRead(writer, align, alignedRead, alignQuality, reverseComplemented);

    // Write aligned read.
//    _writeAlignedRead(writer.store, writer._stream, writer.context,
//                      alignedRead, alignQuality, alignedTags,
//                      alignedMate, align, secondary, Sam());

    // Fill record.
    BamAlignmentRecord record;
    _fillRecord(writer.store, record, alignedRead, alignQuality, alignedTags,
                alignedMate, align, secondary, writer._writeCigar);
    
    //writer.reads.readsAlignments.at(readId).push_back(record);

    // Write record to target.
    #pragma omp critical
    {
		write2(writer._stream, record, writer._context, Sam());
	}
}











// ----------------------------------------------------------------------------
// Function _fillAlignedRead()
// ----------------------------------------------------------------------------

template <typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void _fillAlignedRead(TAlignedReadStoreElement & alignedRead,
                             TAlignQualityStoreElement & alignQuality,
                             TContigId contigId,
                             TContigPos beginPos,
                             TContigPos endPos,
                             TReadId readId,
                             TErrors errors,
                             bool reverseComplemented)
{
    alignedRead.readId   = readId;
    alignedRead.contigId = contigId;

    if (reverseComplemented)
    {
        alignedRead.beginPos = endPos;
        alignedRead.endPos   = beginPos;
    }
    else
    {
        alignedRead.beginPos = beginPos;
        alignedRead.endPos   = endPos;
    }

    alignQuality.errors = errors;
}











// ----------------------------------------------------------------------------
// Function _checkSecondary()                                     [Writer<Sam>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
inline bool _checkSecondary(Writer<TGenome, TReads, Sam, TDistance, TSpec> const & writer,
                            TAlignedReadStoreElement & alignedRead)
{
    TWord mask     = (TWord)1 << (alignedRead.readId % TWordLen);
    bool secondary = (writer.reads._readAligned[alignedRead.readId / TWordLen] & mask) != 0;

    writer.reads._readAligned[alignedRead.readId / TWordLen] |= mask;

    return secondary;
}











// ----------------------------------------------------------------------------
// Function _alignRead()                                 [Writer<EditDistance>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, TFormat, EditDistance, TSpec> const & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement & alignQuality,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer.store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += writer.reads.readsCount;

    TReadSeq const & readSeq = writer.store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);

    // In this case no indels are possible.
    if ((alignQuality.errors <= 1) && (length(row(align, 0)) == length(row(align, 1))))
        return;

    globalAlignment(align, Score<short, EditDistance>(),
                    (short)-alignQuality.errors, (short)alignQuality.errors,
                    NeedlemanWunsch());
}











// ----------------------------------------------------------------------------
// Function open()                                                     [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec, typename TString>
bool open(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & writer, TString const & fileName)
{
    if (writer.disabled)
        return true;

    if (!open(writer._stream, toCString(fileName), OPEN_RDWR | OPEN_CREATE))
        return false;
        
    _writeHeader(writer);

    return true;
}











// ----------------------------------------------------------------------------
// Function _writeHeader()                                             [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void _writeHeader(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer)
{
    _writeHeader(writer.store, writer._stream, writer._context, Sam());
}






template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void _writeHeader(Writer<TGenome, TReads, Raw, TDistance, TSpec> & writer)
{
}











// ----------------------------------------------------------------------------
// Function close()                                                    [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
bool close(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & writer)
{
    if (writer.disabled)
        return true;

    return close(writer._stream);
}











// ----------------------------------------------------------------------------
// Function writeAlignments()                                          [Writer]
// ----------------------------------------------------------------------------
/*
template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void writeAlignments(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer)
{
    if (writer.disabled)
        return;

	//#pragma omp parallel for
	for (TReadSeqStoreSize read_ID = 0; read_ID < writer.reads.readsCount; ++read_ID)
	{
		//#pragma omp critical
		{
			for (TAlignmentRecords::iterator it = writer.reads.readsAlignments.at(read_ID).begin(); it != writer.reads.readsAlignments.at(read_ID).end(); ++it)
			{
				write2(writer._stream, *it, writer._context, Sam());
			}
		}
	}
}
*/





#endif  // #ifndef SEQAN_EXTRAS_MASAI_WRITER_H_
