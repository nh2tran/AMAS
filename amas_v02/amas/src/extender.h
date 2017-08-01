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
// This file contains the Extender class.
// ==========================================================================






// ==========================================================================
// nhTran's notes:
//
// The code below is part of the file "extender.h" in MASAI package (see above copyright).
// 
// We only use the extension for EditDistance.
//
// We slightly modified a few lines to fit with adaptive seeds.
// ==========================================================================





















#ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
#define SEQAN_EXTRAS_MASAI_EXTENDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

#include "sequence_index.h"

using namespace seqan;





















// ============================================================================
// Types for Myers bit-vector algorithm of alignment
// ============================================================================
    typedef Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;

    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;

    typedef PatternState_<TReadInfix, TAlgorithmSpec>       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithmSpec>    TPatternStateRev;











// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TWriter /*TManager*/, typename TDistance/* = HammingDistance, typename TSpec = void*/>
struct Extender
{
    TFragmentStore 	& store;
    TGenome 		& genome;
    TReads 			& reads;
    //TManager 		& manager;
    TWriter 		& writer;

    bool 			disabled;

    Extender(TFragmentStore & store, TGenome & genome, TReads & reads, TWriter & writer /*TManager & manager*/, bool disabled = false) :
        store(store),
        genome(genome),
        reads(reads),
        //manager(manager),
        writer(writer),
        disabled(disabled)
    {}

};






/*
template <typename TReads, typename TGenome, typename TMatchesDelegate, typename TSpec>
struct Extender<TReads, TGenome, TMatchesDelegate, EditDistance, TSpec>:
    public Extender<TReads, TGenome, TMatchesDelegate>
{
    typedef Extender<TReads, TGenome, TMatchesDelegate>                      TBase;
    typedef Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;

    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;

    typedef PatternState_<TReadInfix, TAlgorithmSpec>       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithmSpec>    TPatternStateRev;

    TPatternState patternState;
    TPatternStateRev patternStateRev;

    Extender(TFragmentStore & store, 
			TReads & reads,
			TGenome & genome,
			TMatchesDelegate & matchesDelegate,
			bool disabled = false) :
        TBase(store, reads, genome, matchesDelegate, disabled)
    {}
};
*/






// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function onSeedHit() - EditDistance version                        [Extender]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TWriter /*TManager*/>
inline bool onSeedHit(Extender<TGenome, TReads, TWriter /*TManager*/, EditDistance> const & extender,
                      TReadSeqStoreSize readId,
                      TContigStoreSize contigId,
                      TContigSeqSize contigBegin,
                      TReadSeqSize seedBegin,
                      TReadSeqSize seedLength,
                      TReadSeqSize maxErrors,
                      TReadSeqSize minErrors)

{
    if (extender.disabled)
        return false;
        
        
        
        
        
        
    // nhTran's note: MASAI uses approximate seed matching
    //TReadSeqSize errors = seedErrors;
    
    // nhTran's note: AMAS uses exact seed matching
    TReadSeqSize errors = 0;
    
    
    
    
    
    
    TContigSeq & contig = extender.store.contigStore[contigId].seq;
    TReadSeq & read = extender.store.readSeqStore[readId];
    TReadSeqSize readLength = length(read);
    
    
    
    
    

    // Extend left.
    TPatternStateRev patternStateRev;

    TContigSeqSize matchBegin = contigBegin;

    if (seedBegin > 0)
    {
        TContigSeqSize contigLeftBegin = 0;
        if (contigBegin > seedBegin + maxErrors - errors)
            contigLeftBegin = contigBegin - (seedBegin + maxErrors - errors);

        TContigInfix contigLeft(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft(read, 0, seedBegin);

        if (!_extendLeft(extender, patternStateRev, contigLeft, readLeft, maxErrors, errors, matchBegin))
            return false;
    }

    
    
    
    
    
    // This removes some duplicates.
    if (errors < minErrors)
        return false;

    
    
    
    
    
    // Extend right.
    TPatternState patternState;   

    TContigSeqSize matchEnd = contigBegin + seedLength;

    if (seedBegin + seedLength < readLength)
    {
        TContigSeqSize contigRightEnd = extender.genome.contigs_Size[contigId];
        if (contigRightEnd > contigBegin + readLength - seedBegin + maxErrors - errors)
            contigRightEnd = contigBegin + readLength - seedBegin + maxErrors - errors;

        if (contigBegin + seedLength >= contigRightEnd)
            return false;

        TContigInfix contigRight(contig, contigBegin + seedLength, contigRightEnd);
        TReadInfix readRight(read, seedBegin + seedLength, readLength);

        if (!_extendRight(extender, patternState, contigRight, readRight, maxErrors, errors, matchEnd))
            return false;
    }

    
    
    
    
    
    // This removes some duplicates.
    if (errors < minErrors)
        return false;

    
    
    
    
    
    // nhTran's note: record the number of matches for each read
    extender.reads.matchesCount.at(readId)++;
    





	// writing???
	bool reverseComplemented = false;
	if (readId >= extender.reads.readsCount) // reverse complement
    {
		reverseComplemented = true;
		readId = readId - extender.reads.readsCount;
    }
	onMatch(extender.writer /*extender.manager*/, contigId, matchBegin, matchEnd, readId, errors, reverseComplemented);

    return true;
}






template <typename TGenome, typename TReads, typename TWriter /*TManager*/>
inline bool _extendLeft(Extender<TGenome, TReads, TWriter /*TManager*/, EditDistance> const & extender,
                        TPatternStateRev & patternStateRev,
                        TContigInfix & contigInfix,
                        TReadInfix & readInfix,
                        TReadSeqSize maxErrors,
                        TReadSeqSize & errors,
                        TContigSeqSize & matchBegin)
{
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;
    typedef ModifiedString<TContigInfix, ModReverse>        TContigInfixRev;
    typedef Finder<TContigInfixRev>                         TFinder;

    
    
    
    
    
    // Lcp trick.
    TContigSeqSize lcp = 0;
    {  // TODO(holtgrew): Workaround to storing and returning copies in host() for nested infixes/modified strings. This is ugly and should be fixed later.
        TReadInfixRev readInfixRev(readInfix);
        TContigInfixRev contigInfixRev(contigInfix);
        lcp = lcpLength(contigInfixRev, readInfixRev);
    }
    if (lcp == length(readInfix))
    {
        matchBegin -= lcp;
        return true;
    }
    
    
    
    
    
    
    setEndPosition(contigInfix, endPosition(contigInfix) - lcp);
    setEndPosition(readInfix, endPosition(readInfix) - lcp);

    TReadSeqSize remainingErrors = maxErrors - errors;
    TReadSeqSize minErrors = remainingErrors + 1;

    // Stop seed extension.
    if (!remainingErrors)
        return false;

    TContigSeqSize endPos = 0;

    // Align.
    TReadInfixRev readInfixRev(readInfix);
    TContigInfixRev contigInfixRev(contigInfix);
    TFinder finder(contigInfixRev);
    patternStateRev.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readInfixRev, patternStateRev, -static_cast<int>(remainingErrors)))
    {
        TReadSeqSize currentErrors = -getScore(patternStateRev);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = position(finder) + 1;
        }
    }

    errors += minErrors;
    matchBegin -= endPos + lcp;

    return errors <= maxErrors;
}






template <typename TGenome, typename TReads, typename TWriter /*TManager*/>
inline bool _extendRight(Extender<TGenome, TReads, TWriter /*TManager*/, EditDistance> const & extender,
                         TPatternState & patternState,
                         TContigInfix & contigInfix,
                         TReadInfix & readInfix,
                         TReadSeqSize maxErrors,
                         TReadSeqSize & errors,
                         TContigSeqSize & matchEnd)
{
    typedef Finder<TContigInfix>    TFinder;

    // Lcp trick.
    TContigSeqSize lcp = lcpLength(contigInfix, readInfix);
    if (lcp == length(readInfix))
    {
        matchEnd += lcp;
        return true;
    }
    else if (lcp == length(contigInfix))
    {
        errors += length(readInfix) - length(contigInfix);
        matchEnd += length(readInfix);
        return errors <= maxErrors;
    }
    setBeginPosition(contigInfix, beginPosition(contigInfix) + lcp);
    setBeginPosition(readInfix, beginPosition(readInfix) + lcp);

    // NOTE Uncomment this to disable lcp trick.
//    TContigSeqSize lcp = 0;

    
    
    
    
    
    TReadSeqSize remainingErrors = maxErrors - errors;
    TReadSeqSize minErrors = remainingErrors + 1;
    TContigSeqSize endPos = 0;

    // NOTE Comment this to disable lcp trick.
    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Remove last base.
    TContigInfix contigPrefix(contigInfix);
    TReadInfix readPrefix(readInfix);
    setEndPosition(contigPrefix, endPosition(contigPrefix) - 1);
    setEndPosition(readPrefix, endPosition(readPrefix) - 1);

    // Align.
    TFinder finder(contigPrefix);
    patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readPrefix, patternState, -static_cast<int>(remainingErrors)))
    {
        TContigSeqSize currentEnd = position(finder) + 1;
        TReadSeqSize currentErrors = -getScore(patternState);

        // Compare last base.
        if (contigInfix[currentEnd] != back(readInfix))
            if (++currentErrors > remainingErrors)
                continue;

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    errors += minErrors;
    matchEnd += endPos + lcp + 1;

    return errors <= maxErrors;
}






#endif  // #ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
