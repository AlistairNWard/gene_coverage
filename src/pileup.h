// ***************************************************************************
// Modified from bamtools_pileup_engine.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// Alistair Ward
// USTAR Center for Genetic Discovery
// University of Utah School of Medicine
// ---------------------------------------------------------------------------
// Last modified: 12 September 2015
// ---------------------------------------------------------------------------
// Provides pileup at position functionality for various tools.
// ***************************************************************************

#ifndef PILEUP_H
#define PILEUP_H

#include "utils/utils_global.h"
#include <api/BamAlignment.h>
#include <vector>

namespace BamTools {

// Contains auxiliary data about a single BamAlignment at current position considered.
struct PileupAlignment {
  
  // Data members.
  BamAlignment Alignment;
  int32_t PositionInAlignment;
  bool IsCurrentDeletion;
  bool IsNextDeletion;
  bool IsNextInsertion;
  int DeletionLength;
  int InsertionLength;
  bool IsSegmentBegin;
  bool IsSegmentEnd;
  
  // Constructor.
  PileupAlignment(const BamAlignment& al)
    : Alignment(al)
    , PositionInAlignment(-1)
    , IsCurrentDeletion(false)
    , IsNextDeletion(false)
    , IsNextInsertion(false)
    , DeletionLength(0)
    , InsertionLength(0)
    , IsSegmentBegin(false)
    , IsSegmentEnd(false)
  {}
};
  
// Contains all data at a position.
struct PileupPosition {
  
  // Data members.
  int RefId;
  int Position;
  std::vector<PileupAlignment> PileupAlignments;

  // Constructor.
  PileupPosition(const int& refId = 0,
    const int& position = 0, 
    const std::vector<PileupAlignment>& alignments = std::vector<PileupAlignment>())
    : RefId(refId)
    , Position(position)
    , PileupAlignments(alignments)
  {}
};

class PileupEngine {
  public:
    PileupEngine(void);
    ~PileupEngine(void);
      
  public:
    bool AddAlignment(const BamAlignment& al, const BamRegion& region, coverageData& cov);
    void Flush(const BamRegion& region, coverageData& cov);

  private:
    struct PileupEnginePrivate;
    PileupEnginePrivate* d;
  };

} // namespace BamTools

#endif // PILEUP_H
