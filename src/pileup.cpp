// ***************************************************************************
// Modified bamtools_pileup_engine.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// Alistair Ward
// USTAR Center for Genetic Discovery
// University of Utah School of Medicine
// ---------------------------------------------------------------------------
// Last modified: 12 September 2015
// ---------------------------------------------------------------------------
// Provides pileup at position functionality for various tools.
// ***************************************************************************

#include "dataProcessing.h"
#include "mean.h"
#include "pileup.h"
using namespace BamTools;

#include <iostream>
using namespace std;

// PileupEnginePrivate implementation
struct PileupEngine::PileupEnginePrivate {
  
  // Data members.
  int CurrentId;
  int CurrentPosition;
  vector<BamAlignment> CurrentAlignments;
  PileupPosition CurrentPileupData;
  
  bool IsFirstAlignment;

  // Store the pileup data.
  int firstPosition;
  vector<int> pileupCounts;

  // Constructor & destuctor
  PileupEnginePrivate(void)
    : CurrentId(-1)
    , CurrentPosition(-1)
    , IsFirstAlignment(true)
  { }
  ~PileupEnginePrivate(void) { }
  
  // Public methods.
  bool AddAlignment(const BamAlignment& al, const BamRegion& region, coverageData& cov);
  void Flush(const BamRegion& region, coverageData& cov);
  
  // Internal methods.
  private:
    void ClearOldData(void);
    void CreatePileupData(void);
    void ParseAlignmentCigar(const BamAlignment& al);
};

bool PileupEngine::PileupEnginePrivate::AddAlignment(const BamAlignment& al, const BamRegion& region, coverageData& cov) {
  
  // If this is the first alignment.
  if ( IsFirstAlignment ) {
    
    // Set initial markers.
    CurrentId       = al.RefID;
    CurrentPosition = al.Position;
    firstPosition   = al.Position;
    
    // Store the first entry.
    CurrentAlignments.clear();
    CurrentAlignments.push_back(al);
    
    // Set the flag & return.
    IsFirstAlignment = false;
    return true;
  }
  
  // If the alignment has the same reference.
  if ( al.RefID == CurrentId ) {
    
    // If the alignment has the same position, store and move on.
    if ( al.Position == CurrentPosition ) {
      CurrentAlignments.push_back(al);
    
    // If the alignment position is less than CurrentPosition, there is a sorting error => ABORT.
    } else if ( al.Position < CurrentPosition ) {
      cerr << "ERROR: BAM file is not in coordinate sort order as required." << endl;
      return false;
    
    // Otherwise, store the pileup data until 'catching up' to CurrentPosition.
    } else {
      while ( al.Position > CurrentPosition ) {

        // Remove stored alignments that end before the current position and then parse the CIGAR
        // data in the remaining alignments to build up the pileup data.
        CreatePileupData();
        if ( CurrentPosition + 1 >= region.LeftPosition && CurrentPosition + 1< region.RightPosition + 1) {
          cov.update(CurrentPileupData.PileupAlignments.size());
        }
        ++CurrentPosition;
      }
      CurrentAlignments.push_back(al);
    }

  // If the reference ID less than CurrentId, there is a sorting error => ABORT.
  } else if ( al.RefID < CurrentId ) {
      cerr << "ERROR: BAM file is not in coordinate sort order as required." << endl;
      return false;

  // Otherwise, move forward onto the next reference.
  } else {
      
    // Store any remaining pileup data from previous reference.
    while ( !CurrentAlignments.empty() ) {

      // Remove stored alignments that end before the current position and then parse the CIGAR
      // data in the remaining alignments to build up the pileup data.
      CreatePileupData();
      if ( CurrentPosition + 1 >= region.LeftPosition && CurrentPosition + 1< region.RightPosition + 1) {
        cov.update(CurrentPileupData.PileupAlignments.size());
      }
      ++CurrentPosition;
    }
    
    // store first entry on this new reference, update markers
    CurrentAlignments.clear();
    CurrentAlignments.push_back(al);
    CurrentId       = al.RefID;
    CurrentPosition = al.Position;
  }

  return true;
}

// Build up pileup data.
void PileupEngine::PileupEnginePrivate::CreatePileupData(void) {
  
  // Remove any non-overlapping alignments.
  ClearOldData();

  // Set the pileup refId and position to the current markers.
  CurrentPileupData.RefId    = CurrentId;
  CurrentPileupData.Position = CurrentPosition;
  CurrentPileupData.PileupAlignments.clear();
  
  // Parse the CIGAR data in the remaining alignments.
  vector<BamAlignment>::const_iterator alIter = CurrentAlignments.begin();
  vector<BamAlignment>::const_iterator alEnd  = CurrentAlignments.end(); 
  for ( ; alIter != alEnd; ++alIter ) { ParseAlignmentCigar( (*alIter) ); }
}

// Remove any alignments that end before our CurrentPosition.
// N.B. - BAM positions are 0-based, half-open. GetEndPosition() returns a 1-based position,
//        while our CurrentPosition is 0-based. For example, an alignment with 'endPosition' of
//        100 does not overlap a 'CurrentPosition' of 100, and should be discarded.
void PileupEngine::PileupEnginePrivate::ClearOldData(void) {
  size_t i = 0;
  size_t j = 0;
  const size_t numAlignments = CurrentAlignments.size();
  while ( i < numAlignments ) {

    // Skip over alignment if its (1-based) endPosition is <= to (0-based) CurrentPosition
    // i.e. this entry will not be saved upon vector resize (at the end of the routine).
    const int endPosition = CurrentAlignments[i].GetEndPosition();
    if ( endPosition <= CurrentPosition ) {
      ++i;
      continue;
    }

    // Otherwise, the alignment ends after CurrentPosition and so move it towards vector beginning,
    // at index j.
    if ( i != j ) {
      CurrentAlignments[j] = CurrentAlignments[i];
    }

    // Increment the indices.
    ++i;
    ++j;
  }

  // Squeeze the vector to size j, discarding all remaining alignments in the container.
  CurrentAlignments.resize(j);
}

void PileupEngine::PileupEnginePrivate::ParseAlignmentCigar(const BamAlignment& al) {
  
  // Skip if the read is unmapped.
  if ( !al.IsMapped() ) return;
  
  // Intialize local variables.
  int  genomePosition      = al.Position;
  int  positionInAlignment = 0;
  bool isNewReadSegment    = true;
  bool saveAlignment       = true;    
  PileupAlignment pileupAlignment(al);
  
  // Iterate over the CIGAR operations.
  const int numCigarOps = (const int)al.CigarData.size();
  for (int i = 0; i < numCigarOps; ++i ) { 
    const CigarOp& op = al.CigarData.at(i);
  
    // If the op is a MATCH.
    if ( op.Type == 'M' ) {
    
      // If the match op overlaps current position.
      if ( genomePosition + (int)op.Length > CurrentPosition ) {
        
        // Set pileup data.
        pileupAlignment.IsCurrentDeletion   = false;
        pileupAlignment.IsNextDeletion      = false;
        pileupAlignment.IsNextInsertion     = false;
        pileupAlignment.PositionInAlignment = positionInAlignment + (CurrentPosition - genomePosition);
        
        // Check for the beginning of the read segment.
        if ( genomePosition == CurrentPosition && isNewReadSegment ) {
          pileupAlignment.IsSegmentBegin = true;
        }
        
        // If we're at the end of a match op.
        if ( genomePosition + (int)op.Length - 1 == CurrentPosition ) {
            
          // If this is not the last operation.
          if ( i < numCigarOps - 1 ) {
              
            // Check the next CIGAR op.
            const CigarOp& nextOp = al.CigarData.at(i+1);
            
            // If the next CIGAR op is a DELETION.
            if ( nextOp.Type == 'D') {
              pileupAlignment.IsNextDeletion = true;
              pileupAlignment.DeletionLength = nextOp.Length;
            
            // If the next CIGAR op is an INSERTION.
            } else if ( nextOp.Type == 'I' ) {
              pileupAlignment.IsNextInsertion = true;
              pileupAlignment.InsertionLength = nextOp.Length;
            }
                
            // If the next CIGAR op is either a DELETION or an INSERTION.
            if ( nextOp.Type == 'D' || nextOp.Type == 'I' ) {

              // If there is a CIGAR op after the DEL/INS.
              if ( i < numCigarOps - 2 ) {
                const CigarOp& nextNextOp = al.CigarData.at(i+2);
                
                // If the next CIGAR op is a clipping or a ref_skip.
                if ( nextNextOp.Type == 'S' || nextNextOp.Type == 'N' || nextNextOp.Type == 'H' ) {
                  pileupAlignment.IsSegmentEnd = true;
                }
              } else {
                pileupAlignment.IsSegmentEnd = true;
                
                // if next CIGAR op is clipping or ref_skip
                if ( nextOp.Type == 'S' || nextOp.Type == 'N' || nextOp.Type == 'H' ) {
                  pileupAlignment.IsSegmentEnd = true;
                }
              }
          
            // otherwise
            } else { 
          
              // if next CIGAR op is clipping or ref_skip
              if ( nextOp.Type == 'S' || nextOp.Type == 'N' || nextOp.Type == 'H' ) {
                pileupAlignment.IsSegmentEnd = true;
              }
            }
          
          // else this is last operation
          } else {
            pileupAlignment.IsSegmentEnd = true;
          }
        }
      }
      
      // increment markers
      genomePosition      += op.Length;
      positionInAlignment += op.Length;
      
    // if op is DELETION
    } else if ( op.Type == 'D' ) {
    
      // if deletion op overlaps current position
      if ( genomePosition + (int)op.Length > CurrentPosition ) {
        
        // set pileup data
        pileupAlignment.IsCurrentDeletion   = true;
        pileupAlignment.IsNextDeletion      = false;
        pileupAlignment.IsNextInsertion     = true;
        pileupAlignment.PositionInAlignment = positionInAlignment + (CurrentPosition - genomePosition);
      }
      
      // increment marker
      genomePosition += op.Length;

    // if op is REF_SKIP
    } else if ( op.Type == 'N' ) {
      genomePosition += op.Length;
    
    // if op is INSERTION or SOFT_CLIP
    } else if ( op.Type == 'I' || op.Type == 'S' ) {
      positionInAlignment += op.Length;
    }
    
    // Check for beginning of new read segment.
    if ( op.Type == 'N' || op.Type == 'S' || op.Type == 'H' ) {
      isNewReadSegment = true;
    } else {
      isNewReadSegment = false;
    }
  
    // if we've moved beyond current position
    if ( genomePosition > CurrentPosition ) {
      if ( op.Type == 'N' ) saveAlignment = false; // ignore alignment if REF_SKIP
      break;
    }
  }

  // Save the pileup position if the flag is true.
  if ( saveAlignment ) {
    CurrentPileupData.PileupAlignments.push_back( pileupAlignment );
  }
}

// Clear out remaining pileup information.
void PileupEngine::PileupEnginePrivate::Flush(const BamRegion& region, coverageData& cov) {
  while ( !CurrentAlignments.empty() ) {
    CreatePileupData();
    if ( CurrentPosition + 1 >= region.LeftPosition && CurrentPosition + 1< region.RightPosition + 1) {
      cov.update(CurrentPileupData.PileupAlignments.size());
    }
    ++CurrentPosition;
  }
}

// ---------------------------------------------
// PileupEngine implementation

PileupEngine::PileupEngine(void)
    : d( new PileupEnginePrivate )
{}

PileupEngine::~PileupEngine(void) {
    delete d;
    d = 0;
}

bool PileupEngine::AddAlignment(const BamAlignment& al, const BamRegion& region, coverageData& cov) { return d->AddAlignment(al, region, cov); }
void PileupEngine::Flush(const BamRegion& region, coverageData& cov) { d->Flush(region, cov); }
