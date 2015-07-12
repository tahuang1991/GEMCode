#ifndef GEMCode_GEMValidation_L1TrackFinderCandidateMatcher_h
#define GEMCode_GEMValidation_L1TrackFinderCandidateMatcher_h

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

class L1CSCTrackCollection;

class L1TrackFinderCandidateMatcher //: public BaseMatcher
{
 public:
  /// constructor
  L1TrackFinderCandidateMatcher();
  /// destructor
  ~L1TrackFinderCandidateMatcher();
  
 private:
  
  void clear();
  void init(); 
  
  void matchCSCTfCandToSimTrack(const L1CSCTrackCollection&); 
  void matchDTTfCandToSimTrack(const L1CSCTrackCollection&); 

  std::vector<edm::InputTag> cscTfCandInputLabel_; 
  std::vector<edm::InputTag> dtTfCandInputLabel_; 
};

#endif
