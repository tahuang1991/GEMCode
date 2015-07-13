#ifndef GEMCode_GEMValidation_L1TrackFinderCandidateMatcher_h
#define GEMCode_GEMValidation_L1TrackFinderCandidateMatcher_h

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

typedef std::vector<L1MuRegionalCand> L1MuRegionalCandCollection;

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
  
  void matchCSCTfCandToSimTrack(const L1MuRegionalCandCollection&); 
  void matchDTTfCandToSimTrack(const L1MuRegionalCandCollection&); 

  std::vector<edm::InputTag> cscTfCandInputLabel_; 
  std::vector<edm::InputTag> dtTfCandInputLabel_; 
};

#endif
