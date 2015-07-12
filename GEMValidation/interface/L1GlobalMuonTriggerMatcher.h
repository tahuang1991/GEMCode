#ifndef GEMCode_GEMValidation_L1GlobalMuonTriggerMatcher_h
#define GEMCode_GEMValidation_L1GlobalMuonTriggerMatcher_h

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

class L1CSCTrackCollection;

class L1GlobalMuonTriggerMatcher //: public BaseMatcher
{
 public:
  /// constructor
  L1GlobalMuonTriggerMatcher();
  /// destructor
  ~L1GlobalMuonTriggerMatcher();
  
 private:
  
  void clear();
  void init(); 
  
  void matchRegionalCandCSCToSimTrack(const L1CSCTrackCollection&); 
  void matchRegionalCandDTToSimTrack(const L1CSCTrackCollection&); 
  void matchRegionalCandRPCbToSimTrack(const L1CSCTrackCollection&); 
  void matchRegionalCandRPCfToSimTrack(const L1CSCTrackCollection&); 
  void matchGMTCandToSimTrack(const L1CSCTrackCollection&); 
};

#endif
