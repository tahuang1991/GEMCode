
#ifndef GEMCode_GEMValidation_SimTrackMatchManager_h
#define GEMCode_GEMValidation_SimTrackMatchManager_h

/**\class SimTrackMatchManager

 Description: Matching of SIM and Trigger info for a SimTrack in Muon subsystems

 It's a manager-matcher class, as it uses specialized matching classes to match SimHits, various digis and stubs.

 Original Author:  "Vadim Khotilovich"
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"
#include "GEMCode/GEMValidation/interface/GEMDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/GEMRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/ME0DigiMatcher.h"
//#include "GEMCode/GEMValidation/interface/ME0RecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/RPCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/RPCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/CSCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/CSCStubMatcher.h"
#include "GEMCode/GEMValidation/interface/CSCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/DTDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/DTStubMatcher.h"
#include "GEMCode/GEMValidation/interface/DTRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/L1TrackMatcher.h"
#include "GEMCode/GEMValidation/interface/HLTTrackMatcher.h"

class SimTrackMatchManager
{
public:
  
  SimTrackMatchManager(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es);
  
  ~SimTrackMatchManager();

  const SimHitMatcher& simhits() const {return simhits_;}
  const GEMDigiMatcher& gemDigis() const {return gem_digis_;}
  const GEMRecHitMatcher& gemRecHits() const {return gem_rechits_;}
  const ME0DigiMatcher& me0Digis() const {return me0_digis_;}
  //  const ME0RecHitMatcher& me0RecHits() const {return me0_rechits_;}
  const RPCDigiMatcher& rpcDigis() const {return rpc_digis_;}
  const RPCRecHitMatcher& rpcRecHits() const {return rpc_rechits_;}
  const CSCDigiMatcher& cscDigis() const {return csc_digis_;}
  const CSCStubMatcher& cscStubs() const {return csc_stubs_;}
  const CSCRecHitMatcher& cscRecHits() const {return csc_rechits_;}
  const DTDigiMatcher& dtDigis() const {return dt_digis_;}
  const DTStubMatcher& dtStubs() const {return dt_stubs_;}
  const DTRecHitMatcher& dtRecHits() const {return dt_rechits_;} 
  const L1TrackMatcher& l1Tracks() const {return l1_tracks_;}
  const HLTTrackMatcher& hltTracks() const {return hlt_tracks_;}
  
private:

  SimHitMatcher simhits_;
  GEMDigiMatcher gem_digis_;
  GEMRecHitMatcher gem_rechits_;
  ME0DigiMatcher me0_digis_;
  //  ME0RecHitMatcher me0_rechits_;
  RPCDigiMatcher rpc_digis_;
  RPCRecHitMatcher rpc_rechits_;
  CSCDigiMatcher csc_digis_;
  CSCStubMatcher csc_stubs_;
  CSCRecHitMatcher csc_rechits_;
  DTDigiMatcher dt_digis_;
  DTStubMatcher dt_stubs_;
  DTRecHitMatcher dt_rechits_; 
  L1TrackMatcher l1_tracks_;
  HLTTrackMatcher hlt_tracks_;
};

#endif
