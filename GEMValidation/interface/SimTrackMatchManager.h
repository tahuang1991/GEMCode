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
#include "GEMCode/GEMValidation/interface/RPCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/CSCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/CSCStubMatcher.h"
#include "GEMCode/GEMValidation/interface/GEMRecHitMatcher.h"
//#include "GEMCode/GEMValidation/interface/DTDigiMatcher.h"
//#include "GEMCode/GEMValidation/interface/DTRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/TrackMatcher.h"

class SimTrackMatchManager
{
public:
  
  SimTrackMatchManager(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es);
  
  ~SimTrackMatchManager();

  const SimHitMatcher& simhits() const {return simhits_;}
  const GEMDigiMatcher& gemDigis() const {return gem_digis_;}
  const RPCDigiMatcher& rpcDigis() const {return rpc_digis_;}
  const CSCDigiMatcher& cscDigis() const {return csc_digis_;}
  const CSCStubMatcher& cscStubs() const {return stubs_;}
  const GEMRecHitMatcher& gemRecHits() const {return gem_rechits_;}
  // Add matcher for DT digis & 4D segments
  /* const DTDigiMatcher& dtDigis() const {return dt_digis_;} */
  /* const DTRecHitMatcher& dtRecHits() const {return dt_rechits_;} */
  const TrackMatcher& tracks() const {return tracks_;}
  
private:

  SimHitMatcher simhits_;
  GEMDigiMatcher gem_digis_;
  RPCDigiMatcher rpc_digis_;
  CSCDigiMatcher csc_digis_;
  CSCStubMatcher stubs_;
  GEMRecHitMatcher gem_rechits_;
  /* DTDigiMatcher dt_digis_; */
  /* DTRecHitMatcher dt_rechits_; */
  TrackMatcher tracks_;
};

#endif
