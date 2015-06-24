#ifndef GEMCode_GEMValidation_HLTTrackMatcher_h
#define GEMCode_GEMValidation_HLTTrackMatcher_h

/**\class HLTTrackMatcher

 Description: Matching of tracks to SimTrack

 Original Author:  "Sven Dildick"
*/

#include "GEMCode/GEMValidation/interface/CSCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/DTRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/RPCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/GEMRecHitMatcher.h"
//#include "GEMCode/GEMValidation/interface/ME0RecHitMatcher.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

class HLTTrackMatcher : public BaseMatcher
{
 public:
  /// constructor
  HLTTrackMatcher(CSCRecHitMatcher&, DTRecHitMatcher&, 
		  RPCRecHitMatcher&, GEMRecHitMatcher&);
  /// destructor
  ~HLTTrackMatcher();

 private:

  void init();
  void clear();

  void matchTrackExtraToSimTrack(const reco::TrackExtraCollection&);
  void matchRecoChargedCandidateToSimTrack(const reco::RecoChargedCandidateCollection&);

  const GEMRecHitMatcher* gem_rechit_matcher_;
  const DTRecHitMatcher* dt_rechit_matcher_;
  const RPCRecHitMatcher* rpc_rechit_matcher_;
  const CSCRecHitMatcher* csc_rechit_matcher_;

  std::vector<edm::InputTag> trackExtraInputLabel_;
  std::vector<edm::InputTag> recoChargedCandidateInputLabel_;

  int minBXTrackExtra_, maxBXTrackExtra_;
  int minBXRecoChargedCandidate_, maxBXRecoChargedCandidate_;

  int verboseTrackExtra_;
  int verboseRecoChargedCandidate_;

  double deltaRTrackExtra_;
  double deltaRRecoChargedCandidate_;
};

#endif
