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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class HLTTrackMatcher : public BaseMatcher
{
 public:
  /// constructor
  HLTTrackMatcher(CSCRecHitMatcher&, DTRecHitMatcher&, 
                  RPCRecHitMatcher&, GEMRecHitMatcher&,
                  edm::EDGetTokenT<reco::TrackExtraCollection>&,
                  edm::EDGetTokenT<reco::TrackCollection>&,
                  edm::EDGetTokenT<reco::RecoChargedCandidateCollection>&);
  /// destructor
  ~HLTTrackMatcher();

  const reco::TrackExtraCollection& getMatchedRecoTrackExtras() const {return matchedRecoTrackExtras_;}
  const reco::TrackCollection& getMatchedRecoTracks() const {return matchedRecoTracks_;}
  const reco::RecoChargedCandidateCollection& getMatchedRecoChargedCandidates() const {return matchedRecoChargedCandidates_;}

 private:

  void init();
  void clear();

  void matchRecoTrackExtraToSimTrack(const reco::TrackExtraCollection&);
  void matchRecoTrackToSimTrack(const reco::TrackCollection&);
  void matchRecoChargedCandidateToSimTrack(const reco::RecoChargedCandidateCollection&);

  template<typename T>
  bool areRecoTrackSame(const T&, const T&) const;

  const GEMRecHitMatcher* gem_rechit_matcher_;
  const DTRecHitMatcher* dt_rechit_matcher_;
  const RPCRecHitMatcher* rpc_rechit_matcher_;
  const CSCRecHitMatcher* csc_rechit_matcher_;

  int minBXRecoTrackExtra_, maxBXRecoTrackExtra_;
  int verboseRecoTrackExtra_;
  bool runRecoTrackExtra_;

  int minBXRecoTrack_, maxBXRecoTrack_;
  int verboseRecoTrack_;
  bool runRecoTrack_;

  int minBXRecoChargedCandidate_, maxBXRecoChargedCandidate_;
  int verboseRecoChargedCandidate_;
  bool runRecoChargedCandidate_;

  reco::TrackExtraCollection matchedRecoTrackExtras_;
  reco::TrackCollection matchedRecoTracks_;
  reco::RecoChargedCandidateCollection matchedRecoChargedCandidates_;
};

#endif
