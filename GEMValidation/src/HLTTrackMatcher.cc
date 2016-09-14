#include "GEMCode/GEMValidation/interface/HLTTrackMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TLorentzVector.h"
#include <map>

HLTTrackMatcher::HLTTrackMatcher(CSCRecHitMatcher& csc, 
                                 DTRecHitMatcher& dt, 
                                 RPCRecHitMatcher& rpc, 
                                 GEMRecHitMatcher& gem,  
                                 edm::EDGetTokenT<reco::TrackExtraCollection> &recoTrackExtraInputLabel_,
                                 edm::EDGetTokenT<reco::TrackCollection> &recoTrackInputLabel_,
                                 edm::EDGetTokenT<reco::RecoChargedCandidateCollection> &recoChargedCandidateInputLabel_
                                 )
  : BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
, gem_rechit_matcher_(&gem)
, dt_rechit_matcher_(&dt)
, rpc_rechit_matcher_(&rpc)
, csc_rechit_matcher_(&csc)
{
  auto recoTrackExtra = conf().getParameter<edm::ParameterSet>("recoTrackExtra");
  minBXRecoTrackExtra_ = recoTrackExtra.getParameter<int>("minBX");
  maxBXRecoTrackExtra_ = recoTrackExtra.getParameter<int>("minBX");
  verboseRecoTrackExtra_ = recoTrackExtra.getParameter<int>("verbose");
  runRecoTrackExtra_ = recoTrackExtra.getParameter<bool>("run");

  auto recoTrack = conf().getParameter<edm::ParameterSet>("recoTrack");
  minBXRecoTrack_ = recoTrack.getParameter<int>("minBX");
  maxBXRecoTrack_ = recoTrack.getParameter<int>("minBX");
  verboseRecoTrack_ = recoTrack.getParameter<int>("verbose");
  runRecoTrack_ = recoTrack.getParameter<bool>("run");

  auto recoChargedCandidate = conf().getParameter<edm::ParameterSet>("recoChargedCandidate");
  minBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  maxBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  verboseRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("verbose");
  runRecoChargedCandidate_ = recoChargedCandidate.getParameter<bool>("run");

  // RecoTrackExtra 
  edm::Handle<reco::TrackExtraCollection> recoTrackExtras;
  if (gemvalidation::getByToken(recoTrackExtraInputLabel_, recoTrackExtras, event())) if (runRecoTrackExtra_) matchRecoTrackExtraToSimTrack(*recoTrackExtras.product());
  
  // RecoTrack 
  edm::Handle<reco::TrackCollection> recoTracks;
  if (gemvalidation::getByToken(recoTrackInputLabel_, recoTracks, event())) if (runRecoTrack_) matchRecoTrackToSimTrack(*recoTracks.product());

  // RecoChargedCandidate
  edm::Handle<reco::RecoChargedCandidateCollection> recoChargedCandidates;
  if (gemvalidation::getByToken(recoChargedCandidateInputLabel_, recoChargedCandidates, event())) if (runRecoChargedCandidate_) matchRecoChargedCandidateToSimTrack(*recoChargedCandidates.product());
}


HLTTrackMatcher::~HLTTrackMatcher() 
{
}

void 
HLTTrackMatcher::clear()
{
}

void 
HLTTrackMatcher::matchRecoTrackExtraToSimTrack(const reco::TrackExtraCollection& tracks)
{
  if (verboseRecoTrackExtra_) std::cout << "Number of RecoTrackExtras: " <<tracks.size() << std::endl;
  for(auto& track: tracks) {
    // do not anlyze tracsks with large deltaR
    if (reco::deltaR(track.innerPosition(), trk().momentum()) > 0.5) continue;
    if (verboseRecoTrackExtra_) {
      std::cout << "RecoTrackExtra " << std::endl
		<< "\tpT_inner: "<<track.innerMomentum().Rho()
	        << ", eta_inner: "<<track.innerPosition().eta()
	        << ", phi_inner: "<<track.innerPosition().phi()
		<< "\tpT_outer: "<<track.outerMomentum().Rho()
	        << ", eta_outer: "<<track.outerPosition().eta()
	        << ", phi_outer: "<<track.outerPosition().phi()
		<<std::endl;  
      std::cout << "\tDeltaR(SimTrack, RecoTrackExtra): " << reco::deltaR(track.innerPosition(), trk().momentum()) << std::endl;
      std::cout << "\tDeltaPt(SimTrack, RecoTrackExtra): " << std::fabs(track.innerMomentum().Rho()-trk().momentum().pt()) << std::endl;     
      std::cout << "\tRechits/Segments: " << track.recHitsSize()<< std::endl;
    }
    int matchingCSCSegments(0);
    int matchingRPCSegments(0);
    int matchingGEMSegments(0);
    int matchingDTSegments(0);
    int matchingSegments(0);
    int nValidSegments(0);
    for(auto rh = track.recHitsBegin(); rh != track.recHitsEnd(); rh++) {
      if (!(**rh).isValid()) continue;
      ++nValidSegments;
      auto id((**rh).rawId());
      if (gemvalidation::is_dt(id)) {
	const DTRecSegment4D *seg = dynamic_cast<const DTRecSegment4D*>(*rh);
	if (verboseRecoTrackExtra_) {
	  std::cout << "\t\tDT  :: id :: " << DTChamberId(id) << std::endl;
	  std::cout << "\t\t    :: segment :: " << *seg << std::endl;
	}
	if (dt_rechit_matcher_->isDTRecSegment4DMatched(*seg)) {
	  if (verboseRecoTrackExtra_) std::cout << "\t\t    :: MATCHED!" << std::endl;
	  ++matchingDTSegments;
	  ++matchingSegments;
	}
      }
      if (gemvalidation::is_rpc(id)) {
	const RPCRecHit* rpcrh = dynamic_cast<const RPCRecHit*>(*rh);
	if (verboseRecoTrackExtra_) {
	  std::cout << "\t\tRPC :: id :: " << RPCDetId(id) << std::endl;
	  std::cout << "\t\t    :: rechit :: " << *rpcrh << std::endl;
	}
	if (rpc_rechit_matcher_->isRPCRecHitMatched(*rpcrh)) {
	  if (verboseRecoTrackExtra_) std::cout << "\t\t    :: MATCHED!" << std::endl;
	  ++matchingRPCSegments;
	  ++matchingSegments;
	}
      }
      if (gemvalidation::is_gem(id)) {
	const GEMRecHit* gemrh = dynamic_cast<const GEMRecHit*>(*rh);
	if (verboseRecoTrackExtra_) {
	  std::cout << "\t\tGEM :: id :: " << GEMDetId(id) << std::endl;
	  std::cout << "\t\t    :: rechit :: " << *gemrh << std::endl;
	}
	if (gem_rechit_matcher_->isGEMRecHitMatched(*gemrh)) {
	  if (verboseRecoTrackExtra_) std::cout << "\t\t    :: MATCHED!" << std::endl;
	  ++matchingGEMSegments;
	  ++matchingSegments;
	}
      }
      if (gemvalidation::is_csc(id)) {
	const CSCSegment *seg = dynamic_cast<const CSCSegment*>(*rh);
	if (verboseRecoTrackExtra_) {
	  std::cout << "\t\tCSC :: id :: " << CSCDetId(id) << std::endl;
	  std::cout << "\t\t    :: segment :: " << *seg << std::endl;
	}
	if (csc_rechit_matcher_->isCSCSegmentMatched(*seg)) {
	  if (verboseRecoTrackExtra_) std::cout << "\t\t    :: MATCHED!" << std::endl;
	  ++matchingCSCSegments;
	  ++matchingSegments;
	}
      }
    }
    if (verboseRecoTrackExtra_) {
      std::cout << "\tValid Segments:    " << nValidSegments << std::endl;
      std::cout << "\tMatching Segments: " << matchingSegments << std::endl;
      std::cout << "\t              RPC: " << matchingRPCSegments << std::endl;
      std::cout << "\t              CSC: " << matchingCSCSegments << std::endl;
      std::cout << "\t              GEM: " << matchingGEMSegments << std::endl;
      std::cout << "\t               DT: " << matchingDTSegments << std::endl;
    }
    // store matching L1RecoTrackExtra
    if (matchingDTSegments + matchingCSCSegments>=2) {
      if (verboseRecoTrackExtra_) {
	std::cout << "\tRecoTrackExtra was matched! (deltaR = " << reco::deltaR(track.innerPosition(), trk().momentum()) << ") " << std::endl;
      }
      matchedRecoTrackExtras_.push_back(track);
    }
  }
}


void 
HLTTrackMatcher::matchRecoTrackToSimTrack(const reco::TrackCollection& tracks)
{
  if (verboseRecoTrack_) std::cout << "Number of RecoTracks: " <<tracks.size() << std::endl;
  int i=0;
  for(auto& track: tracks) {
    const double deltaR(reco::deltaR(track.outerEta(), track.outerPhi(), trk().momentum().eta(), trk().momentum().phi()));
    if (verboseRecoTrack_) {
      std::cout<< "RecoTrack " << i+1 << " - pT: "<<track.outerPt()
	       <<", eta: "<<track.outerEta()
	       <<", phi: "<<track.outerPhi()
	       <<", deltaR: "<< deltaR << std::endl;
    }
    // check if the associated RecoTrackExtra was matched!
    for (auto &otherTrackExtra: getMatchedRecoTrackExtras()) {
      if (areRecoTrackSame(*(track.extra()), otherTrackExtra)) {
     	if (verboseRecoTrack_) std::cout << "\tRecoTrack was matched!" << std::endl;
     	matchedRecoTracks_.push_back(track);
      }
    }
    ++i;
  }
}


void 
HLTTrackMatcher::matchRecoChargedCandidateToSimTrack(const reco::RecoChargedCandidateCollection& candidates)
{
  if (verboseRecoTrack_) std::cout << "Number of RecoChargedCandidates: " <<candidates.size() << std::endl;
  int i=0;
  for(auto& candidate: candidates) {
    const double deltaR(reco::deltaR(candidate.eta(), candidate.phi(), trk().momentum().eta(), trk().momentum().phi()));
    if (verboseRecoChargedCandidate_) {
      std::cout<< "RecoCandidate " << i+1 << " - pT: "<<candidate.pt()
	       <<", eta: "<<candidate.eta()
	       <<", phi: "<<candidate.phi()
	       <<", deltaR: "<< deltaR << std::endl;
    }
    // get the RecoTrack
    for (auto &otherTrack: getMatchedRecoTracks()) {
      if (areRecoTrackSame(*(candidate.track()), otherTrack)) {
     	if (verboseRecoTrack_) std::cout << "\tRecoChargedCandidate was matched!" << std::endl;
	matchedRecoChargedCandidates_.push_back(candidate);
      }
    }
    ++i;    
  }
}


template<typename T>
bool
HLTTrackMatcher::areRecoTrackSame(const T& l, const T& r) const
{
  return (l.outerPosition()==r.outerPosition() and 
	  l.outerMomentum()==r.outerMomentum() and
	  l.outerDetId()==r.outerDetId() and 
	  l.outerDetId()==r.outerDetId() and
	  l.innerPosition()==r.innerPosition() and 
	  l.innerMomentum()==r.innerMomentum() and
	  l.innerOk()==r.innerOk() and
	  l.outerDetId()==r.outerDetId());
}

