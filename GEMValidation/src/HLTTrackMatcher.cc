#include "GEMCode/GEMValidation/interface/HLTTrackMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TLorentzVector.h"
#include <map>

HLTTrackMatcher::HLTTrackMatcher(CSCRecHitMatcher& csc, DTRecHitMatcher& dt, 
				 RPCRecHitMatcher& rpc, GEMRecHitMatcher& gem)
: BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
, gem_rechit_matcher_(&gem)
, dt_rechit_matcher_(&dt)
, rpc_rechit_matcher_(&rpc)
, csc_rechit_matcher_(&csc)
{
  auto trackExtra = conf().getParameter<edm::ParameterSet>("trackExtra");
  trackExtraInputLabel_ = trackExtra.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBXTrackExtra_ = trackExtra.getParameter<int>("minBX");
  maxBXTrackExtra_ = trackExtra.getParameter<int>("minBX");
  verboseTrackExtra_ = trackExtra.getParameter<int>("verbose");
  runTrackExtra_ = trackExtra.getParameter<bool>("run");
  deltaRTrackExtra_ = trackExtra.getParameter<double>("deltaR");

  auto recoChargedCandidate = conf().getParameter<edm::ParameterSet>("recoChargedCandidate");
  recoChargedCandidateInputLabel_ = recoChargedCandidate.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  maxBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  verboseRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("verbose");
  runRecoChargedCandidate_ = recoChargedCandidate.getParameter<bool>("run");
  deltaRRecoChargedCandidate_ = recoChargedCandidate.getParameter<double>("deltaR");

  clear();
  init();
}

HLTTrackMatcher::~HLTTrackMatcher() 
{
}

void 
HLTTrackMatcher::clear()
{
}

void 
HLTTrackMatcher::init()
{  
  // TrackExtra 
  edm::Handle<reco::TrackExtraCollection> l2TrackExtras;
  if (gemvalidation::getByLabel(trackExtraInputLabel_, l2TrackExtras, event())) if (runTrackExtra_) matchTrackExtraToSimTrack(*l2TrackExtras.product());
  
  // RecoChargedCandidate
  edm::Handle<reco::RecoChargedCandidateCollection> recoChargedCandidates;
  if (gemvalidation::getByLabel(recoChargedCandidateInputLabel_, recoChargedCandidates, event())) if (runRecoChargedCandidate_) matchRecoChargedCandidateToSimTrack(*recoChargedCandidates.product());
}


void 
HLTTrackMatcher::matchTrackExtraToSimTrack(const reco::TrackExtraCollection& tracks)
{
  if (verboseTrackExtra_) std::cout << "Number of L1ExtraTracks: " <<tracks.size() << std::endl;
  for(auto& track: tracks) {
    // do not anlyze tracsks with large deltaR
    if (reco::deltaR(track.innerPosition(), trk().momentum()) > 0.5) continue;
    if (verboseTrackExtra_) {
      std::cout<<"L2 TrackExtra pT: "<<track.innerMomentum().Rho()
	       <<", eta: "<<track.innerPosition().eta()
	       <<", phi: "<<track.innerPosition().phi()<<std::endl;  
      std::cout<< "\tDeltaR(SimTrack, L2TrackExtra): " << reco::deltaR(track.innerPosition(), trk().momentum()) << std::endl;
      std::cout<< "\tDeltaPt(SimTrack, L2TrackExtra): " << std::fabs(track.innerMomentum().Rho()-trk().momentum().pt()) << std::endl;     
      std::cout << "\tRechits/Segments: " << track.recHitsSize()<< std::endl;
    }
    int matchingCSCSegments(0);
    int matchingRPCSegments(0);
    int matchingDTSegments(0);
    int matchingSegments(0);
    int nValidSegments(0);
    for(auto rh = track.recHitsBegin(); rh != track.recHitsEnd(); rh++) {
      if (!(**rh).isValid()) continue;
      ++nValidSegments;
      auto id((**rh).rawId());
      if (is_dt(id)) {
	const DTRecSegment4D *seg = dynamic_cast<const DTRecSegment4D*>(*rh);
	if (verboseTrackExtra_) {
	  std::cout << "\t\tDT :: id :: " << DTChamberId(id) << std::endl;
	  std::cout << "\t\t   :: segment :: " << *seg << std::endl;
	}
	if (dt_rechit_matcher_->isDTRecSegment4DMatched(*seg)) {
	  ++matchingDTSegments;
	  ++matchingSegments;
	}
      }
      if (is_rpc(id)) {
	const RPCRecHit* rpcrh = dynamic_cast<const RPCRecHit*>(*rh);
	if (verboseTrackExtra_) {
	  std::cout << "\t\tRPC :: id :: " << RPCDetId(id) << std::endl;
	  std::cout << "\t\t    :: rechit :: " << *rpcrh << std::endl;
	}
	if (rpc_rechit_matcher_->isRPCRecHitMatched(*rpcrh)) {
	  ++matchingRPCSegments;
	  ++matchingSegments;
	}
      }
      if (is_csc(id)) {
	const CSCSegment *seg = dynamic_cast<const CSCSegment*>(*rh);
	if (verboseTrackExtra_) {
	  std::cout << "\t\tCSC " << CSCDetId(id) << std::endl;
	  std::cout << "\t\t    :: segment :: " << *seg << std::endl;
	}
	if (csc_rechit_matcher_->isCSCSegmentMatched(*seg)) {
	  ++matchingCSCSegments;
	  ++matchingSegments;
	}
      }
    }
    if (verboseTrackExtra_) {
      std::cout << "\tValid Segments:    " << nValidSegments << std::endl;
      std::cout << "\tMatching Segments: " << matchingSegments << std::endl;
      std::cout << "\t              RPC: " << matchingRPCSegments << std::endl;
      std::cout << "\t              CSC: " << matchingCSCSegments << std::endl;
      std::cout << "\t               DT: " << matchingDTSegments << std::endl;
    }
    // store matching L1TrackExtra
    if (matchingCSCSegments>=2 or matchingDTSegments>=2) {
      if (verboseTrackExtra_) std::cout << "\tTrackExtra was matched!" << std::endl;
      matchedTrackExtras_.push_back(track);
    }
  }
}


void 
HLTTrackMatcher::matchRecoChargedCandidateToSimTrack(const reco::RecoChargedCandidateCollection& candidates)
{
  // for(auto& muon: candidates) {
  //   std::cout<<" L1 Muon Particle PT: "<<muon.pt()
  //            <<", eta: "<<muon.eta()
  //            <<", charge: "<<muon.charge()
  //            <<", phi: "<<muon.phi()<<std::endl;  
  // }
}


