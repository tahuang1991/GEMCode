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
  deltaRTrackExtra_ = trackExtra.getParameter<double>("deltaR");

  auto recoChargedCandidate = conf().getParameter<edm::ParameterSet>("recoChargedCandidate");
  recoChargedCandidateInputLabel_ = recoChargedCandidate.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  maxBXRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("minBX");
  verboseRecoChargedCandidate_ = recoChargedCandidate.getParameter<int>("verbose");
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
  if (gemvalidation::getByLabel(trackExtraInputLabel_, l2TrackExtras, event())) matchTrackExtraToSimTrack(*l2TrackExtras.product());
  
  // RecoChargedCandidate
  edm::Handle<reco::RecoChargedCandidateCollection> recoChargedCandidates;
  if (gemvalidation::getByLabel(recoChargedCandidateInputLabel_, recoChargedCandidates, event())) matchRecoChargedCandidateToSimTrack(*recoChargedCandidates.product());
}


void 
HLTTrackMatcher::matchTrackExtraToSimTrack(const reco::TrackExtraCollection& tracks)
{
  for(auto& track: tracks) {
    std::cout<<"L2 TrackExtra pT: "<<track.innerMomentum().Rho()
	     <<", eta: "<<track.innerPosition().eta()
     	     <<", phi: "<<track.innerPosition().phi()<<std::endl;  
    std::cout<< "\tDeltaR(SimTrack, L2TrackExtra): " << reco::deltaR(track.innerPosition(), trk().momentum()) << std::endl;
    std::cout<< "\tDeltaPt(SimTrack, L2TrackExtra): " << std::fabs(track.innerMomentum().Rho()-trk().momentum().pt()) << std::endl;

    std::cout << "Rechits/Segments: " << track.recHitsSize()<< std::endl;
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
	std::cout << "\tDT :: id :: " << DTChamberId(id) << std::endl;
	const DTRecSegment4D *seg = dynamic_cast<const DTRecSegment4D*>(*rh);
	std::cout << "\t   :: segment :: " << *seg << std::endl;
	// if RecoTrack RegSegment was found in matching RecSegments
	if (dt_rechit_matcher_->isDTRecSegment4DMatched(*seg)) {
	  ++matchingDTSegments;
	  ++matchingSegments;
	}
      }
      if (is_rpc(id)) {
	std::cout << "\tRPC :: id :: " << RPCDetId(id) << std::endl;
	const RPCRecHit* rpcrh = dynamic_cast<const RPCRecHit*>(*rh);
	std::cout << "\t    :: rechit :: " << *rpcrh << std::endl;
	++matchingRPCSegments;
	++matchingSegments;
      }
      if (is_csc(id)) {
	std::cout << "\tCSC " << CSCDetId(id) << std::endl;
	const CSCSegment *seg = dynamic_cast<const CSCSegment*>(*rh);
	std::cout << "\t    :: segment :: " << *seg << std::endl;
	++matchingCSCSegments;
	++matchingSegments;
      }
    }
    std::cout << "Valid Segments:    " << nValidSegments << std::endl << std::endl;
    std::cout << "Matching Segments: " << matchingSegments << std::endl << std::endl;
    std::cout <<"\tRPC:   " << matchingRPCSegments << std::endl;
    std::cout <<"\tCSC:   " << matchingCSCSegments << std::endl;
    std::cout <<"\tDT:    " << matchingDTSegments << std::endl;
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


