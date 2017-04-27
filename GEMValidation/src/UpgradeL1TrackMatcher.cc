#include "GEMCode/GEMValidation/interface/UpgradeL1TrackMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TLorentzVector.h"
#include <map>

UpgradeL1TrackMatcher::UpgradeL1TrackMatcher(CSCStubMatcher& csc,
                                             edm::EDGetTokenT<l1t::EMTFTrackCollection> &emtfTrackInputLabel_,
                                             edm::EDGetTokenT< BXVector<l1t::RegionalMuonCand> > & gmtInputLabel_)
  : BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
  , csc_stub_matcher_(&csc)
{
  auto tfTrack = conf().getParameter<edm::ParameterSet>("upgradeEmtfTrack");
  minBXEMTFTrack_ = tfTrack.getParameter<int>("minBX");
  maxBXEMTFTrack_ = tfTrack.getParameter<int>("maxBX");
  verboseEMTFTrack_ = tfTrack.getParameter<int>("verbose");
  deltaREMTFTrack_ = tfTrack.getParameter<double>("deltaR");

  auto gmt = conf().getParameter<edm::ParameterSet>("upgradeGMT");
  minBXGMT_ = gmt.getParameter<int>("minBX");
  maxBXGMT_ = gmt.getParameter<int>("maxBX");
  verboseGMT_ = gmt.getParameter<int>("verbose");
  deltaRGMT_ = gmt.getParameter<double>("deltaR");

  //std::cout<<" UpgradeL1TrackMatcher constructor" <<std::endl;
  clear();

  simPt = trk().momentum().pt();
  simEta = trk().momentum().eta();
  simPhi = trk().momentum().phi();
  simE = trk().momentum().E();
  simCharge = trk().charge();

  // tracks produced by EMEMTF
  edm::Handle<l1t::EMTFTrackCollection> hl1Tracks;
  if (gemvalidation::getByToken(emtfTrackInputLabel_,hl1Tracks, event())) matchEmtfTrackToSimTrack(*hl1Tracks.product());
  else 
      std::cout  <<"failed readout EMTFTracks " << std::endl;

  edm::Handle<BXVector<l1t::RegionalMuonCand>> hGMT;
  if (gemvalidation::getByToken(gmtInputLabel_,hGMT, event())) matchGMTToSimTrack(*hGMT.product());
  else 
      std::cout  <<"failed readout RegionalMuonCand " << std::endl;
}

UpgradeL1TrackMatcher::~UpgradeL1TrackMatcher()
{
}

void
UpgradeL1TrackMatcher::clear()
{
  tfTracks_.clear();
}

void
UpgradeL1TrackMatcher::matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection& tracks)
{
  GlobalPoint gp_st2(propagatedPositionSt2());
  mindREMTFTrack = 10.0;
  if (verboseEMTFTrack_)
      std::cout <<"propaget position to st2 eta "<< float(gp_st2.eta()) <<" phi "<< float(gp_st2.phi()) << std::endl;
  for (const auto& trk : tracks) {
    if (verboseEMTFTrack_)
	std::cout <<"track BX "<< trk.BX() <<  " pt "<< trk.Pt() <<" eta "<< trk.Eta() <<" phi "<< trk.Phi_glob_rad()<< std::endl;
    if (trk.BX() < minBXEMTFTrack_ or trk.BX() > maxBXEMTFTrack_) continue;
    float dR = 10.0;
    dR = deltaR(float(gp_st2.eta()), float(gp_st2.phi()), trk.Eta(), trk.Phi_glob_rad());
    if (verboseEMTFTrack_)
	std::cout <<"dR (track, sim) "<< dR <<" deltaREMTFTrack_ "<< deltaREMTFTrack_ << std::endl;
    if (dR < deltaREMTFTrack_){
	TFTrack* track  = new TFTrack(&trk);
	track->setDR(dR);
	if (verboseEMTFTrack_){
	    std::cout <<"bestTrack cand "<< std::endl;
	    track->print();
	}
	tfTracks_.push_back(track);
	if (dR < mindREMTFTrack)
	    bestTrack = track;
    }
  //   // check the matching CSC stubs
  //   const auto sim_stubs = csc_stub_matcher_->allLctsMatched2SimMuon();
  //   const l1t::EMTFHitCollection l1_stubs = trk.Hits();
  //   for (const auto& l1_stub: l1_stubs){
  //     CSCCorrelatedLCTDigi csc_stub = l1_stub.CSC_LCTDigi();

  //     for (const auto& sim_stub: sim_stubs){
  //       if (csc_stub==sim_stub.
  //     }
   }
}

void UpgradeL1TrackMatcher::matchGMTToSimTrack(const BXVector<l1t::RegionalMuonCand>& gmtCands)
{
   if (tfTracks_.size()  ==  0) return;
   bestTrack->print();
   for (int bx = gmtCands.getFirstBX(); bx != gmtCands.getLastBX(); bx++ ){
       if ( bx < minBXGMT_ or bx > maxBXGMT_) continue;
       for (std::vector<l1t::RegionalMuonCand>::const_iterator cand = gmtCands.begin(bx); cand != gmtCands.end(bx); ++cand ){
	   float pt = cand->hwPt() * 0.5;
	   float phi = cand->hwPhi() * 2.0/576.0 ;
	   float eta = cand->hwEta() * 0.010875;
	   int charge = cand->hwSign() == 0? 0 : -1;
	   //int quality  = cand->hwQual();
	   std::cout <<"RegionalMuonCand hwPt "<< cand->hwPt()<<" pt "<< pt <<" hwPhi " << cand->hwPhi() <<" phi "<< phi <<" hwEta "<< cand->hwEta() <<" eta "<< eta <<" charge "<< charge << std::endl;
       	   
       }
   }

}

