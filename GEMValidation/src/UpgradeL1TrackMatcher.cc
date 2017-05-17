#include "GEMCode/GEMValidation/interface/UpgradeL1TrackMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TLorentzVector.h"
#include <map>

UpgradeL1TrackMatcher::UpgradeL1TrackMatcher(CSCStubMatcher& csc,
                                             edm::EDGetTokenT<l1t::EMTFTrackCollection> &emtfTrackInputLabel_,
                                             edm::EDGetTokenT< BXVector<l1t::RegionalMuonCand> > & regionalMuonCandInputLabel_,
                                             edm::EDGetTokenT< BXVector<l1t::Muon> > & gmtInputLabel_)
  : BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
  , csc_stub_matcher_(&csc)
{
  auto tfTrack = conf().getParameter<edm::ParameterSet>("upgradeEmtfTrack");
  minBXEMTFTrack_ = tfTrack.getParameter<int>("minBX");
  maxBXEMTFTrack_ = tfTrack.getParameter<int>("maxBX");
  verboseEMTFTrack_ = tfTrack.getParameter<int>("verbose");
  deltaREMTFTrack_ = tfTrack.getParameter<double>("deltaR");

  auto regionalMuonCand = conf().getParameter<edm::ParameterSet>("upgradeEmtfCand");
  minBXRegMuCand_ = regionalMuonCand.getParameter<int>("minBX");
  maxBXRegMuCand_ = regionalMuonCand.getParameter<int>("maxBX");
  verboseRegMuCand_ = regionalMuonCand.getParameter<int>("verbose");
  deltaRRegMuCand_ = regionalMuonCand.getParameter<double>("deltaR");

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

  edm::Handle<BXVector<l1t::RegionalMuonCand>> hRegMuonCand;
  if (gemvalidation::getByToken(regionalMuonCandInputLabel_,hRegMuonCand, event())) matchRegionalMuonCandToSimTrack(*hRegMuonCand.product());
  else 
      std::cout  <<"failed readout RegionalMuonCand " << std::endl;

  edm::Handle<BXVector<l1t::Muon>> hGMT;
  if (gemvalidation::getByToken(gmtInputLabel_,hGMT, event())) matchGMTToSimTrack(*hGMT.product());
  else 
      std::cout  <<"failed readout GMT " << std::endl;
}

UpgradeL1TrackMatcher::~UpgradeL1TrackMatcher()
{
}

void
UpgradeL1TrackMatcher::clear()
{
  bestTrack = NULL;
  bestGMT = NULL;
  tfTracks_.clear();
}

void
UpgradeL1TrackMatcher::matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection& tracks)
{
  GlobalPoint gp_st2(propagatedPositionSt2());
  mindREMTFTrack = deltaREMTFTrack_;
  if (verboseEMTFTrack_)
      std::cout <<"propaget position to st2 eta "<< float(gp_st2.eta()) <<" phi "<< float(gp_st2.phi()) << std::endl;
  for (const auto& trk : tracks) {
    if (verboseEMTFTrack_)
	std::cout <<"track BX "<< trk.BX() <<  " pt "<< trk.Pt() <<" eta "<< trk.Eta() <<" phi "<< trk.Phi_glob_rad()<<" phi_local "<< trk.Phi_loc_rad() << std::endl;
    if (trk.BX() < minBXEMTFTrack_ or trk.BX() > maxBXEMTFTrack_) continue;
    float dR = 10.0;
    dR = deltaR(float(gp_st2.eta()), float(gp_st2.phi()), trk.Eta(), trk.Phi_glob_rad());
    if (verboseEMTFTrack_)
	std::cout <<"dR (track, sim) "<< dR <<" deltaREMTFTrack_ "<< deltaREMTFTrack_ << std::endl;
    if (dR < deltaREMTFTrack_){
	TFTrack* track  = new TFTrack(&trk);
	track->setDR(dR);
	if (verboseEMTFTrack_){
	    track->print();
	}
	tfTracks_.push_back(track);
	if (dR < mindREMTFTrack){
	    mindREMTFTrack = dR;
	    bestTrack = track;
	}
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
   if (verboseGMT_ and bestTrack){
       std::cout <<"all matched TFTRack size "<< tfTracks_.size() << std::endl;
       std::cout <<"best TFTrack ";  bestTrack->print();
   }
}

void UpgradeL1TrackMatcher::matchRegionalMuonCandToSimTrack(const BXVector<l1t::RegionalMuonCand>& regMuCands)
{
   if (tfTracks_.size()  ==  0) return;
   float mindPtRel = 0.5;
   mindRRegMuCand = deltaRRegMuCand_;
   for (int bx = regMuCands.getFirstBX(); bx <= regMuCands.getLastBX(); bx++ ){
       if ( bx < minBXRegMuCand_ or bx > maxBXRegMuCand_) continue;
       for (std::vector<l1t::RegionalMuonCand>::const_iterator cand = regMuCands.begin(bx); cand != regMuCands.end(bx); ++cand ){
	   TFCand *L1Mu = new TFCand(&(*cand));
	   L1Mu->setBx(bx);
	   float pt = L1Mu->pt();
	   float phi = L1Mu->phi_local() ;
	   float eta = L1Mu->eta();
           for (auto trk : tfTracks_){
	       float dR = deltaR(trk->eta(), trk->phi_local(), eta, phi);
	       float dPtRel = std::fabs(trk->pt() - pt)/pt;
	       if (dR < deltaRRegMuCand_ and dPtRel < mindPtRel){
		   L1Mu->setDR( dR );
		   L1Mu->setGlobalPhi(trk->phi());
		   L1Mu->setMatchedTFTrack( trk );
		   regMuCands_.push_back(L1Mu);
	       }
	   }
	   if (verboseRegMuCand_)
	       L1Mu->print();
       }
   }
   for (auto cand : regMuCands_){
       float phi = cand->phi_local();
       float eta = cand->eta();
       float dR = deltaR(bestTrack->eta(), bestTrack->phi_local(), eta, phi);
       if (dR < mindRRegMuCand){
	   mindRRegMuCand = dR;
	   bestRegMuCand = cand;
	   if (verboseRegMuCand_){
	       std::cout <<"bestRegMuCand "; bestRegMuCand->print();
	   }
       }
   }

}

void UpgradeL1TrackMatcher::matchGMTToSimTrack(const BXVector<l1t::Muon>& gmtCands)
{
   std::cout <<"matchGMTToSimTrack, tftracks size "<< tfTracks_.size() << std::endl;
   if (tfTracks_.size()  ==  0) return;
   float mindPtRel = 0.5;
   mindRGMT = deltaRGMT_;
   for (int bx = gmtCands.getFirstBX(); bx <= gmtCands.getLastBX(); bx++ ){
       std::cout <<"matching L1Mu to EMTF track,  bx "<< bx << std::endl;
       if ( bx < minBXGMT_ or bx > maxBXGMT_) continue;
       for (std::vector<l1t::Muon>::const_iterator cand = gmtCands.begin(bx); cand != gmtCands.end(bx); ++cand ){
	   TFCand *L1Mu = new TFCand(&(*cand));
	   L1Mu->setBx(bx);
	   float pt = L1Mu->pt();
	   float phi = L1Mu->phi() ;
	   float eta = L1Mu->eta();
           for (auto trk : tfTracks_){
	       float dR = deltaR(trk->eta(), trk->phi(), eta, phi);
	       float dPtRel = std::fabs(trk->pt() - pt)/pt;
	       if (dR < deltaRGMT_ and dPtRel < mindPtRel){
		   L1Mu->setDR( dR );
		   L1Mu->setMatchedTFTrack( trk );
		   gmts_.push_back(L1Mu);
	       }
	   }
	   if (verboseGMT_)
	       L1Mu->print();
       }
   }
   for (auto cand : gmts_){
       float phi = cand->phi();
       float eta = cand->eta();
       float dR = deltaR(bestTrack->eta(), bestTrack->phi(), eta, phi);
       if (dR < mindRGMT){
	   mindRGMT = dR;
	   bestGMT = cand;
	   if (verboseGMT_){
	       std::cout <<"bestGMT "; bestGMT->print();
	   }
       }
   }

}

