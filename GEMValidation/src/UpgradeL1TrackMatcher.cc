#include "GEMCode/GEMValidation/interface/UpgradeL1TrackMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TLorentzVector.h"
#include <map>

UpgradeL1TrackMatcher::UpgradeL1TrackMatcher(CSCStubMatcher& csc,
                                             edm::EDGetTokenT<l1t::EMTFTrackCollection> &emtfTrackInputLabel_)
  : BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
  , csc_stub_matcher_(&csc)
{
  auto tfTrack = conf().getParameter<edm::ParameterSet>("cscTfTrack");
  minBXEMTFTrack_ = tfTrack.getParameter<int>("minBX");
  maxBXEMTFTrack_ = tfTrack.getParameter<int>("minBX");
  verboseEMTFTrack_ = tfTrack.getParameter<int>("verbose");
  deltaREMTFTrack_ = tfTrack.getParameter<double>("deltaR");

  //std::cout<<" UpgradeL1TrackMatcher constructor" <<std::endl;
  clear();
  init();

  // tracks produced by EMEMTF
  edm::Handle<l1t::EMTFTrackCollection> hl1Tracks;
  if (gemvalidation::getByToken(emtfTrackInputLabel_,hl1Tracks, event())) matchEmtfTrackToSimTrack(*hl1Tracks.product());
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
UpgradeL1TrackMatcher::init()
{
}


void
UpgradeL1TrackMatcher::matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection& tracks)
{
  for (auto& trk : tracks) {

    // check the matching CSC stubs
    auto sim_stubs = csc_stub_matcher_->allLctsMatched2SimMuon();
    // EMTFHitCollection l1_stubs = trk.Hits();
  }
  /*
  const float simPt = trk().momentum().pt();
  const float simEta = trk().momentum().eta();
  const float simPhi = trk().momentum().phi();
  const float simE = trk().momentum().E();
  const float simCharge = trk().charge();

  for (const auto& trk = tracks.begin(); trk != tracks.end(); trk++) {
     verboseTFTrack_ = 0;
    TFTrack *track = new TFTrack(&trk->first,&trk->second);
    //track->init(ptLUT_, muScalesHd_, muPtScaleHd_);
    track->init(muScalesHd_, muPtScaleHd_);

    float dR = 10.0;
    TLorentzVector templ1muon;
    templ1muon.SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), 0.1057);


    if (simEta*track->eta() < 0) continue;
    // calculate the deltaR using the simhits in the 2nd CSC station -- reference station
    auto p(rpc_digi_matcher_->simHitMatcher()->chamberIdsCSC(CSC_ME21));
    auto pp(rpc_digi_matcher_->simHitMatcher()->chamberIdsCSC(CSC_ME22));
    p.insert(pp.begin(),pp.end());

    TLorentzVector simmuon;

    if (verboseTFTrack_ > 1) std::cout << "----------------------------------------------------------" << std::endl;
    if (verboseTFTrack_ > 1) std::cout << "detids " << std::endl;

    if (p.size()!=0) {
      if (verboseTFTrack_ > 1) std::cout << CSCDetId(*p.begin()) << std::endl;
      auto hits(rpc_digi_matcher_->simHitMatcher()->hitsInChamber(*p.begin()));
      auto gp(rpc_digi_matcher_->simHitMatcher()->simHitsMeanPosition(hits));
      if (verboseTFTrack_ > 1) std::cout << gp << std::endl;
      // simmuon.SetPtEtaPhiM(simPt, gp.eta(), gp.phi(), 1.057);
    }

    //propagate simtrack to station2
    float phi_cor = phiHeavyCorr(simPt, simEta, simPhi, simCharge);
    simmuon.SetPtEtaPhiM(simPt, simEta, phi_cor, 0.1057);
    dR = simmuon.DeltaR(templ1muon);

   // auto lcts2(lctsInStation(2));
   //  auto lcts3(lctsInStation(3));
   //  lcts2.insert(lcts2.end(),lcts3.begin(),lcts3.end());

   //  if (lcts2.size() != 0) {

   // 	TLorentzVector simmuon;
   //      for (auto stub : lcts2)
   //      {
   // 	    if (verboseTFTrack_ > 1)   std::cout <<"stubs in st2,3 "<< stub << std::endl;
   // 	    auto EtaPhi(intersectionEtaPhi(digi_id(stub), digi_wg(stub), digi_channel(stub)));
   //          //simmuon.SetPtEtaPhiE(simPt, simEta, simPhi, simE);
   //          simmuon.SetPtEtaPhiE(simPt, EtaPhi.first, EtaPhi.second, simE);
   // 	    if (simmuon.DeltaR(templ1muon) < dR)  dR = simmuon.DeltaR(templ1muon);
   // 	}

   //  }
   //  else continue;//if no stubs available in station2,or 3 that can match to simtrack, then continue;???

    if (verboseTFTrack_ > 1) std::cout << "----------------------------------------------------------" << std::endl;


    track->setDR(dR);

    // debugging
    if (track->nStubs()==0)
    {
        std::cout << "nstubs == 0" << std::endl;
	verboseTFTrack_ = 2;
    }
    //debug
    //if (!(track->hasStubEndcap(1) and track->hasStubEndcap(2)) and track->dr()< deltaRTFTrack_)
   // {
    //     std::cout <<"no stubs in station 1 or 2" << std::endl;
//	 verboseTFTrack_ = 2;
  //  }

    if (verboseTFTrack_){
      std::cout << "\tL1CSC TFTrack "<< trk-tracks.begin() << " information:" << std::endl;
      std::cout << "\tpt (GeV/c) = " << track->pt() << ", eta = " << track->eta()
                << "\t, phi = " << track->phi() << ", dR(sim,L1) = " << track->dr()
		<<" nStubs = "<< track->nStubs() << ", deltaphi12 = "<< track->dPhi12() <<", deltaphi23 = "<<track->dPhi23() <<std::endl;
      std::cout << " pt_packed " << track->ptPacked()  << " eta_packed " << track->etaPacked() << " phi_packed " << track->phiPacked() << std::endl;
      std::cout << "simTrack \t simpt " << simPt << " simeta "<< simEta << " simPhi "<< simPhi <<" simenergy "<< simE << std::endl;
    }

    // check the stubs
    if (verboseTFTrack_ > 1){
      std::cout << " \t\tStubs:" << std::endl;
    }
    for (auto  detUnitIt = trk->second.begin();
         detUnitIt != trk->second.end(); detUnitIt++) {
      const CSCDetId& id = (*detUnitIt).first;
      if (verboseTFTrack_ > 1){
        std::cout << "\t\tDetId: " << id << std::endl;
      }
      const CSCCorrelatedLCTDigiCollection::Range& range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        auto lct(*digiIt);
        // check for valid stubs
        if (!(lct.isValid())) continue;
	//  track.addTriggerDigi(&lct);
	//  track.addTriggerDigiId(id);
	auto EtaPhi(intersectionEtaPhi(id, lct.getKeyWG(), lct.getStrip()));
        track->addTriggerEtaPhi(EtaPhi);
        track->addTriggerStub(buildTrackStub(lct, id));
        if (verboseTFTrack_ > 1 ) {
	  //    auto EtaPhi(intersectionEtaPhi(id, lct.getKeyWG(), lct.getStrip()));
          std::cout << "\t\tLCT" << digiIt-range.first<<" eta:"<< EtaPhi.first
		    <<" phi:"<< EtaPhi.second << ": " << lct << std::endl;
        }
      }
    }

    //check the dR of simtrack and tftrack
    if (track->dr() < deltaRTFTrack_) {
      tfTracks_.push_back(track);
      if (verboseTFTrack_) std::cout<<" dR(sim, L1) is ok "<<" deltaRTFTrack"<< deltaRTFTrack_ <<" real dr"<<track->dr()<< std::endl;
    }

  }
  int i=0;
  for (auto tftrk : tfTracks_)
    {
      if (verboseTFTrack_)  std::cout<<" track that can match to simtrack "<< i << " track pt " << tftrk->pt() <<" track eta " << tftrk->eta()
				     <<" nstubs:"<< tftrk->nStubs()<<" dr:" << tftrk->dr() <<std::endl;
      i++;
    }
*/
}

