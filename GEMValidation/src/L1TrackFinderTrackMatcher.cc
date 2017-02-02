#include "GEMCode/GEMValidation/interface/L1TrackFinderTrackMatcher.h"

L1TrackFinderTrackMatcher::L1TrackFinderTrackMatcher(SimHitMatcher& sh,
                                                     edm::EDGetTokenT<L1CSCTrackCollection> &cscTfTrackInputLabel_)
: BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
{
  auto cscTfTrack = conf().getParameter<edm::ParameterSet>("cscTfTrack");
  auto dtTfTrack = conf().getParameter<edm::ParameterSet>("dtTfTrack");
  auto rpcTfTrack = conf().getParameter<edm::ParameterSet>("rpcTfTrack");

  //dtTfTrackInputLabel_ = dtTfTrack.getParameter<std::vector<edm::InputTag>>("validInputTags");
  //rpcTfTrackInputLabel_ = rpcTfTrack.getParameter<std::vector<edm::InputTag>>("validInputTags");
  
  verboseCscTfTrack_ = cscTfTrack.getParameter<int>("verbose");
  verboseDtTfTrack_ = dtTfTrack.getParameter<int>("verbose");
  verboseRpcTfTrack_ = rpcTfTrack.getParameter<int>("verbose");

  runCscTfTrack_ = cscTfTrack.getParameter<bool>("run");
  runDtTfTrack_ = dtTfTrack.getParameter<bool>("run");
  runRpcTfTrack_ = rpcTfTrack.getParameter<bool>("run");

  minBXCscTfTrack_ = cscTfTrack.getParameter<int>("minBX");
  minBXDtTfTrack_ = dtTfTrack.getParameter<int>("minBX");
  minBXRpcTfTrack_ = rpcTfTrack.getParameter<int>("minBX");

  maxBXCscTfTrack_ = cscTfTrack.getParameter<int>("maxBX");
  maxBXDtTfTrack_ = dtTfTrack.getParameter<int>("maxBX");
  maxBXRpcTfTrack_ = rpcTfTrack.getParameter<int>("maxBX");

  edm::Handle<L1CSCTrackCollection> hCscTfTrack;
  if (runCscTfTrack_) if (gemvalidation::getByToken(cscTfTrackInputLabel_, hCscTfTrack, event())) matchCSCTfTrackToSimTrack(*hCscTfTrack.product());

  /*
  edm::Handle<L1CSCTrackCollection> hDtTfTrack;
  if (runDtTfTrack_) if (gemvalidation::getByToken(dtTfTrackInputLabel_, hDtTfTrack, event())) matchDTTfTrackToSimTrack(*hDtTfTrack.product());

  edm::Handle<L1CSCTrackCollection> hRpcTfTrack;
  if (runRpcTfTrack_) if (gemvalidation::getByToken(rpcTfTrackInputLabel_, hRpcTfTrack, event())) matchRPCTfTrackToSimTrack(*hRpcTfTrack.product());
   */
  init();
}

L1TrackFinderTrackMatcher::~L1TrackFinderTrackMatcher()
{}

void 
L1TrackFinderTrackMatcher::init()
{

}

void 
L1TrackFinderTrackMatcher::matchCSCTfTrackToSimTrack(const L1CSCTrackCollection& tracks)
{

  return;
  for (auto trk : tracks) {
    auto tfTrack(trk.first);
    std::cout << "Stubs in CSCTF Track" << std::endl;
    for (auto  detUnitIt = trk.second.begin(); detUnitIt != trk.second.end(); detUnitIt++) {
      const CSCDetId& id = (*detUnitIt).first;
      // if (verboseTFTrack_ > 1){
      std::cout << "\tDetId: " << id << std::endl;
      // }
      const CSCCorrelatedLCTDigiCollection::Range& range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
	auto lct(*digiIt);
	std::cout << "\t\tStub: " << lct << std::endl;
      }
    }
    //    matchedTfTracks_.push_back(tfTrack)
  }
  /*
  const float simPt = trk().momentum().pt();
  const float simEta = trk().momentum().eta();
  const float simPhi = trk().momentum().phi();
  const float simE = trk().momentum().E();
  const float simCharge = trk().charge();
  
  for (auto trk = tracks.begin(); trk != tracks.end(); trk++) {
    verboseTFTrack_ = 0;
    TFTrack *track = new TFTrack(&trk->first,&trk->second);
    track->init(ptLUT_, muScalesHd_, muPtScaleHd_);
    
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
    
    // TLorentzVector simmuon;
    //      for (auto stub : lcts2)
    //      {
    //     if (verboseTFTrack_ > 1)   std::cout <<"stubs in st2,3 "<< stub << std::endl;
    //     auto EtaPhi(intersectionEtaPhi(digi_id(stub), digi_wg(stub), digi_channel(stub)));
    //          //simmuon.SetPtEtaPhiE(simPt, simEta, simPhi, simE);
    //          simmuon.SetPtEtaPhiE(simPt, EtaPhi.first, EtaPhi.second, simE);
    //     if (simmuon.DeltaR(templ1muon) < dR)  dR = simmuon.DeltaR(templ1muon);
    // }
    
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
    // verboseTFTrack_ = 2;
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
      const CSCCorrelated
  */
}

void 
L1TrackFinderTrackMatcher::matchDTTfTrackToSimTrack(const L1CSCTrackCollection& tracks)
{}

void 
L1TrackFinderTrackMatcher::matchRPCTfTrackToSimTrack(const L1CSCTrackCollection& tracks)
{}

