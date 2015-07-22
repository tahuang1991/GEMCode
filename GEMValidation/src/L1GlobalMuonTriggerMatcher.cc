#include "GEMCode/GEMValidation/interface/L1GlobalMuonTriggerMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

using namespace std;

L1GlobalMuonTriggerMatcher::L1GlobalMuonTriggerMatcher(SimHitMatcher& sh)
: BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
, simhit_matcher_(&sh)
{
  auto gmtRegCandCSC = conf().getParameter<edm::ParameterSet>("gmtRegCandCSC");
  auto gmtRegCandDT = conf().getParameter<edm::ParameterSet>("gmtRegCandDT");
  auto gmtRegCandRPCf = conf().getParameter<edm::ParameterSet>("gmtRegCandRPCf");
  auto gmtRegCandRPCb = conf().getParameter<edm::ParameterSet>("gmtRegCandRPCb");
  auto gmtCand = conf().getParameter<edm::ParameterSet>("gmtCand");
  auto l1ExtraMuonParticle = conf().getParameter<edm::ParameterSet>("l1ExtraMuonParticle");

  gmtRegCandCSCInputLabel_ = gmtRegCandCSC.getParameter<vector<edm::InputTag>>("validInputTags");
  gmtRegCandDTInputLabel_ = gmtRegCandDT.getParameter<vector<edm::InputTag>>("validInputTags");
  gmtRegCandRPCfInputLabel_ = gmtRegCandRPCf.getParameter<vector<edm::InputTag>>("validInputTags");
  gmtRegCandRPCbInputLabel_ = gmtRegCandRPCb.getParameter<vector<edm::InputTag>>("validInputTags");
  gmtCandInputLabel_ = gmtCand.getParameter<vector<edm::InputTag>>("validInputTags");
  l1ExtraMuonInputLabel_ = l1ExtraMuonParticle.getParameter<vector<edm::InputTag>>("validInputTags");

  verboseGmtRegCandCSC_ = gmtRegCandCSC.getParameter<int>("verbose");
  verboseGmtRegCandDT_ = gmtRegCandDT.getParameter<int>("verbose");
  verboseGmtRegCandRPCf_ = gmtRegCandRPCf.getParameter<int>("verbose");
  verboseGmtRegCandRPCb_ = gmtRegCandRPCb.getParameter<int>("verbose");
  verboseGmtCand_ = gmtCand.getParameter<int>("verbose");
  verboseL1ExtraMuon_ = l1ExtraMuonParticle.getParameter<int>("verbose");

  runGmtRegCandCSC_ = gmtRegCandCSC.getParameter<bool>("run");
  runGmtRegCandDT_ = gmtRegCandDT.getParameter<bool>("run");
  runGmtRegCandRPCf_ = gmtRegCandRPCf.getParameter<bool>("run");
  runGmtRegCandRPCb_ = gmtRegCandRPCb.getParameter<bool>("run");
  runGmtCand_ = gmtCand.getParameter<bool>("run");
  runL1ExtraMuon_ = l1ExtraMuonParticle.getParameter<bool>("run");

  minBXGmtRegCandCSC_ = gmtRegCandCSC.getParameter<int>("minBX");
  minBXGmtRegCandDT_ = gmtRegCandDT.getParameter<int>("minBX");
  minBXGmtRegCandRPCf_ = gmtRegCandRPCf.getParameter<int>("minBX");
  minBXGmtRegCandRPCb_ = gmtRegCandRPCb.getParameter<int>("minBX");
  minBXGmtCand_ = gmtCand.getParameter<int>("minBX");
  minBXL1ExtraMuon_ = l1ExtraMuonParticle.getParameter<int>("minBX");

  maxBXGmtRegCandCSC_ = gmtRegCandCSC.getParameter<int>("maxBX");
  maxBXGmtRegCandDT_ = gmtRegCandDT.getParameter<int>("maxBX");
  maxBXGmtRegCandRPCf_ = gmtRegCandRPCf.getParameter<int>("maxBX");
  maxBXGmtRegCandRPCb_ = gmtRegCandRPCb.getParameter<int>("maxBX");
  maxBXGmtCand_ = gmtCand.getParameter<int>("maxBX");
  maxBXL1ExtraMuon_ = l1ExtraMuonParticle.getParameter<int>("maxBX");

  deltaRGmtRegCandCSC_ = gmtRegCandCSC.getParameter<double>("deltaR");
  deltaRGmtRegCandDT_ = gmtRegCandDT.getParameter<double>("deltaR");
  deltaRGmtRegCandRPCf_ = gmtRegCandRPCf.getParameter<double>("deltaR");
  deltaRGmtRegCandRPCb_ = gmtRegCandRPCb.getParameter<double>("deltaR");
  deltaRGmtCand_ = gmtCand.getParameter<double>("deltaR");
  deltaRL1ExtraMuon_ = l1ExtraMuonParticle.getParameter<double>("deltaR");

  init();
}

L1GlobalMuonTriggerMatcher::~L1GlobalMuonTriggerMatcher()
{}


void 
L1GlobalMuonTriggerMatcher::clear()
{}


void 
L1GlobalMuonTriggerMatcher::init()
{
  edm::Handle<L1MuRegionalCandCollection> hGmtRegCandCSC;
  if (gemvalidation::getByLabel(gmtRegCandCSCInputLabel_, hGmtRegCandCSC, event())) if (runGmtRegCandCSC_) matchRegionalCandCSCToSimTrack(*hGmtRegCandCSC.product());

  edm::Handle<L1MuRegionalCandCollection> hGmtRegCandRPCf;
  if (gemvalidation::getByLabel(gmtRegCandRPCfInputLabel_, hGmtRegCandRPCf, event())) if (runGmtRegCandRPCf_) matchRegionalCandRPCfToSimTrack(*hGmtRegCandRPCf.product());

  edm::Handle<L1MuRegionalCandCollection> hGmtRegCandRPCb;
  if (gemvalidation::getByLabel(gmtRegCandRPCbInputLabel_, hGmtRegCandRPCb, event())) if (runGmtRegCandRPCb_) matchRegionalCandRPCbToSimTrack(*hGmtRegCandRPCb.product());

  edm::Handle<L1MuRegionalCandCollection> hGmtRegCandDT;
  if (gemvalidation::getByLabel(gmtRegCandDTInputLabel_, hGmtRegCandDT, event())) if (runGmtRegCandDT_) matchRegionalCandDTToSimTrack(*hGmtRegCandDT.product());

  edm::Handle<L1MuGMTCandCollection> hGmtCand;
  if (gemvalidation::getByLabel(gmtCandInputLabel_, hGmtCand, event())) if (runGmtCand_) matchGMTCandToSimTrack(*hGmtCand.product());

  edm::Handle<l1extra::L1MuonParticleCollection> hL1ExtraMuonParticle;
  if (gemvalidation::getByLabel(l1ExtraMuonInputLabel_, hL1ExtraMuonParticle, event())) if (runL1ExtraMuon_) matchL1ExtraMuonParticleToSimTrack(*hL1ExtraMuonParticle.product());
}

void 
L1GlobalMuonTriggerMatcher::matchRegionalCandCSCToSimTrack(const L1MuRegionalCandCollection& cands)
{
  if (verboseGmtRegCandCSC_) cout << "Match SimTrack to CSC GMTCands" << endl;
  int i=0;
  for (auto& cand: cands) {
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), cand.etaValue(), normalizedPhi(cand.phiValue())));
    if (verboseGmtRegCandCSC_) {
      cout << i+1 << ": pT = " << cand.ptValue() << ", eta = " << cand.etaValue() <<  ", phi = " << normalizedPhi(cand.phiValue())
	   << ", bx = " << cand.bx() << ", charge = " << cand.chargeValue() << ", quality = " << cand.quality() << endl;
      cout << "\tDeltaR = " << dR << endl;
    }
    // BX accept
    if (std::abs(cand.bx()) > 1) continue;
    // DeltaR accept
    if (dR > deltaRGmtRegCandCSC_) continue;
    matchedL1GmtCSCCands_.push_back(cand);
    ++i;
  }
}

void 
L1GlobalMuonTriggerMatcher::matchRegionalCandDTToSimTrack(const L1MuRegionalCandCollection& cands) 
{
  if (verboseGmtRegCandDT_) cout << "Match SimTrack to DT GMTCands" << endl;
  int i=0;
  for (auto& cand: cands) {
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), cand.etaValue(), normalizedPhi(cand.phiValue())));
    if (verboseGmtRegCandDT_) {
      cout << i+1 << ": pT = " << cand.ptValue() << ", eta = " << cand.etaValue() <<  ", phi = " << normalizedPhi(cand.phiValue())
	   << ", bx = " << cand.bx() << ", charge = " << cand.chargeValue() << ", quality = " << cand.quality() << endl;
      cout << "\tDeltaR = " << dR << endl;
    }
    // BX accept
    if (std::abs(cand.bx()) > 1) continue;
    // DeltaR accept
    if (dR > deltaRGmtRegCandDT_) continue;
    matchedL1GmtDTCands_.push_back(cand);
    ++i;
  }
}

void 
L1GlobalMuonTriggerMatcher::matchRegionalCandRPCbToSimTrack(const L1MuRegionalCandCollection& cands) 
{
  if (verboseGmtRegCandRPCb_) cout << "Match SimTrack to RPCb GMTCands" << endl;
  int i=0;
  for (auto& cand: cands) {
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), cand.etaValue(), normalizedPhi(cand.phiValue())));
    if (verboseGmtRegCandRPCb_) {
      cout << i+1 << ": pT = " << cand.ptValue() << ", eta = " << cand.etaValue() <<  ", phi = " << normalizedPhi(cand.phiValue())
	   << ", bx = " << cand.bx() << ", charge = " << cand.chargeValue() << ", quality = " << cand.quality() << endl;
      cout << "\tDeltaR = " << dR << endl;   
    }
    // BX accept
    if (std::abs(cand.bx()) > 1) continue;
    // DeltaR accept
    if (dR > deltaRGmtRegCandRPCb_) continue;
    matchedL1GmtRPCbCands_.push_back(cand);
    ++i;
  }
}

void 
L1GlobalMuonTriggerMatcher::matchRegionalCandRPCfToSimTrack(const L1MuRegionalCandCollection& cands) 
{
  if (verboseGmtRegCandRPCf_) cout << "Match SimTrack to RPCf GMTCands" << endl;
  int i=0;
  for (auto& cand: cands) {
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), cand.etaValue(), normalizedPhi(cand.phiValue())));
    if (verboseGmtRegCandRPCf_) {
      cout << i+1 << ": pT = " << cand.ptValue() << ", eta = " << cand.etaValue() <<  ", phi = " << normalizedPhi(cand.phiValue())
	   << ", bx = " << cand.bx() << ", charge = " << cand.chargeValue() << ", quality = " << cand.quality() << endl;
      cout << "\tDeltaR = " << dR << endl;
    }
    // BX accept
    if (std::abs(cand.bx()) > 1) continue;
    // DeltaR accept
    if (dR > deltaRGmtRegCandRPCf_) continue;
    matchedL1GmtRPCfCands_.push_back(cand);
    ++i;
  }
}

void 
L1GlobalMuonTriggerMatcher::matchGMTCandToSimTrack(const L1MuGMTCandCollection& cands) 
{
  if (verboseGmtCand_) cout << "Match SimTrack to GMTCands" << endl;
  int i=0;
  for (auto& cand: cands) {
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), cand.etaValue(), normalizedPhi(cand.phiValue())));
    if (verboseGmtCand_) {
      cout << i+1 << ": pT = " << cand.ptValue() << ", eta = " << cand.etaValue() <<  ", phi = " << normalizedPhi(cand.phiValue())
	   << ", bx = " << cand.bx() << ", charge = " << cand.charge() << ", quality = " << cand.quality() << endl;    
      cout << "\tDeltaR = " << dR << endl;
    }
    matchedL1GmtCands_.push_back(cand);
    ++i;
  }
}

void 
L1GlobalMuonTriggerMatcher::matchL1ExtraMuonParticleToSimTrack(const l1extra::L1MuonParticleCollection& muons) 
{
  if (verboseL1ExtraMuon_) cout << "Match SimTrack to L1ExtraMuonParticle" << endl;

  int i=0;
  int indexFound = -99;
  float dRhitMin = 99;
  for (auto& muon: muons) {
    auto gmt(muon.gmtMuonCand());    
    const float dR(deltaR(trk().momentum().eta(), trk().momentum().phi(), muon.eta(), normalizedPhi(muon.phi())));
    if (verboseL1ExtraMuon_) {
      cout << i+1 << ": pT = " << muon.pt() << ", eta = " << muon.eta() <<  ", phi = " << muon.phi()
	   << ", bx = " << muon.bx() << ", charge = " << muon.charge() << endl;    
      cout << "\tDeltaR = " << dR << endl;
      cout << "\tAssociated GMT " << gmt << endl;
    }

    const double absSimEta(std::abs(trk().momentum().eta()));
    float dRhit = 99;
    // propagate the simtrack to the second station DT or CSC and match to the L1Extra particle
    if (absSimEta < 1.1) {
      auto q(simhit_matcher_->chamberIdsDT(DT_MB22p));  
      auto qq(simhit_matcher_->chamberIdsDT(DT_MB12p));  
      auto qqq(simhit_matcher_->chamberIdsDT(DT_MB02));  
      auto qqqq(simhit_matcher_->chamberIdsDT(DT_MB12n));  
      auto qqqqq(simhit_matcher_->chamberIdsDT(DT_MB22n));  
      q.insert(qq.begin(),qq.end());
      q.insert(qqq.begin(),qqq.end());
      q.insert(qqqq.begin(),qqqq.end());
      q.insert(qqqqq.begin(),qqqqq.end());

      if (q.size()!=0) {
	auto hits(simhit_matcher_->hitsInChamber(*q.begin()));
	auto gp(simhit_matcher_->simHitsMeanPosition(hits));
	dRhit = deltaR(trk().momentum().eta(), trk().momentum().phi(), gp.eta(), normalizedPhi(gp.phi()));
	if (verboseL1ExtraMuon_>1) {
	  std::cout << "DTChamberId " << DTChamberId(*q.begin()) << std::endl;  
	  std::cout << "\t" << gp << std::endl;
	}
	if (verboseL1ExtraMuon_) std::cout << "\tDeltaR = " << dRhit << endl;
      }
    }
    else if (1.1 < absSimEta and absSimEta < 2.4) {
      auto p(simhit_matcher_->chamberIdsCSC(CSC_ME21));  
      auto pp(simhit_matcher_->chamberIdsCSC(CSC_ME22));  
      p.insert(pp.begin(),pp.end());

      if (p.size()!=0) {
	auto hits(simhit_matcher_->hitsInChamber(*p.begin()));
	auto gp(simhit_matcher_->simHitsMeanPosition(hits));
	dRhit = deltaR(trk().momentum().eta(), trk().momentum().phi(), gp.eta(), normalizedPhi(gp.phi()));
	if (verboseL1ExtraMuon_>1) {
	  std::cout << "CSCDetId " << CSCDetId(*p.begin()) << std::endl;  
	  std::cout << "\t" << gp << std::endl;
	}
	if (verboseL1ExtraMuon_) std::cout << "\tDeltaR = " << dRhit << endl;
      }
    }
    if (dRhit < dRhitMin) {
      dRhitMin = dRhit;
      indexFound = i;
    }
    ++i;
  }

  if (dRhitMin < deltaRL1ExtraMuon_) matchedL1MuonParticles_.push_back(muons[indexFound]);    
  if (verboseL1ExtraMuon_) std::cout << "\tdRhitMin = " << dRhitMin << endl;	
}


bool 
L1GlobalMuonTriggerMatcher::gmtCandInContainer(const L1MuGMTCand&, const L1MuGMTCandCollection&) const
{
  return false;
}

bool 
L1GlobalMuonTriggerMatcher::isGmtCandMatched(const L1MuGMTCand&) const
{
  return false;
}
