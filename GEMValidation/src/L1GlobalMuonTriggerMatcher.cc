#include "GEMCode/GEMValidation/interface/L1GlobalMuonTriggerMatcher.h"

L1GlobalMuonTriggerMatcher::L1GlobalMuonTriggerMatcher(SimHitMatcher& sh)
: BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
{
  auto gmtRegCandCSC = conf().getParameter<edm::ParameterSet>("gmtRegCandCSC");
  auto gmtRegCandDT = conf().getParameter<edm::ParameterSet>("gmtRegCandDT");
  auto gmtRegCandRPCf = conf().getParameter<edm::ParameterSet>("gmtRegCandRPCf");
  auto gmtRegCandRPCb = conf().getParameter<edm::ParameterSet>("gmtRegCandRPCb");
  auto gmtCand = conf().getParameter<edm::ParameterSet>("gmtCand");
  auto l1ExtraMuonParticle = conf().getParameter<edm::ParameterSet>("l1ExtraMuonParticle");

  gmtRegCandCSCInputLabel_ = gmtRegCandCSC.getParameter<std::vector<edm::InputTag>>("validInputTags");
  gmtRegCandDTInputLabel_ = gmtRegCandDT.getParameter<std::vector<edm::InputTag>>("validInputTags");
  gmtRegCandRPCfInputLabel_ = gmtRegCandRPCf.getParameter<std::vector<edm::InputTag>>("validInputTags");
  gmtRegCandRPCbInputLabel_ = gmtRegCandRPCb.getParameter<std::vector<edm::InputTag>>("validInputTags");
  gmtCandInputLabel_ = gmtCand.getParameter<std::vector<edm::InputTag>>("validInputTags");
  l1ExtraMuonInputLabel_ = l1ExtraMuonParticle.getParameter<std::vector<edm::InputTag>>("validInputTags");

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

  edm::Handle<L1MuGMTExtendedCandCollection> hGmtCand;
  if (gemvalidation::getByLabel(gmtCandInputLabel_, hGmtCand, event())) if (runGmtCand_) matchGMTCandToSimTrack(*hGmtCand.product());

  edm::Handle<l1extra::L1MuonParticleCollection> hL1ExtraMuonParticle;
  if (gemvalidation::getByLabel(l1ExtraMuonInputLabel_, hL1ExtraMuonParticle, event())) if (runL1ExtraMuon_) matchL1ExtraMuonParticleToSimTrack(*hL1ExtraMuonParticle.product());
}

void L1GlobalMuonTriggerMatcher::matchRegionalCandCSCToSimTrack(const L1MuRegionalCandCollection&)
{}

void L1GlobalMuonTriggerMatcher::matchRegionalCandDTToSimTrack(const L1MuRegionalCandCollection&) 
{}

void L1GlobalMuonTriggerMatcher::matchRegionalCandRPCbToSimTrack(const L1MuRegionalCandCollection&) 
{}

void L1GlobalMuonTriggerMatcher::matchRegionalCandRPCfToSimTrack(const L1MuRegionalCandCollection&) 
{}

void L1GlobalMuonTriggerMatcher::matchGMTCandToSimTrack(const L1MuGMTExtendedCandCollection&) 
{}

void L1GlobalMuonTriggerMatcher::matchL1ExtraMuonParticleToSimTrack(const l1extra::L1MuonParticleCollection&) 
{}
