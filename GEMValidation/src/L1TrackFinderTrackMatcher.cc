#include "GEMCode/GEMValidation/interface/L1TrackFinderTrackMatcher.h"

L1TrackFinderTrackMatcher::L1TrackFinderTrackMatcher(SimHitMatcher& sh)
: BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
{
  auto cscTfTrack = conf().getParameter<edm::ParameterSet>("cscTfTrack");
  auto dtTfTrack = conf().getParameter<edm::ParameterSet>("dtTfTrack");
  auto rpcTfTrack = conf().getParameter<edm::ParameterSet>("rpcTfTrack");

  cscTfTrackInputLabel_ = cscTfTrack.getParameter<std::vector<edm::InputTag>>("validInputTags");
  dtTfTrackInputLabel_ = dtTfTrack.getParameter<std::vector<edm::InputTag>>("validInputTags");
  rpcTfTrackInputLabel_ = rpcTfTrack.getParameter<std::vector<edm::InputTag>>("validInputTags");
  
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

  init();
}

L1TrackFinderTrackMatcher::~L1TrackFinderTrackMatcher()
{}

void 
L1TrackFinderTrackMatcher::init()
{
  edm::Handle<L1CSCTrackCollection> hCscTfTrack;
  if (gemvalidation::getByLabel(cscTfTrackInputLabel_, hCscTfTrack, event())) if (runCscTfTrack_) matchCSCTfTrackToSimTrack(*hCscTfTrack.product());

  edm::Handle<L1CSCTrackCollection> hDtTfTrack;
  if (gemvalidation::getByLabel(dtTfTrackInputLabel_, hDtTfTrack, event())) if (runDtTfTrack_) matchDTTfTrackToSimTrack(*hDtTfTrack.product());

  edm::Handle<L1CSCTrackCollection> hRpcTfTrack;
  if (gemvalidation::getByLabel(rpcTfTrackInputLabel_, hRpcTfTrack, event())) if (runRpcTfTrack_) matchRPCTfTrackToSimTrack(*hRpcTfTrack.product());
}

void 
L1TrackFinderTrackMatcher::matchCSCTfTrackToSimTrack(const L1CSCTrackCollection&)
{}

void 
L1TrackFinderTrackMatcher::matchDTTfTrackToSimTrack(const L1CSCTrackCollection&)
{}

void 
L1TrackFinderTrackMatcher::matchRPCTfTrackToSimTrack(const L1CSCTrackCollection&)
{}

