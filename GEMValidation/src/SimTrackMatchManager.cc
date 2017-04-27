#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"

SimTrackMatchManager::SimTrackMatchManager(const SimTrack& t, 
                                           const SimVertex& v,
                                           const edm::ParameterSet& ps, 
                                           const edm::Event& ev, 
                                           const edm::EventSetup& es, 
                                           edm::EDGetTokenT<reco::GenParticleCollection>& genParticleInput_,
                                           edm::EDGetTokenT<edm::SimVertexContainer>& simVertexInput_,
                                           edm::EDGetTokenT<edm::SimTrackContainer>& simTrackInput_,
                                           edm::EDGetTokenT<edm::PSimHitContainer>& gemSimHitInput_,
                                           edm::EDGetTokenT<edm::PSimHitContainer>& cscSimHitInput_,
                                           edm::EDGetTokenT<edm::PSimHitContainer>& rpcSimHitInput_,
                                           edm::EDGetTokenT<edm::PSimHitContainer>& me0SimHitInput_,
                                           edm::EDGetTokenT<edm::PSimHitContainer>& dtSimHitInput_,
                                           edm::EDGetTokenT<GEMDigiCollection>& gemDigiInput_,
                                           edm::EDGetTokenT<GEMPadDigiCollection>& gemPadDigiInput_,
                                           edm::EDGetTokenT<GEMCoPadDigiCollection>& gemCoPadDigiInput_,
                                           edm::EDGetTokenT<GEMRecHitCollection>& gemRecHitInput_,
                                           edm::EDGetTokenT<ME0DigiPreRecoCollection>& me0DigiInput_,
                                           edm::EDGetTokenT<CSCComparatorDigiCollection>& cscComparatorDigiInput_,
                                           edm::EDGetTokenT<CSCWireDigiCollection>& cscWireDigiInput_,
                                           edm::EDGetTokenT<CSCCLCTDigiCollection>& clctInputs_,
                                           edm::EDGetTokenT<CSCALCTDigiCollection>& alctInputs_,
                                           edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection>& lctInputs_,
                                           edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection>& mplctInputs_,
                                           edm::EDGetTokenT<CSCRecHit2DCollection>& cscRecHit2DInput_,
                                           edm::EDGetTokenT<CSCSegmentCollection>& cscSegmentInput_,
                                           edm::EDGetTokenT<DTDigiCollection>& dtDigiInput_,
                                           edm::EDGetTokenT<DTLocalTriggerCollection>& dtStubInput_,
                                           edm::EDGetTokenT<DTRecHitCollection>& dtRecHit1DPairInput_,
                                           edm::EDGetTokenT<DTRecSegment2DCollection>& dtRecSegment2DInput_,
                                           edm::EDGetTokenT<DTRecSegment4DCollection>& dtRecSegment4DInput_,
                                           edm::EDGetTokenT<RPCDigiCollection>& rpcDigiInput_,
                                           edm::EDGetTokenT<RPCRecHitCollection>& rpcRecHitInput_,
                                           //edm::EDGetTokenT<L1CSCTrackCollection>& cscTfTrackInputLabel_, 
                                           //edm::EDGetTokenT<L1MuRegionalCandCollection>& cscTfCandInputLabel_, 
					   edm::EDGetTokenT<l1t::EMTFTrackCollection>& emtfTrackInputLabel_,
					   edm::EDGetTokenT< BXVector<l1t::RegionalMuonCand> >& gmtInputLabel_,
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& dtTfCandInputLabel_, 
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& rpcfTfCandInputLabel_, 
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& rpcbTfCandInputLabel_, 
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& gmtRegCandCSCInputLabel_,
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& gmtRegCandDTInputLabel_,
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& gmtRegCandRPCfInputLabel_,
                                           edm::EDGetTokenT<L1MuRegionalCandCollection>& gmtRegCandRPCbInputLabel_,
                                           edm::EDGetTokenT<L1MuGMTCandCollection>& gmtCandInputLabel_,
                                           edm::EDGetTokenT<l1extra::L1MuonParticleCollection>& l1ExtraMuonInputLabel_,
                                           edm::EDGetTokenT<reco::TrackExtraCollection>& recoTrackExtraInputLabel_,
                                           edm::EDGetTokenT<reco::TrackCollection>& recoTrackInputLabel_,
                                           edm::EDGetTokenT<reco::RecoChargedCandidateCollection>& recoChargedCandidateInputLabel_
					   )
: 
  genMuons_(t, v, ps, ev, es, genParticleInput_)
  ,simhits_(t, v, ps, ev, es, simVertexInput_, simTrackInput_, 
           gemSimHitInput_, cscSimHitInput_, 
           rpcSimHitInput_, me0SimHitInput_, dtSimHitInput_)
  , gem_digis_(simhits_, gemDigiInput_, gemPadDigiInput_, gemCoPadDigiInput_)
  , gem_rechits_(simhits_, gemRecHitInput_)
  , me0_digis_(simhits_, me0DigiInput_)
  // , me0_rechits_(simhits_)
  , rpc_digis_(simhits_, rpcDigiInput_)
  , rpc_rechits_(simhits_, rpcRecHitInput_)
  , csc_digis_(simhits_, cscComparatorDigiInput_, cscWireDigiInput_)
  , csc_stubs_(simhits_, 
               csc_digis_, gem_digis_, 
               rpc_digis_, clctInputs_, 
               alctInputs_, lctInputs_, mplctInputs_)
  , csc_rechits_(simhits_, 
                 cscRecHit2DInput_, cscSegmentInput_)
  , dt_digis_(simhits_, dtDigiInput_)
  , dt_stubs_(simhits_, dtStubInput_)
  , dt_rechits_(simhits_, 
                dtRecHit1DPairInput_, dtRecSegment2DInput_, 
                dtRecSegment4DInput_)
  , l1_tracks_(csc_stubs_, emtfTrackInputLabel_, gmtInputLabel_)
  //, l1_tracks_(csc_stubs_, dt_digis_, rpc_digis_, 
  //	  	cscTfTrackInputLabel_, cscTfCandInputLabel_)
  , l1_tf_tracks_(simhits_)//, cscTfTrackInputLabel_)
  , l1_tf_cands_(simhits_, 
  //               cscTfCandInputLabel_, 
                 dtTfCandInputLabel_, 
                 rpcfTfCandInputLabel_, rpcbTfCandInputLabel_)
  , l1_gmt_cands_(simhits_, 
                  gmtRegCandCSCInputLabel_, gmtRegCandDTInputLabel_, 
                  gmtRegCandRPCfInputLabel_, gmtRegCandRPCbInputLabel_, 
                  gmtCandInputLabel_, l1ExtraMuonInputLabel_)
  , hlt_tracks_(csc_rechits_, dt_rechits_, rpc_rechits_, gem_rechits_, 
                recoTrackExtraInputLabel_, recoTrackInputLabel_, recoChargedCandidateInputLabel_)
		
{
  //  std::cout <<"Constructing new SimTrackMatchManager" << std::endl;
}

SimTrackMatchManager::~SimTrackMatchManager() 
{
  //  std::cout <<"Removing SimTrackMatchManager" << std::endl;
}
