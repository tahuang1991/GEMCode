#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"

SimTrackMatchManager::SimTrackMatchManager(const SimTrack& t, 
                                           const SimVertex& v,
                                           const edm::ParameterSet& ps, 
                                           const edm::Event& ev, 
                                           const edm::EventSetup& es, 
                                           edm::ConsumesCollector && iC)
: 
 //genMuons_(t, v, ps, ev, es)
  simhits_(t, v, ps, ev, es, iC)
, gem_digis_(simhits_, iC)
, gem_rechits_(simhits_, iC)
, me0_digis_(simhits_, iC)
// , me0_rechits_(simhits_)
, rpc_digis_(simhits_, iC)
, rpc_rechits_(simhits_, iC)
, csc_digis_(simhits_, iC)
, csc_stubs_(simhits_, csc_digis_, gem_digis_, rpc_digis_, iC)
  , csc_rechits_(simhits_, iC)
, dt_digis_(simhits_, iC)
, dt_stubs_(simhits_, iC)
, dt_rechits_(simhits_, iC)
  //, l1_tracks_(csc_stubs_, dt_digis_, rpc_digis_)
, l1_tf_tracks_(simhits_)
, l1_tf_cands_(simhits_)
, l1_gmt_cands_(simhits_)
, hlt_tracks_(csc_rechits_, dt_rechits_, rpc_rechits_, gem_rechits_)
{
  //  std::cout <<"Constructing new SimTrackMatchManager" << std::endl;
}

SimTrackMatchManager::~SimTrackMatchManager() 
{
  //  std::cout <<"Removing SimTrackMatchManager" << std::endl;
}
