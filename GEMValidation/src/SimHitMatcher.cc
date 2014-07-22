#include "GEMCode/GEMValidation/src/SimHitMatcher.h"

#include <algorithm>
#include <iomanip>

using namespace std;


SimHitMatcher::SimHitMatcher(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es)
: BaseMatcher(t, v, ps, ev, es)
{
  auto gemSimHit_ = conf().getParameter<edm::ParameterSet>("gemSimHit");
  verboseGEM_ = gemSimHit_.getParameter<int>("verbose");
  gemSimHitInput_ = gemSimHit_.getParameter<edm::InputTag>("input");
  simMuOnlyGEM_ = gemSimHit_.getParameter<bool>("simMuOnly");
  discardEleHitsGEM_ = gemSimHit_.getParameter<bool>("discardEleHits");

  auto cscSimHit_= conf().getParameter<edm::ParameterSet>("cscSimHit");
  verboseCSC_ = cscSimHit_.getParameter<int>("verbose");
  cscSimHitInput_ = cscSimHit_.getParameter<edm::InputTag>("input");
  simMuOnlyCSC_ = cscSimHit_.getParameter<bool>("simMuOnly");
  discardEleHitsCSC_ = cscSimHit_.getParameter<bool>("discardEleHits");

  auto me0SimHit_ = conf().getParameter<edm::ParameterSet>("me0SimHit");
  verboseME0_ = me0SimHit_.getParameter<int>("verbose");
  me0SimHitInput_ = me0SimHit_.getParameter<edm::InputTag>("input");
  simMuOnlyME0_ = me0SimHit_.getParameter<bool>("simMuOnly");
  discardEleHitsME0_ = me0SimHit_.getParameter<bool>("discardEleHits");

  auto rpcSimHit_ = conf().getParameter<edm::ParameterSet>("rpcSimHit");
  verboseRPC_ = rpcSimHit_.getParameter<int>("verbose");
  rpcSimHitInput_ = rpcSimHit_.getParameter<edm::InputTag>("input");
  simMuOnlyRPC_ = rpcSimHit_.getParameter<bool>("simMuOnly");
  discardEleHitsRPC_ = rpcSimHit_.getParameter<bool>("discardEleHits");

  simInputLabel_ = conf().getUntrackedParameter<std::string>("simInputLabel", "g4SimHits");

  setVerbose(conf().getUntrackedParameter<int>("verboseSimHit", 0));

  init();
}


SimHitMatcher::~SimHitMatcher() {}


void 
SimHitMatcher::init()
{
  edm::Handle<edm::PSimHitContainer> csc_hits;
  edm::Handle<edm::PSimHitContainer> gem_hits;
  edm::Handle<edm::PSimHitContainer> rpc_hits;
  edm::Handle<edm::PSimHitContainer> me0_hits;
  edm::Handle<edm::SimTrackContainer> sim_tracks;
  edm::Handle<edm::SimVertexContainer> sim_vertices;

  event().getByLabel(simInputLabel_, sim_tracks);
  event().getByLabel(simInputLabel_, sim_vertices);
  event().getByLabel(cscSimHitInput_, csc_hits);
  event().getByLabel(gemSimHitInput_, gem_hits);
  event().getByLabel(rpcSimHitInput_, rpc_hits);
  event().getByLabel(me0SimHitInput_, me0_hits);

  // fill trkId2Index associoation:
  int no = 0;
  trkid_to_index_.clear();
  for (auto& t: *sim_tracks.product())
  {
    trkid_to_index_[t.trackId()] = no;
    no++;
  }
  vector<unsigned> track_ids = getIdsOfSimTrackShower(trk().trackId(), *sim_tracks.product(), *sim_vertices.product());

  // select CSC simhits
  edm::PSimHitContainer csc_hits_select;
  for (auto& h: *csc_hits.product())
  {
    CSCDetId id(h.detUnitId());
    if ( useCSCChamberType(id.iChamberType()) )  csc_hits_select.push_back(h);
  }

  matchSimHitsToSimTrack(track_ids, csc_hits_select, *gem_hits.product(),
                         *rpc_hits.product(),*me0_hits.product());

  if (verboseCSC_)
  {
    cout<<"sh tn ntids "<<no<<" "<<track_ids.size()<<" "<<csc_hits_select.size()<<endl;
    cout<<"detids "<<detIdsGEM().size()<<" "<<detIdsCSC().size()<<endl;

    auto gem_det_ids = detIdsGEM();
    for (auto id: gem_det_ids)
    {
      //auto& gem_simhits = hitsInDetId(id);
      auto gem_simhits = hitsInDetId(id);
      auto gem_simhits_gp = simHitsMeanPosition(gem_simhits);
      auto strips = hitStripsInDetId(id);
      cout<<"detid "<<GEMDetId(id)<<": "<<gem_simhits.size()<<" "<<gem_simhits_gp.phi()<<" "<< gem_detid_to_hits_[id].size()<<endl;
      cout<<"nstrp "<<strips.size()<<endl;
      cout<<"strps : "; std::copy(strips.begin(), strips.end(), ostream_iterator<int>(cout, " ")); cout<<endl;
    }
  }

  if (verboseGEM_) 
  {
    auto gem_ch_ids = chamberIdsGEM();
    for (auto id: gem_ch_ids)
    {
      auto& gem_simhits = hitsInChamber(id);
      auto gem_simhits_gp = simHitsMeanPosition(gem_simhits);
      cout<<"cchid "<<GEMDetId(id)<<": "<<gem_simhits.size()<<" "<<gem_simhits_gp.phi()<<" "<< gem_chamber_to_hits_[id].size()<<endl;
    }
    auto gem_sch_ids = superChamberIdsGEM();
    for (auto id: gem_sch_ids)
    {
      auto& gem_simhits = hitsInSuperChamber(id);
      auto gem_simhits_gp = simHitsMeanPosition(gem_simhits);
      cout<<"schid "<<GEMDetId(id)<<": "<<nCoincidencePadsWithHits() <<" | "<<gem_simhits.size()<<" "<<gem_simhits_gp.phi()<<" "<< gem_superchamber_to_hits_[id].size()<<endl;
    }
  }
}


std::vector<unsigned int>
SimHitMatcher::getIdsOfSimTrackShower(unsigned int initial_trk_id,
    const edm::SimTrackContainer & sim_tracks, const edm::SimVertexContainer & sim_vertices)
{
  vector<unsigned int> result;
  result.push_back(initial_trk_id);

  if (! (simMuOnlyGEM_ || simMuOnlyCSC_) ) return result;

  for (auto& t: sim_tracks)
  {
    SimTrack last_trk = t;
    bool is_child = 0;
    while (1)
    {
      if ( last_trk.noVertex() ) break;
      if ( sim_vertices[last_trk.vertIndex()].noParent() ) break;
      
      unsigned parentId = sim_vertices[last_trk.vertIndex()].parentIndex();
      if ( parentId == initial_trk_id )
      {
        is_child = 1;
        break;
      }
      
      auto association = trkid_to_index_.find( parentId );
      if ( association == trkid_to_index_.end() ) break;

      last_trk = sim_tracks[ association->second ];
    }
    if (is_child)
    {
      result.push_back(t.trackId());
    }
  }
  return result;
}


void
SimHitMatcher::matchSimHitsToSimTrack(std::vector<unsigned int> track_ids,
                                      const edm::PSimHitContainer& csc_hits, 
                                      const edm::PSimHitContainer& gem_hits,
                                      const edm::PSimHitContainer& rpc_hits, 
                                      const edm::PSimHitContainer& me0_hits)
{
  for (auto& track_id: track_ids)
  {
    for (auto& h: csc_hits)
    {
      if (h.trackId() != track_id) continue;
      int pdgid = h.particleType();
      if (simMuOnlyCSC_ && std::abs(pdgid) != 13) continue;
      // discard electron hits in the CSC chambers
      if (discardEleHitsCSC_ && pdgid == 11) continue;

      csc_detid_to_hits_[ h.detUnitId() ].push_back(h);
      csc_hits_.push_back(h);
      CSCDetId layer_id( h.detUnitId() );
      csc_chamber_to_hits_[ layer_id.chamberId().rawId() ].push_back(h);
    }
    for (auto& h: gem_hits)
    {
      if (h.trackId() != track_id) continue;
      int pdgid = h.particleType();
      if (simMuOnlyGEM_ && std::abs(pdgid) != 13) continue;
      // discard electron hits in the GEM chambers
      if (discardEleHitsGEM_ && pdgid == 11) continue;

      gem_detid_to_hits_[ h.detUnitId() ].push_back(h);
      gem_hits_.push_back(h);
      GEMDetId p_id( h.detUnitId() );
      gem_chamber_to_hits_[ p_id.chamberId().rawId() ].push_back(h);
      GEMDetId superch_id(p_id.region(), p_id.ring(), p_id.station(), 1, p_id.chamber(), 0);
      gem_superchamber_to_hits_[ superch_id() ].push_back(h);
    }
    for (auto& h: rpc_hits)
    {
      if (h.trackId() != track_id) continue;
      int pdgid = h.particleType();
      if (simMuOnlyRPC_ && std::abs(pdgid) != 13) continue;
      // discard electron hits in the RPC chambers
      if (discardEleHitsRPC_ && pdgid == 11) continue;

      rpc_detid_to_hits_[ h.detUnitId() ].push_back(h);
      rpc_hits_.push_back(h);
      RPCDetId layer_id( h.detUnitId() );
      rpc_chamber_to_hits_[ layer_id.chamberId().rawId() ].push_back(h);
    }
    for (auto& h: me0_hits)
    {
      if (h.trackId() != track_id) continue;
      int pdgid = h.particleType();
      if (simMuOnlyME0_ && std::abs(pdgid) != 13) continue;
      // discard electron hits in the ME0 chambers
      if (discardEleHitsME0_ && pdgid == 11) continue;

      me0_detid_to_hits_[ h.detUnitId() ].push_back(h);
      me0_hits_.push_back(h);
      ME0DetId layer_id( h.detUnitId() );
      me0_chamber_to_hits_[ layer_id.chamberId().rawId() ].push_back(h);
    }
  }

  // find pads with hits
  auto detids = detIdsGEM();
  // find 2-layer coincidence pads with hits
  for (auto d: detids)
  {
    GEMDetId id(d);
    auto hits = hitsInDetId(d);
    auto roll = gemGeometry_->etaPartition(id);
    //int max_npads = roll->npads();
    set<int> pads;
    for (auto& h: hits)
    {
      LocalPoint lp = h.entryPoint();
      pads.insert( 1 + static_cast<int>(roll->padTopology().channel(lp)) );
    }
    gem_detids_to_pads_[d] = pads;
  }

  // find 2-layer coincidence pads with hits
  for (auto d: detids)
  {
    GEMDetId id1(d);
    if (id1.layer() != 1) continue;
    GEMDetId id2(id1.region(), id1.ring(), id1.station(), 2, id1.chamber(), id1.roll());
    // does layer 2 has simhits?
    if (detids.find(id2()) == detids.end()) continue;

    // find pads with hits in layer1
    auto hits1 = hitsInDetId(d);
    auto roll1 = gemGeometry_->etaPartition(id1);
    set<int> pads1;
    for (auto& h: hits1)
    {
      LocalPoint lp = h.entryPoint();
      pads1.insert( 1 + static_cast<int>(roll1->padTopology().channel(lp)) );
    }

    // find pads with hits in layer2
    auto hits2 = hitsInDetId(id2());
    auto roll2 = gemGeometry_->etaPartition(id2);
    set<int> pads2;
    for (auto& h: hits2)
    {
      LocalPoint lp = h.entryPoint();
      pads2.insert( 1 + static_cast<int>(roll2->padTopology().channel(lp)) );
    }

    set<int> copads;
    std::set_intersection(pads1.begin(), pads1.end(), pads2.begin(), pads2.end(), std::inserter(copads, copads.begin()));
    if (copads.empty()) continue;
    gem_detids_to_copads_[d] = copads;
  }
}


std::set<unsigned int> 
SimHitMatcher::detIdsGEM() const
{
  std::set<unsigned int> result;
  for (auto& p: gem_detid_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::detIdsRPC() const
{
  std::set<unsigned int> result;
  for (auto& p: rpc_detid_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::detIdsME0() const
{
  std::set<unsigned int> result;
  for (auto& p: me0_detid_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::detIdsCSC(int csc_type) const
{
  std::set<unsigned int> result;
  for (auto& p: csc_detid_to_hits_)
  {
    auto id = p.first;
    if (csc_type > 0)
    {
      CSCDetId detId(id);
      if (detId.iChamberType() != csc_type) continue;
    }
    result.insert(id);
  }
  return result;
}


std::set<unsigned int> 
SimHitMatcher::detIdsGEMCoincidences() const
{
  std::set<unsigned int> result;
  for (auto& p: gem_detids_to_copads_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::detIdsME0Coincidences(int min_n_layers) const
{
  std::set<unsigned int> result;

  //   int result = 0;
  //   auto chamber_ids = chamberIdsCSC();
  //   for (auto id: chamber_ids)
  //   {
  //     if (nLayersWithHitsInSuperChamber(id) >= min_n_layers) result += 1;
  //   }
  //   return result;

  //   for (auto& p: me0_detids_to_copads_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::chamberIdsGEM() const
{
  std::set<unsigned int> result;
  for (auto& p: gem_chamber_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::chamberIdsRPC() const
{
  std::set<unsigned int> result;
  for (auto& p: rpc_chamber_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::chamberIdsME0() const
{
  std::set<unsigned int> result;
  for (auto& p: me0_chamber_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
SimHitMatcher::chamberIdsCSC(int csc_type) const
{
  std::set<unsigned int> result;
  for (auto& p: csc_chamber_to_hits_)
  {
    auto id = p.first;
    if (csc_type > 0)
    {
      CSCDetId detId(id);
      if (detId.iChamberType() != csc_type) continue;
    }
    result.insert(id);
  }
  return result;
}


std::set<unsigned int>
SimHitMatcher::superChamberIdsGEM() const
{
  std::set<unsigned int> result;
  for (auto& p: gem_superchamber_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
SimHitMatcher::superChamberIdsME0() const
{
  std::set<unsigned int> result;
  for (auto& p: me0_superchamber_to_hits_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
SimHitMatcher::superChamberIdsGEMCoincidences() const
{
  std::set<unsigned int> result;
  for (auto& p: gem_detids_to_copads_)
  {
    GEMDetId id(p.first);
    result.insert(id.chamberId().rawId());
  }
  return result;
}


const edm::PSimHitContainer&
SimHitMatcher::hitsInDetId(unsigned int detid) const
{
  if (is_gem(detid))
  {
    if (gem_detid_to_hits_.find(detid) == gem_detid_to_hits_.end()) return no_hits_;
    return gem_detid_to_hits_.at(detid);
  }
  if (is_me0(detid))
  {
    if (me0_detid_to_hits_.find(detid) == me0_detid_to_hits_.end()) return no_hits_;
    return me0_detid_to_hits_.at(detid);
  }
  if (is_csc(detid))
  {
    if (csc_detid_to_hits_.find(detid) == csc_detid_to_hits_.end()) return no_hits_;
    return csc_detid_to_hits_.at(detid);
  }
  if (is_rpc(detid))
  {
    if (rpc_detid_to_hits_.find(detid) == rpc_detid_to_hits_.end()) return no_hits_;
    return rpc_detid_to_hits_.at(detid);
  }
  return no_hits_;
}


const edm::PSimHitContainer&
SimHitMatcher::hitsInChamber(unsigned int detid) const
{
  if (is_gem(detid)) // make sure we use chamber id
  {
    GEMDetId id(detid);
    if (gem_chamber_to_hits_.find(id.chamberId().rawId()) == gem_chamber_to_hits_.end()) return no_hits_;
    return gem_chamber_to_hits_.at(id.chamberId().rawId());
  }
  if (is_me0(detid)) // make sure we use chamber id
  {
    ME0DetId id(detid);
    if (me0_chamber_to_hits_.find(id.chamberId().rawId()) == me0_chamber_to_hits_.end()) return no_hits_;
    return me0_chamber_to_hits_.at(id.chamberId().rawId());
  }
  if (is_csc(detid))
  {
    CSCDetId id(detid);
    if (csc_chamber_to_hits_.find(id.chamberId().rawId()) == gem_chamber_to_hits_.end()) return no_hits_;
    return csc_chamber_to_hits_.at(id.chamberId().rawId());
  }
  if (is_rpc(detid))
  {
    RPCDetId id(detid);
    if (rpc_chamber_to_hits_.find(id.chamberId().rawId()) == gem_chamber_to_hits_.end()) return no_hits_;
    return rpc_chamber_to_hits_.at(id.chamberId().rawId());
  }
  return no_hits_;
}


const edm::PSimHitContainer&
SimHitMatcher::hitsInSuperChamber(unsigned int detid) const
{
  if (is_gem(detid))
  {
    GEMDetId id(detid);
    if (gem_superchamber_to_hits_.find(id.chamberId().rawId()) == gem_superchamber_to_hits_.end()) return no_hits_;
    return gem_superchamber_to_hits_.at(id.chamberId().rawId());
  }
  if (is_csc(detid)) return hitsInChamber(detid);

  return no_hits_;
}


int
SimHitMatcher::nLayersWithHitsInSuperChamber(unsigned int detid) const
{
  set<int> layers_with_hits;
  auto hits = hitsInSuperChamber(detid);
  for (auto& h: hits)
  {
    if (is_gem(detid))
    {
      GEMDetId idd(h.detUnitId());
      layers_with_hits.insert(idd.layer());
    }
    if (is_me0(detid))
    {
      ME0DetId idd(h.detUnitId());
      layers_with_hits.insert(idd.layer());
    }
    if (is_csc(detid))
    {
      CSCDetId idd(h.detUnitId());
      layers_with_hits.insert(idd.layer());
    }
    if (is_rpc(detid))
    {
      RPCDetId idd(h.detUnitId());
      layers_with_hits.insert(idd.layer());
    }
  }
  return layers_with_hits.size();
}


GlobalPoint
SimHitMatcher::simHitsMeanPosition(const edm::PSimHitContainer& sim_hits) const
{
  if (sim_hits.empty()) return GlobalPoint(); // point "zero"

  float sumx, sumy, sumz;
  sumx = sumy = sumz = 0.f;
  size_t n = 0;
  for (auto& h: sim_hits)
  {
    LocalPoint lp = h.entryPoint();
    GlobalPoint gp;
    if ( is_gem(h.detUnitId()) )
    {
      gp = gemGeometry_->idToDet(h.detUnitId())->surface().toGlobal(lp);
  //    std::cout<<"gem simhit eta:"<<gp.eta()<<"  phi:"<<gp.phi()<<std::endl;
    }
    else if ( is_me0(h.detUnitId()) )
    {
      gp = me0Geometry_->idToDet(h.detUnitId())->surface().toGlobal(lp);
    }
    else if (is_csc(h.detUnitId()))
    {
      gp = cscGeometry_->idToDet(h.detUnitId())->surface().toGlobal(lp);
    }
    else if (is_rpc(h.detUnitId()))
    {
      gp = rpcGeometry_->idToDet(h.detUnitId())->surface().toGlobal(lp);
    }
    else continue;
    sumx += gp.x();
    sumy += gp.y();
    sumz += gp.z();
    ++n;
  }
  if (n == 0) return GlobalPoint();
  return GlobalPoint(sumx/n, sumy/n, sumz/n);
}


float 
SimHitMatcher::simHitsMeanStrip(const edm::PSimHitContainer& sim_hits) const
{
  if (sim_hits.empty()) return -1.f;

  float sums = 0.f;
  size_t n = 0;
  for (auto& h: sim_hits)
  {
    LocalPoint lp = h.entryPoint();
    float s;
    auto d = h.detUnitId();
    if ( is_gem(d) )
    {
      s = gemGeometry_->etaPartition(d)->strip(lp);
    }
    if ( is_me0(d) )
    {
      s = me0Geometry_->etaPartition(d)->strip(lp);
    }
    else if (is_csc(d))
    {
      s = cscGeometry_->layer(d)->geometry()->strip(lp);
      // convert to half-strip:
      s *= 2.;
    }
    else if (is_rpc(d))
    {
      s = rpcGeometry_->roll(d)->strip(lp);
    }
    else continue;
    sums += s;
    ++n;
  }
  if (n == 0) return -1.f;
  return sums/n;
}


float 
SimHitMatcher::simHitsMeanWG(const edm::PSimHitContainer& sim_hits) const
{
  if (sim_hits.empty()) return -1.f;

  float sums = 0.f;
  size_t n = 0;
  for (auto& h: sim_hits)
  {
    LocalPoint lp = h.entryPoint();
    float s;
    auto d = h.detUnitId();
    if (is_csc(d))
    {
      // find nearest wire
      int nearestWire(cscGeometry_->layer(d)->geometry()->nearestWire(lp));
      // then find the corresponding wire group
      s = cscGeometry_->layer(d)->geometry()->wireGroup(nearestWire);
    }
    else continue;
    sums += s;
    ++n;
  }
  if (n == 0) return -1.f;
  return sums/n;
}



std::set<int> 
SimHitMatcher::hitStripsInDetId(unsigned int detid, int margin_n_strips) const
{
  set<int> result;
  auto simhits = hitsInDetId(detid);
  if ( is_gem(detid) )
  {
    GEMDetId id(detid);
    int max_nstrips = gemGeometry_->etaPartition(id)->nstrips();
    for (auto& h: simhits)
    {
      LocalPoint lp = h.entryPoint();
      int central_strip = 1 + static_cast<int>(gemGeometry_->etaPartition(id)->topology().channel(lp));
      int smin = central_strip - margin_n_strips;
      smin = (smin > 0) ? smin : 1;
      int smax = central_strip + margin_n_strips;
      smax = (smax <= max_nstrips) ? smax : max_nstrips;
      for (int ss = smin; ss <= smax; ++ss) result.insert(ss);
    }
  }
  else if ( is_me0(detid) )
  {
    ME0DetId id(detid);
    int max_nstrips = me0Geometry_->etaPartition(id)->nstrips();
    for (auto& h: simhits)
    {
      LocalPoint lp = h.entryPoint();
      int central_strip = 1 + static_cast<int>(me0Geometry_->etaPartition(id)->topology().channel(lp));
      int smin = central_strip - margin_n_strips;
      smin = (smin > 0) ? smin : 1;
      int smax = central_strip + margin_n_strips;
      smax = (smax <= max_nstrips) ? smax : max_nstrips;
      for (int ss = smin; ss <= smax; ++ss) result.insert(ss);
    }
  }
  else if ( is_csc(detid) )
  {
    CSCDetId id(detid);
    int max_nstrips = cscGeometry_->layer(id)->geometry()->numberOfStrips();
    for (auto& h: simhits)
    {
      LocalPoint lp = h.entryPoint();
      int central_strip = cscGeometry_->layer(id)->geometry()->nearestStrip(lp);
      int smin = central_strip - margin_n_strips;
      smin = (smin > 0) ? smin : 1;
      int smax = central_strip + margin_n_strips;
      smax = (smax <= max_nstrips) ? smax : max_nstrips;
      for (int ss = smin; ss <= smax; ++ss) result.insert(ss);
    }
  }
  else if ( is_rpc(detid) )
  {
    RPCDetId id(detid); 
    int max_nstrips = rpcGeometry_->roll(id)->nstrips();
    for (auto& h: simhits)
    {
      LocalPoint lp = h.entryPoint();
      int central_strip = 1 + static_cast<int>(rpcGeometry_->roll(id)->topology().channel(lp));
    //  int central_strip2 = 1 + static_cast<int>(rpcGeometry_->roll(id)->strip(lp));
    //  std::cout <<"strip from topology"<< central_strip <<" strip from roll" << central_strip2 <<std::endl; 
      int smin = central_strip - margin_n_strips;
      smin = (smin > 0) ? smin : 1;
      int smax = central_strip + margin_n_strips;
      smax = (smax <= max_nstrips) ? smax : max_nstrips;
      for (int ss = smin; ss <= smax; ++ss) result.insert(ss);
    }
  
  }
  return result;
}


std::set<int> 
SimHitMatcher::hitWiregroupsInDetId(unsigned int detid, int margin_n_wg) const
{
  set<int> result;
  if ( !is_csc(detid) ) return result;

  auto simhits = hitsInDetId(detid);
  CSCDetId id(detid);
  int max_n_wg = cscGeometry_->layer(id)->geometry()->numberOfWireGroups();
  for (auto& h: simhits)
  {
    LocalPoint lp = h.entryPoint();
    auto layer_geo = cscGeometry_->layer(id)->geometry();
    int central_wg = layer_geo->wireGroup(layer_geo->nearestWire(lp));
    int wg_min = central_wg - margin_n_wg;
    wg_min = (wg_min > 0) ? wg_min : 1;
    int wg_max = central_wg + margin_n_wg;
    wg_max = (wg_max <= max_n_wg) ? wg_max : max_n_wg;
    for (int wg = wg_min; wg <= wg_max; ++wg) result.insert(wg);

  }
  return result;
}


std::set<int> 
SimHitMatcher::hitPadsInDetId(unsigned int detid) const
{
  set<int> none;
  if (gem_detids_to_pads_.find(detid) == gem_detids_to_pads_.end()) return none;
  return gem_detids_to_pads_.at(detid);
}


std::set<int>
SimHitMatcher::hitCoPadsInDetId(unsigned int detid) const
{
  set<int> none;
  if (gem_detids_to_copads_.find(detid) == gem_detids_to_copads_.end()) return none;
  return gem_detids_to_copads_.at(detid);
}


std::set<int> 
SimHitMatcher::hitPartitions() const
{
  std::set<int> result;

  auto detids = detIdsGEM();
  for (auto id: detids)
  {
    GEMDetId idd(id);
    result.insert( idd.roll() );
  }
  return result;
}


int
SimHitMatcher::nPadsWithHits() const
{
  int result = 0;
  auto pad_ids = detIdsGEM();
  for (auto id: pad_ids)
  {
    result += hitPadsInDetId(id).size();
  }
  return result;
}


int
SimHitMatcher::nCoincidencePadsWithHits() const
{
  int result = 0;
  auto copad_ids = detIdsGEMCoincidences();
  for (auto id: copad_ids)
  {
    result += hitCoPadsInDetId(id).size();
  }
  return result;
}


int
SimHitMatcher::nCoincidenceCSCChambers(int min_n_layers) const
{
  int result = 0;
  auto chamber_ids = chamberIdsCSC();
  for (auto id: chamber_ids)
  {
    if (nLayersWithHitsInSuperChamber(id) >= min_n_layers) result += 1;
  }
  return result;
}


