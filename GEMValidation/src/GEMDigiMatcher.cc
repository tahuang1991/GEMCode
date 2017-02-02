#include "GEMCode/GEMValidation/interface/GEMDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;
using namespace matching;

GEMDigiMatcher::GEMDigiMatcher(SimHitMatcher& sh, 
                               edm::EDGetTokenT<GEMDigiCollection> &gemDigiInput_, 
                               edm::EDGetTokenT<GEMPadDigiCollection> &gemPadDigiInput_, 
                               edm::EDGetTokenT<GEMCoPadDigiCollection> &gemCoPadDigiInput_)
: DigiMatcher(sh)
{
  auto gemDigi_= conf().getParameter<edm::ParameterSet>("gemStripDigi");
  minBXGEMDigi_ = gemDigi_.getParameter<int>("minBX");
  maxBXGEMDigi_ = gemDigi_.getParameter<int>("maxBX");
  matchDeltaStrip_ = gemDigi_.getParameter<int>("matchDeltaStrip");
  verboseDigi_ = gemDigi_.getParameter<int>("verbose");
  runGEMDigi_ = gemDigi_.getParameter<bool>("run");

  auto gemPad_= conf().getParameter<edm::ParameterSet>("gemPadDigi");
  minBXGEMPad_ = gemPad_.getParameter<int>("minBX");
  maxBXGEMPad_ = gemPad_.getParameter<int>("maxBX");
  verbosePad_ = gemPad_.getParameter<int>("verbose");
  runGEMPad_ = gemPad_.getParameter<bool>("run");

  auto gemCoPad_= conf().getParameter<edm::ParameterSet>("gemCoPadDigi");
  minBXGEMCoPad_ = gemCoPad_.getParameter<int>("minBX");
  maxBXGEMCoPad_ = gemCoPad_.getParameter<int>("maxBX");
  verboseCoPad_ = gemCoPad_.getParameter<int>("verbose");
  runGEMCoPad_ = gemCoPad_.getParameter<bool>("run");
  if (hasGEMGeometry_) {
    edm::Handle<GEMDigiCollection> gem_digis;
    if (gemvalidation::getByToken(gemDigiInput_, gem_digis, event())) if (runGEMDigi_) matchDigisToSimTrack(*gem_digis.product());
    
    edm::Handle<GEMPadDigiCollection> gem_pads;
    if (gemvalidation::getByToken(gemPadDigiInput_, gem_pads, event())) if (runGEMPad_) matchPadsToSimTrack(*gem_pads.product());
    
    edm::Handle<GEMCoPadDigiCollection> gem_co_pads;
    if (gemvalidation::getByToken(gemCoPadDigiInput_, gem_co_pads, event())) if (runGEMCoPad_) matchCoPadsToSimTrack(*gem_co_pads.product());
  }
}

GEMDigiMatcher::~GEMDigiMatcher() {}


void
GEMDigiMatcher::matchDigisToSimTrack(const GEMDigiCollection& digis)
{
  if (verboseDigi_) cout << "Matching simtrack to GEM digis" << endl;
  auto det_ids = simhit_matcher_->detIdsGEM();
  for (auto id: det_ids)
  {
    GEMDetId p_id(id);
    GEMDetId superch_id(p_id.region(), p_id.ring(), p_id.station(), 0, p_id.chamber(), 0);
    auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    if (verboseDigi_)
    {
      cout<<"hit_strips_fat ";
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto digis_in_det = digis.get(GEMDetId(id));

    for (auto d = digis_in_det.first; d != digis_in_det.second; ++d)
    {
      if (verboseDigi_) cout<<"GEMDigi "<<p_id<<" "<<*d<<endl;
      // check that the digi is within BX range
      if (d->bx() < minBXGEMDigi_ || d->bx() > maxBXGEMDigi_) continue;
      // check that it matches a strip that was hit by SimHits from our track
      if (hit_strips.find(d->strip()) == hit_strips.end()) continue;
      if (verboseDigi_) cout<<"...was matched!"<<endl;

      auto mydigi = make_digi(id, d->strip(), d->bx(), GEM_STRIP);
      detid_to_digis_[id].push_back(mydigi);
      chamber_to_digis_[ p_id.chamberId().rawId() ].push_back(mydigi);
      superchamber_to_digis_[ superch_id() ].push_back(mydigi);

      //std::cout <<" strip "<< d->strip()<<" 2-strip pad "<<(d->strip()+1)/2 << " bx "<< d->bx() << std::endl;
      detid_to_gemdigis_[id].push_back(*d);
      GEMPadDigi pad = GEMPadDigi((d->strip()+1)/2, d->bx());
      if (std::find(detid_to_gempads_2strip_[id].begin(), detid_to_gempads_2strip_[id].end(), pad) == detid_to_gempads_2strip_[id].end())
      	detid_to_gempads_2strip_[id].push_back(pad);
      chamber_to_gemdigis_[ p_id.chamberId().rawId() ].push_back(*d);
      superchamber_to_gemdigis_[ superch_id() ].push_back(*d);
      //int pad_num = 1 + static_cast<int>( roll->padOfStrip(d->strip()) ); // d->strip() is int
      //digi_map[ make_pair(pad_num, d->bx()) ].push_back( d->strip() );
    }
  }
}


void
GEMDigiMatcher::matchPadsToSimTrack(const GEMPadDigiCollection& pads)
{
  auto det_ids = simhit_matcher_->detIdsGEM();
  for (auto id: det_ids)
  {
    GEMDetId p_id(id);
    GEMDetId superch_id(p_id.region(), p_id.ring(), p_id.station(), 0, p_id.chamber(), 0);

    auto hit_pads = simhit_matcher_->hitPadsInDetId(id);
    auto pads_in_det = pads.get(p_id);

    if (verbosePad_)
    {
      cout<<"checkpads "<<hit_pads.size()<<" "<<std::distance(pads_in_det.first, pads_in_det.second)<<" hit_pads: ";
      copy(hit_pads.begin(), hit_pads.end(), ostream_iterator<int>(cout," "));
      cout<<endl;
    }

    for (auto pad = pads_in_det.first; pad != pads_in_det.second; ++pad)
    {
      if (verbosePad_) cout<<"chp "<<*pad<<endl;
      // check that the pad BX is within the range
      if (pad->bx() < minBXGEMPad_ || pad->bx() > maxBXGEMPad_) continue;
      if (verbosePad_) cout<<"chp1"<<endl;
      // check that it matches a pad that was hit by SimHits from our track
      if (hit_pads.find(pad->pad()) == hit_pads.end()) continue;
      if (verbosePad_) cout<<"chp2"<<endl;
      auto mydigi = make_digi(id, pad->pad(), pad->bx(), GEM_PAD);
      detid_to_pads_[id].push_back(mydigi);
      chamber_to_pads_[ p_id.chamberId().rawId() ].push_back(mydigi);
      superchamber_to_pads_[ superch_id() ].push_back(mydigi);

      detid_to_gempads_[id].push_back(*pad);
      chamber_to_gempads_[ p_id.chamberId().rawId() ].push_back(*pad);
      superchamber_to_gempads_[ superch_id() ].push_back(*pad);
    }
  }
}


void
GEMDigiMatcher::matchCoPadsToSimTrack(const GEMCoPadDigiCollection& co_pads)
{
  auto det_ids = simhit_matcher_->detIdsGEMCoincidences();
  for (auto id: det_ids)
  {
    GEMDetId p_id(id);
    GEMDetId superch_id(p_id.region(), p_id.ring(), p_id.station(), 0, p_id.chamber(), 0);

    auto hit_co_pads = simhit_matcher_->hitCoPadsInDetId(id);
    auto co_pads_in_det = co_pads.get(p_id);

    if (verboseCoPad_)
    {
      cout<<"checkcopads "<<hit_co_pads.size()<<" "<<std::distance(co_pads_in_det.first, co_pads_in_det.second)<<" hit_pads: ";
      copy(hit_co_pads.begin(), hit_co_pads.end(), ostream_iterator<int>(cout," "));
      cout<<endl;
    }

    for (auto pad = co_pads_in_det.first; pad != co_pads_in_det.second; ++pad)
    {
      if (verbosePad_) cout<<"CoPad: chp "<<*pad<<endl;
      // check that the pad BX is within the range
      if (pad->bx(1) < minBXGEMCoPad_ || pad->bx(1) > maxBXGEMCoPad_) continue;
      if (verboseCoPad_) cout<<"CoPad: chp1 "<<*pad<<endl;
      // check that it matches a pad that was hit by SimHits from our track
      if (hit_co_pads.find(pad->pad(1)) == hit_co_pads.end()) continue;
      if (verboseCoPad_) cout<<"CoPad: chp2 "<<*pad<<endl;
      // check that it matches a coincidence pad that was hit by SimHits from our track

      auto mydigi = make_digi(id, pad->pad(1), pad->bx(1), GEM_COPAD);
      //detid_to_copads_[id].push_back(mydigi);
      superchamber_to_copads_[ superch_id() ].push_back(mydigi);

      superchamber_to_gemcopads_[ superch_id() ].push_back(*pad);
    }
  }
}


std::set<unsigned int>
GEMDigiMatcher::selectDetIds(const Id2DigiContainer &digis, int gem_type) const
{
  std::set<unsigned int> result;
  for (auto& p: digis)
  {
    auto id = p.first;
    if (gem_type > 0)
    {
      GEMDetId detId(id);
      if (gemvalidation::toGEMType(detId.station(),detId.ring()) != gem_type) continue;
    }
    result.insert(p.first);
  }
  return result;
}


std::set<unsigned int>
GEMDigiMatcher::detIdsDigi(int gem_type) const
{
  return selectDetIds(detid_to_digis_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::detIdsPad(int gem_type) const
{
  return selectDetIds(detid_to_pads_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::chamberIdsDigi(int gem_type) const
{
  return selectDetIds(chamber_to_digis_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::chamberIdsPad(int gem_type) const
{
  return selectDetIds(chamber_to_pads_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::superChamberIdsDigi(int gem_type) const
{
  return selectDetIds(superchamber_to_digis_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::superChamberIdsPad(int gem_type) const
{
  return selectDetIds(superchamber_to_pads_, gem_type);
}


std::set<unsigned int>
GEMDigiMatcher::superChamberIdsCoPad(int gem_type) const
{
  return selectDetIds(superchamber_to_copads_, gem_type);
}


const matching::DigiContainer&
GEMDigiMatcher::digisInDetId(unsigned int detid) const
{
  if (detid_to_digis_.find(detid) == detid_to_digis_.end()) return no_digis_;
  return detid_to_digis_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::digisInChamber(unsigned int detid) const
{
  if (chamber_to_digis_.find(detid) == chamber_to_digis_.end()) return no_digis_;
  return chamber_to_digis_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::digisInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_digis_.find(detid) == superchamber_to_digis_.end()) return no_digis_;
  return superchamber_to_digis_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::padsInDetId(unsigned int detid) const
{
  if (detid_to_pads_.find(detid) == detid_to_pads_.end()) return no_digis_;
  return detid_to_pads_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::padsInChamber(unsigned int detid) const
{
  if (chamber_to_pads_.find(detid) == chamber_to_pads_.end()) return no_digis_;
  return chamber_to_pads_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::padsInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_pads_.find(detid) == superchamber_to_pads_.end()) return no_digis_;
  return superchamber_to_pads_.at(detid);
}


const matching::DigiContainer&
GEMDigiMatcher::coPadsInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_copads_.find(detid) == superchamber_to_copads_.end()) return no_digis_;
  return superchamber_to_copads_.at(detid);
}


const GEMDigiContainer&
GEMDigiMatcher::gemDigisInDetId(unsigned int detid) const
{
  if (detid_to_gemdigis_.find(detid) == detid_to_gemdigis_.end()) return no_gem_digis_;
  return detid_to_gemdigis_.at(detid);
}


const GEMDigiContainer&
GEMDigiMatcher::gemDigisInChamber(unsigned int detid) const
{
  if (chamber_to_gemdigis_.find(detid) == chamber_to_gemdigis_.end()) return no_gem_digis_;
  return chamber_to_gemdigis_.at(detid);
}


const GEMDigiContainer&
GEMDigiMatcher::gemDigisInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_gemdigis_.find(detid) == superchamber_to_gemdigis_.end()) return no_gem_digis_;
  return superchamber_to_gemdigis_.at(detid);
}


const GEMPadDigiContainer&
GEMDigiMatcher::gemPadsInDetId(unsigned int detid) const
{
  if (detid_to_gempads_.find(detid) == detid_to_gempads_.end()) return no_gem_pads_;
  return detid_to_gempads_.at(detid);
}


const GEMPadDigiContainer&
GEMDigiMatcher::gemPadsInChamber(unsigned int detid) const
{
  if (chamber_to_gempads_.find(detid) == chamber_to_gempads_.end()) return no_gem_pads_;
  return chamber_to_gempads_.at(detid);
}


const GEMPadDigiContainer&
GEMDigiMatcher::gemPadsInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_gempads_.find(detid) == superchamber_to_gempads_.end()) return no_gem_pads_;
  return superchamber_to_gempads_.at(detid);
}


const GEMCoPadDigiContainer&
GEMDigiMatcher::gemCoPadsInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_gemcopads_.find(detid) == superchamber_to_gemcopads_.end()) return no_gem_copads_;
  return superchamber_to_gemcopads_.at(detid);
}


int
GEMDigiMatcher::nLayersWithDigisInSuperChamber(unsigned int detid) const
{
  set<int> layers;
  auto digis = digisInSuperChamber(detid);
  for (auto& d: digis)
  {
    GEMDetId idd(digi_id(d));
    layers.insert(idd.layer());
  }
  return layers.size();
}


int
GEMDigiMatcher::nLayersWithPadsInSuperChamber(unsigned int detid) const
{
  set<int> layers;
  auto digis = padsInSuperChamber(detid);
  for (auto& d: digis)
  {
    GEMDetId idd(digi_id(d));
    layers.insert(idd.layer());
  }
  return layers.size();
}


int
GEMDigiMatcher::nPads() const
{
  int n = 0;
  auto ids = superChamberIdsPad();
  for (auto id: ids)
  {
    n += padsInSuperChamber(id).size();
  }
  return n;
}


int
GEMDigiMatcher::nCoPads() const
{
  int n = 0;
  auto ids = superChamberIdsCoPad();
  for (auto id: ids)
  {
    n += coPadsInSuperChamber(id).size();
  }
  return n;
}


std::set<int>
GEMDigiMatcher::stripNumbersInDetId(unsigned int detid) const
{
  set<int> result;
  auto digis = digisInDetId(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  return result;
}


std::set<int>
GEMDigiMatcher::padNumbersInDetId(unsigned int detid) const
{
  set<int> result;
  auto digis = padsInDetId(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  return result;
}


std::set<int>
GEMDigiMatcher::partitionNumbers() const
{
  std::set<int> result;

  auto detids = detIdsDigi();
  for (auto id: detids)
  {
    GEMDetId idd(id);
    result.insert( idd.roll() );
  }
  return result;
}


std::set<int>
GEMDigiMatcher::partitionNumbersWithCoPads() const
{
  std::set<int> result;

  auto detids = superChamberIdsCoPad();
  for (auto id: detids)
  {
    GEMDetId idd(id);
    result.insert( idd.roll() );
  }
  return result;
}


int 
GEMDigiMatcher::extrapolateHsfromGEMPad(unsigned int id, int gempad) const
{
  int result = -1 ;
  
  GEMDetId gem_id(id);
  int endcap = (gem_id.region()>0 ? 1 : 2);
  int station = gem_id.station();
  CSCDetId csc_id(endcap, station, gem_id.ring(), gem_id.chamber(), 0);

  std::cout <<"extrapolateHsfromGEMPad gemid "<< gem_id <<" cscid "<< csc_id << std::endl;
  const CSCChamber* cscChamber(getCSCGeometry()->chamber(csc_id));
  const CSCLayer* cscKeyLayer(cscChamber->layer(3));
  const CSCLayerGeometry* cscKeyLayerGeometry(cscKeyLayer->geometry());

  const GEMSuperChamber* gemSuperChamber(getGEMGeometry()->superChamber(id));
  const GEMChamber* gemChamber(gemSuperChamber->chamber(1));
  auto gemRoll(gemChamber->etaPartition(2));//any roll
  const int nGEMPads(gemRoll->npads());
  std::cout <<"total GEMPads in roll 2 "<< nGEMPads << std::endl;
  if (gempad > nGEMPads or gempad < 0) result = -1;

  const LocalPoint lpGEM(gemRoll->centreOfPad(gempad));
  const GlobalPoint gp(gemRoll->toGlobal(lpGEM));
  const LocalPoint lpCSC(cscKeyLayer->toLocal(gp));
  const float strip(cscKeyLayerGeometry->strip(lpCSC));
  // HS are wrapped-around
  result = (int) (strip - 0.25)/0.5;
  return result;
}


int 
GEMDigiMatcher::extrapolateHsfromGEMStrip(unsigned int id, int gemstrip) const
{
  int result = -1 ;
  
  GEMDetId gem_id(id);//chamberid
  int endcap = (gem_id.region()>0 ? 1 : 2);
  int station = gem_id.station();
  CSCDetId csc_id(endcap, station, gem_id.ring(), gem_id.chamber(), 0);
  std::cout <<"extrapolateHsfromGEMStrip gemid "<< gem_id <<" cscid "<< csc_id << std::endl;

  const CSCChamber* cscChamber(getCSCGeometry()->chamber(csc_id));
  const CSCLayer* cscKeyLayer(cscChamber->layer(3));
  const CSCLayerGeometry* cscKeyLayerGeometry(cscKeyLayer->geometry());

  const GEMSuperChamber* gemSuperChamber(getGEMGeometry()->superChamber(id));
  const GEMChamber* gemChamber(gemSuperChamber->chamber(1));
  auto gemRoll(gemChamber->etaPartition(2));//any roll
  const int nGEMStrips(gemRoll->nstrips());
  if (gemstrip > nGEMStrips or gemstrip < 0) result = -1;

  const LocalPoint lpGEM(gemRoll->centreOfStrip(gemstrip));
  const GlobalPoint gp(gemRoll->toGlobal(lpGEM));
  const LocalPoint lpCSC(cscKeyLayer->toLocal(gp));
  const float strip(cscKeyLayerGeometry->strip(lpCSC));
  // HS are wrapped-around
  result = (int) (strip - 0.25)/0.5;
  return result;
}

std::vector<GlobalPoint> 
GEMDigiMatcher::positionPad1InDetId(unsigned int id) const
{
  std::vector<GlobalPoint> result;
  bool verbose = false;
  GEMDetId gem_id(id);
  if (verbose) std::cout << "In function positionPad2InDetId gem_id " << gem_id << std::endl;
  for (auto p: gemDigisInDetId(id)){
    LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfStrip(p.strip());
    GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
    // check if GP is already in container
    if (std::find(result.begin(), result.end(), gem_gp) == result.end()) result.push_back(gem_gp);
  }
  return result;
}

std::vector<GlobalPoint> 
GEMDigiMatcher::positionPad2InDetId(unsigned int id) const
{
  bool verbose = false;
  std::vector<GlobalPoint> result;
  GEMDetId gem_id(id);
  if (verbose) std::cout << "In function positionPad2InDetId gem_id " << gem_id << std::endl;
  for (auto p: gemDigisInDetId(id)){
    if (verbose) std::cout << "Strip " << p.strip() << std::endl;
    float middleStripOfPad = 0;
    if (p.strip()%2==0){
      middleStripOfPad = p.strip() - 1.;
    }
    else{
      middleStripOfPad = p.strip() + 0.;
    }
    if (verbose) std::cout << "middle strip " << middleStripOfPad << std::endl;
    LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfStrip(middleStripOfPad);
    GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
    if (std::find(result.begin(), result.end(), gem_gp) == result.end()) result.push_back(gem_gp);
  }
  return result;
}

std::vector<GlobalPoint> 
GEMDigiMatcher::positionPad4InDetId(unsigned int id) const
{
  bool verbose = false;
  std::vector<GlobalPoint> result;
  GEMDetId gem_id(id);
  if (verbose) std::cout << "In function positionPad4InDetId gem_id " << gem_id << std::endl;
  for (auto p: gemDigisInDetId(id)){
    if (verbose) std::cout << "Strip " << p.strip() << std::endl;
    float middleStripOfPad = 0;
    if (p.strip()%4==0){
      middleStripOfPad = p.strip() - 2.;
    }
    else if (p.strip()%4==3){
      middleStripOfPad = p.strip() - 1.;
    }
    else if (p.strip()%4==2){
      middleStripOfPad = p.strip() + 0.;
    }
    else{
      middleStripOfPad = p.strip() + 1.;
    }
    if (verbose) std::cout << "middle strip " << middleStripOfPad << std::endl;
    LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfStrip(middleStripOfPad);
    GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
    if (std::find(result.begin(), result.end(), gem_gp) == result.end()) result.push_back(gem_gp);
  }
  return result;
}

std::vector<GlobalPoint> 
GEMDigiMatcher::positionPad8InDetId(unsigned int id) const
{
  bool verbose = false;
  std::vector<GlobalPoint> result;
  GEMDetId gem_id(id);
  if (verbose) std::cout << "In function positionPad8InDetId gem_id " << gem_id << std::endl;
  for (auto p: gemDigisInDetId(id)){
    if (verbose) std::cout << "Strip " << p.strip() << std::endl;
    float middleStripOfPad = 0;
    if (p.strip()%8==0){
      middleStripOfPad = p.strip() - 4.;
    }
    else if (p.strip()%8==7){
      middleStripOfPad = p.strip() - 3.;
    }
    else if (p.strip()%8==6){
      middleStripOfPad = p.strip() - 2.;
    }
    else if (p.strip()%8==5){
      middleStripOfPad = p.strip() - 1.;
    }
    else if (p.strip()%8==4){
      middleStripOfPad = p.strip() + 0.;
    }
    else if (p.strip()%8==3){
      middleStripOfPad = p.strip() + 1.;
    }
    else if (p.strip()%8==2){
      middleStripOfPad = p.strip() + 2.;
    }
    else{
      middleStripOfPad = p.strip() + 3.;
    }
    if (verbose) std::cout << "middle strip " << middleStripOfPad << std::endl;
    LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfStrip(middleStripOfPad);
    GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
    if (std::find(result.begin(), result.end(), gem_gp) == result.end()) result.push_back(gem_gp);
  }
  return result;
}

GlobalPoint 
GEMDigiMatcher::getGlobalPointDigi(unsigned int rawId, const GEMDigi& d) const
{
  GEMDetId gem_id(rawId);
  LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfStrip(d.strip());
  GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
  return gem_gp;
}

GlobalPoint 
GEMDigiMatcher::getGlobalPointPad(unsigned int rawId, const GEMPadDigi& tp) const
{
  GEMDetId gem_id(rawId);
  LocalPoint gem_lp = getGEMGeometry()->etaPartition(gem_id)->centreOfPad(tp.pad());
  GlobalPoint gem_gp = getGEMGeometry()->idToDet(gem_id)->surface().toGlobal(gem_lp);
  return gem_gp;
}
