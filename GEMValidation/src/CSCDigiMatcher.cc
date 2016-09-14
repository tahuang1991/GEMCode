#include "GEMCode/GEMValidation/interface/CSCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;
using namespace matching;


CSCDigiMatcher::CSCDigiMatcher(SimHitMatcher& sh, edm::ConsumesCollector && iC)
: DigiMatcher(sh)
{
  auto cscWireDigi_ = conf().getParameter<edm::ParameterSet>("cscWireDigi");
  cscWireDigiInput_ = iC.consumes<CSCWireDigiCollection>(cscWireDigi_.getParameter<edm::InputTag>("validInputTags"));
  verboseWG_ = cscWireDigi_.getParameter<int>("verbose");
  minBXCSCWire_ = cscWireDigi_.getParameter<int>("minBX");
  maxBXCSCWire_ = cscWireDigi_.getParameter<int>("maxBX");
  matchDeltaWG_ = cscWireDigi_.getParameter<int>("matchDeltaWG");
  runWG_ = cscWireDigi_.getParameter<bool>("run");

  auto cscComparatorDigi_ = conf().getParameter<edm::ParameterSet>("cscStripDigi");
  cscComparatorDigiInput_ = iC.consumes<CSCComparatorDigiCollection>(cscComparatorDigi_.getParameter<edm::InputTag>("validInputTags"));
  verboseStrip_ = cscComparatorDigi_.getParameter<int>("verbose");
  minBXCSCComp_ = cscComparatorDigi_.getParameter<int>("minBX");
  maxBXCSCComp_ = cscComparatorDigi_.getParameter<int>("maxBX");
  matchDeltaStrip_ = cscComparatorDigi_.getParameter<int>("matchDeltaStrip");
  runStrip_ = cscComparatorDigi_.getParameter<bool>("run");

  if (hasCSCGeometry_) {
    edm::Handle<CSCComparatorDigiCollection> comp_digis;
    if(gemvalidation::getByToken(cscComparatorDigiInput_, comp_digis, event())) if (runWG_) matchStripsToSimTrack(*comp_digis.product());
    
    edm::Handle<CSCWireDigiCollection> wire_digis;
    if (gemvalidation::getByToken(cscWireDigiInput_, wire_digis, event())) if (runStrip_) matchWiresToSimTrack(*wire_digis.product());
  }
}


CSCDigiMatcher::~CSCDigiMatcher() {}


void
CSCDigiMatcher::matchStripsToSimTrack(const CSCComparatorDigiCollection& comparators)
{
  auto det_ids = simhit_matcher_->detIdsCSC(0);
  for (auto id: det_ids)
  {
    CSCDetId layer_id(id);

    auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    if (verboseStrip_)
    {
      cout<<"hit_strips_fat ";
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto comp_digis_in_det = comparators.get(layer_id);
    for (auto c = comp_digis_in_det.first; c != comp_digis_in_det.second; ++c)
    {
      if (verboseStrip_) cout<<"sdigi "<<layer_id<<" "<<*c<<endl;

      // check that the first BX for this digi wasn't too early or too late
      if (c->getTimeBin() < minBXCSCComp_ || c->getTimeBin() > maxBXCSCComp_) continue;

      int strip = c->getStrip(); // strips are counted from 1
      // check that it matches a strip that was hit by SimHits from our track
      if (hit_strips.find(strip) == hit_strips.end()) continue;
      if (verboseStrip_) cout<<"oki"<<endl;

      // get half-strip, counting from 1
      int half_strip = 2*strip - 1 + c->getComparator();

      auto mydigi = make_digi(id, half_strip, c->getTimeBin(), CSC_STRIP);
      detid_to_halfstrips_[id].push_back(mydigi);
      chamber_to_halfstrips_[ layer_id.chamberId().rawId() ].push_back(mydigi);

      detid_to_cschalfstrips_[id].push_back(*c);
      chamber_to_cschalfstrips_[ layer_id.chamberId().rawId() ].push_back(*c);
    }
  }
}


void
CSCDigiMatcher::matchWiresToSimTrack(const CSCWireDigiCollection& wires)
{
  auto det_ids = simhit_matcher_->detIdsCSC(0);
  for (auto id: det_ids)
  {
    CSCDetId layer_id(id);

    auto hit_wires = simhit_matcher_->hitWiregroupsInDetId(id, matchDeltaWG_);
    if (verboseWG_)
    {
      cout<<"hit_wires ";
      copy(hit_wires.begin(), hit_wires.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto wire_digis_in_det = wires.get(layer_id);
    for (auto w = wire_digis_in_det.first; w != wire_digis_in_det.second; ++w)
    {
      // check that the first BX for this digi wasn't too early or too late
      if (w->getTimeBin() < minBXCSCWire_ || w->getTimeBin() > maxBXCSCWire_) continue;

      int wg = w->getWireGroup(); // wiregroups are counted from 1
      // check that it matches a strip that was hit by SimHits from our track
      if (hit_wires.find(wg) == hit_wires.end()) continue;

      auto mydigi = make_digi(id, wg, w->getTimeBin(), CSC_WIRE);
      detid_to_wires_[id].push_back(mydigi);
      chamber_to_wires_[ layer_id.chamberId().rawId() ].push_back(mydigi);

      detid_to_cscwires_[id].push_back(*w);
      chamber_to_cscwires_[ layer_id.chamberId().rawId() ].push_back(*w);
    }
  }
}


std::set<unsigned int>
CSCDigiMatcher::selectDetIds(const Id2DigiContainer &digis, int csc_type) const
{
  std::set<unsigned int> result;
  for (auto& p: digis)
  {
    auto id = p.first;
    if (csc_type > 0)
    {
      CSCDetId detId(id);
      if (gemvalidation::toCSCType(detId.station(), detId.ring()) != csc_type) continue;
    }
    result.insert(p.first);
  }
  return result;
}


std::set<unsigned int>
CSCDigiMatcher::detIdsStrip(int csc_type) const
{
  return selectDetIds(detid_to_halfstrips_, csc_type);
}


std::set<unsigned int>
CSCDigiMatcher::detIdsWire(int csc_type) const
{
  return selectDetIds(detid_to_wires_, csc_type);
}


std::set<unsigned int>
CSCDigiMatcher::chamberIdsStrip(int csc_type) const
{
  return selectDetIds(chamber_to_halfstrips_, csc_type);
}


std::set<unsigned int>
CSCDigiMatcher::chamberIdsWire(int csc_type) const
{
  return selectDetIds(chamber_to_wires_, csc_type);
}


const matching::DigiContainer&
CSCDigiMatcher::stripDigisInDetId(unsigned int detid) const
{
  if (detid_to_halfstrips_.find(detid) == detid_to_halfstrips_.end()) return no_digis_;
  return detid_to_halfstrips_.at(detid);
}


const matching::DigiContainer&
CSCDigiMatcher::stripDigisInChamber(unsigned int detid) const
{
  if (chamber_to_halfstrips_.find(detid) == chamber_to_halfstrips_.end()) return no_digis_;
  return chamber_to_halfstrips_.at(detid);
}


const matching::DigiContainer&
CSCDigiMatcher::wireDigisInDetId(unsigned int detid) const
{
  if (detid_to_wires_.find(detid) == detid_to_wires_.end()) return no_digis_;
  return detid_to_wires_.at(detid);
}


const matching::DigiContainer&
CSCDigiMatcher::wireDigisInChamber(unsigned int detid) const
{
  if (chamber_to_wires_.find(detid) == chamber_to_wires_.end()) return no_digis_;
  return chamber_to_wires_.at(detid);
}


const CSCStripDigiContainer&
CSCDigiMatcher::cscStripDigisInDetId(unsigned int detid) const
{
  if (detid_to_cschalfstrips_.find(detid) == detid_to_cschalfstrips_.end()) return no_csc_strips_;
  return detid_to_cschalfstrips_.at(detid);
}


const CSCStripDigiContainer&
CSCDigiMatcher::cscStripDigisInChamber(unsigned int detid) const
{
  if (chamber_to_cschalfstrips_.find(detid) == chamber_to_cschalfstrips_.end()) return no_csc_strips_;
  return chamber_to_cschalfstrips_.at(detid);
}


const CSCWireDigiContainer&
CSCDigiMatcher::cscWireDigisInDetId(unsigned int detid) const
{
  if (detid_to_cscwires_.find(detid) == detid_to_cscwires_.end()) return no_csc_wires_;
  return detid_to_cscwires_.at(detid);
}


const CSCWireDigiContainer&
CSCDigiMatcher::cscWireDigisInChamber(unsigned int detid) const
{
  if (chamber_to_cscwires_.find(detid) == chamber_to_cscwires_.end()) return no_csc_wires_;
  return chamber_to_cscwires_.at(detid);
}


int
CSCDigiMatcher::nLayersWithStripInChamber(unsigned int detid) const
{
  set<int> layers_with_hits;
  auto digis = stripDigisInChamber(detid);
  for (auto& d: digis)
  {
    CSCDetId idd(digi_id(d));
    layers_with_hits.insert(idd.layer());
  }
  return layers_with_hits.size();
}


int
CSCDigiMatcher::nLayersWithWireInChamber(unsigned int detid) const
{
  set<int> layers_with_hits;
  auto digis = wireDigisInChamber(detid);
  for (auto& d: digis)
  {
    CSCDetId idd(digi_id(d));
    layers_with_hits.insert(idd.layer());
  }
  return layers_with_hits.size();
}


int
CSCDigiMatcher::nCoincidenceStripChambers(int min_n_layers) const
{
  int result = 0;
  auto chamber_ids = chamberIdsStrip();
  for (auto id: chamber_ids)
  {
    if (nLayersWithStripInChamber(id) >= min_n_layers) result += 1;
  }
  return result;
}


int
CSCDigiMatcher::nCoincidenceWireChambers(int min_n_layers) const
{
  int result = 0;
  auto chamber_ids = chamberIdsWire();
  for (auto id: chamber_ids)
  {
    if (nLayersWithWireInChamber(id) >= min_n_layers) result += 1;
  }
  return result;
}


std::set<int>
CSCDigiMatcher::stripsInDetId(unsigned int detid) const
{
  set<int> result;
  auto digis = stripDigisInDetId(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  return result;
}


std::set<int>
CSCDigiMatcher::wiregroupsInDetId(unsigned int detid) const
{
  set<int> result;
  auto digis = wireDigisInDetId(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  return result;
}


std::set<int>
CSCDigiMatcher::stripsInChamber(unsigned int detid, int max_gap_to_fill) const
{
  set<int> result;
  auto digis = stripDigisInChamber(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  if (max_gap_to_fill > 0)
  {
    int prev = -111;
    for (auto s: result)
    {
      //cout<<"gap "<<s<<" - "<<prev<<" = "<<s - prev<<"  added 0";
      if (s - prev > 1 && s - prev - 1 <= max_gap_to_fill)
      {
        //int sz = result.size();
        for (int fill_s = prev+1; fill_s < s; ++fill_s) result.insert(fill_s);
        //cout<<result.size() - sz;
      }
      //cout<<" elems"<<endl;
      prev = s;
    }
  }

  return result;
}

std::set<int>
CSCDigiMatcher::wiregroupsInChamber(unsigned int detid, int max_gap_to_fill) const
{
  set<int> result;
  auto digis = wireDigisInChamber(detid);
  for (auto& d: digis)
  {
    result.insert( digi_channel(d) );
  }
  if (max_gap_to_fill > 0)
  {
    int prev = -111;
    for (auto w: result)
    {
      if (w - prev > 1 && w - prev - 1 <= max_gap_to_fill)
      {
        for (int fill_w = prev+1; fill_w < w; ++fill_w) result.insert(fill_w);
      }
      prev = w;
    }
  }
  return result;
}
