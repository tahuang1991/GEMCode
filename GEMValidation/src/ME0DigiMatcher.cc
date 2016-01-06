#include "GEMCode/GEMValidation/interface/ME0DigiMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;
using namespace matching;

ME0DigiMatcher::ME0DigiMatcher(SimHitMatcher& sh)
: DigiMatcher(sh)
{
  auto me0Digi_= conf().getParameter<edm::ParameterSet>("me0DigiPreReco");
  me0DigiInput_ = me0Digi_.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBXME0_ = me0Digi_.getParameter<int>("minBX");
  maxBXME0_ = me0Digi_.getParameter<int>("maxBX");
  matchDeltaStrip_ = me0Digi_.getParameter<int>("matchDeltaStrip");
  verboseDigi_ = me0Digi_.getParameter<int>("verbose");
  runME0Digi_ = me0Digi_.getParameter<bool>("run");

  if (hasME0Geometry_) {
    edm::Handle<ME0DigiPreRecoCollection> me0_digis;
    if (gemvalidation::getByLabel(me0DigiInput_, me0_digis, event())) if (runME0Digi_) matchPreRecoDigisToSimTrack(*me0_digis.product());    
  }
}

ME0DigiMatcher::~ME0DigiMatcher() 
{}


void
ME0DigiMatcher::matchPreRecoDigisToSimTrack(const ME0DigiPreRecoCollection& digis)
{
  auto det_ids = simhit_matcher_->detIdsME0();
  for (auto id: det_ids)
  {
    ME0DetId p_id(id);

    auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    if (verboseDigi_)
    {
      cout<<"hit_strips_fat ";
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto digis_in_det = digis.get(ME0DetId(id));

    for (auto d = digis_in_det.first; d != digis_in_det.second; ++d)
    {
      if (verboseDigi_) cout<<"gdigi "<<p_id<<" "<<*d<<endl;
      // check that the digi is within BX range
      //      if (d->bx() < minBXME0_ || d->bx() > maxBXME0_) continue;
      // check that it matches a strip that was hit by SimHits from our track
      //      if (hit_strips.find(d->strip()) == hit_strips.end()) continue;
      if (verboseDigi_) cout<<"oki"<<endl;

      detid_to_digis_[id].push_back(*d);
      chamber_to_digis_[ p_id.chamberId().rawId() ].push_back(*d);

      //int pad_num = 1 + static_cast<int>( roll->padOfStrip(d->strip()) ); // d->strip() is int
      //digi_map[ make_pair(pad_num, d->bx()) ].push_back( d->strip() );
    }
  }
}


std::set<unsigned int>
ME0DigiMatcher::detIds() const
{
  std::set<unsigned int> result;
  for (auto& p: detid_to_digis_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
ME0DigiMatcher::chamberIds() const
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_digis_) result.insert(p.first);
  return result;
}

std::set<unsigned int>
ME0DigiMatcher::superChamberIds() const
{
  std::set<unsigned int> result;
  for (auto& p: superchamber_to_digis_) result.insert(p.first);
  return result;
}


const ME0DigiPreRecoContainer&
ME0DigiMatcher::digisInDetId(unsigned int detid) const
{
  if (detid_to_digis_.find(detid) == detid_to_digis_.end()) return no_me0_digis_;
  return detid_to_digis_.at(detid);
}

const ME0DigiPreRecoContainer&
ME0DigiMatcher::digisInChamber(unsigned int detid) const
{
  if (chamber_to_digis_.find(detid) == chamber_to_digis_.end()) return no_me0_digis_;
  return chamber_to_digis_.at(detid);
}

const ME0DigiPreRecoContainer&
ME0DigiMatcher::digisInSuperChamber(unsigned int detid) const
{
  if (superchamber_to_digis_.find(detid) == superchamber_to_digis_.end()) return no_me0_digis_;
  return superchamber_to_digis_.at(detid);
}

int
ME0DigiMatcher::nLayersWithDigisInSuperChamber(unsigned int detid) const
{
  set<int> layers;
  return layers.size();
}


std::set<int>
ME0DigiMatcher::stripNumbersInDetId(unsigned int detid) const
{
  set<int> result;
  // auto digis = digisInDetId(detid);
  // for (auto& d: digis)
  // {
  //   result.insert( digi_channel(d) );
  // }
  return result;
}

std::set<int>
ME0DigiMatcher::partitionNumbers() const
{
  std::set<int> result;

  auto detids = detIds();
  for (auto id: detids)
  {
    ME0DetId idd(id);
    result.insert( idd.roll() );
  }
  return result;
}
