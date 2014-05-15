#include "RPCDigiMatcher.h"
#include "SimHitMatcher.h"

using namespace std;
using namespace matching;


RPCDigiMatcher::RPCDigiMatcher(SimHitMatcher& sh)
: DigiMatcher(sh)
{
  auto rpcDigi_= conf().getParameter<edm::ParameterSet>("rpcStripDigi");
  rpcDigiInput_ = rpcDigi_.getParameter<edm::InputTag>("input");
  minBXRPC_ = rpcDigi_.getParameter<int>("minBX");
  maxBXRPC_ = rpcDigi_.getParameter<int>("maxBX");
  matchDeltaStrip_ = rpcDigi_.getParameter<int>("matchDeltaStrip");
  verboseDigi_ = rpcDigi_.getParameter<int>("verbose");

  matchDeltaStrip_ = conf().getUntrackedParameter<int>("matchDeltaStripRPC", 1);

  setVerbose(conf().getUntrackedParameter<int>("verboseRPCDigi", 0));

  if (!rpcDigiInput_.label().empty())
  {
    init();
  }
}

RPCDigiMatcher::~RPCDigiMatcher() {}


void
RPCDigiMatcher::init()
{
  edm::Handle<RPCDigiCollection> rpc_digis;
  event().getByLabel(rpcDigiInput_, rpc_digis);
  matchDigisToSimTrack(*rpc_digis.product());
}


void
RPCDigiMatcher::matchDigisToSimTrack(const RPCDigiCollection& digis)
{
  auto det_ids = simhit_matcher_->detIdsRPC();
  for (auto id: det_ids)
  {
    RPCDetId p_id(id);

    auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    if (verboseDigi_)
    {
      cout<<"hit_strips_fat ";
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto digis_in_det = digis.get(RPCDetId(id));

    for (auto d = digis_in_det.first; d != digis_in_det.second; ++d)
    {
      if (verboseDigi_) cout<<"gdigi "<<p_id<<" "<<*d<<endl;
      // check that the digi is within BX range
      if (d->bx() < minBXRPC_ || d->bx() > maxBXRPC_) continue;
      // check that it matches a strip that was hit by SimHits from our track
      if (hit_strips.find(d->strip()) == hit_strips.end()) continue;
      if (verboseDigi_) cout<<"oki"<<endl;

      auto mydigi = make_digi(id, d->strip(), d->bx(), RPC_STRIP);
      detid_to_digis_[id].push_back(mydigi);
      chamber_to_digis_[ p_id.chamberId().rawId() ].push_back(mydigi);
    }
  }
}


std::set<unsigned int>
RPCDigiMatcher::detIds() const
{
  std::set<unsigned int> result;
  for (auto& p: detid_to_digis_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
RPCDigiMatcher::chamberIds() const
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_digis_) result.insert(p.first);
  return result;
}

const matching::DigiContainer&
RPCDigiMatcher::digisInDetId(unsigned int detid) const
{
  if (detid_to_digis_.find(detid) == detid_to_digis_.end()) return no_digis_;
  return detid_to_digis_.at(detid);
}

const matching::DigiContainer&
RPCDigiMatcher::digisInChamber(unsigned int detid) const
{
  if (chamber_to_digis_.find(detid) == chamber_to_digis_.end()) return no_digis_;
  return chamber_to_digis_.at(detid);
}


std::set<int>
RPCDigiMatcher::stripNumbersInDetId(unsigned int detid) const
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
RPCDigiMatcher::partitionNumbers() const
{
  std::set<int> result;

  auto detids = detIds();
  for (auto id: detids)
  {
    RPCDetId idd(id);
    result.insert( idd.roll() );
  }
  return result;
}

