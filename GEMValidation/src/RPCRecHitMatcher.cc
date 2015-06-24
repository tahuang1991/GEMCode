#include "GEMCode/GEMValidation/interface/RPCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;
using namespace matching;

RPCRecHitMatcher::RPCRecHitMatcher(SimHitMatcher& sh)
  : BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
  , simhit_matcher_(&sh)
{
  auto rpcRecHit_= conf().getParameter<edm::ParameterSet>("rpcRecHit");
  rpcRecHitInput_ = rpcRecHit_.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBXRPC_ = rpcRecHit_.getParameter<int>("minBX");
  maxBXRPC_ = rpcRecHit_.getParameter<int>("maxBX");
  matchDeltaStrip_ = rpcRecHit_.getParameter<int>("matchDeltaStrip");
  verboseRPCRecHit_ = rpcRecHit_.getParameter<int>("verbose");
  runRPCRecHit_ = rpcRecHit_.getParameter<bool>("run");

  if (hasRPCGeometry_) {
    edm::Handle<RPCRecHitCollection> rpc_rechits;
    if (gemvalidation::getByLabel(rpcRecHitInput_, rpc_rechits, event())) if (runRPCRecHit_) matchRecHitsToSimTrack(*rpc_rechits.product());
  }
}

RPCRecHitMatcher::~RPCRecHitMatcher() {}


void
RPCRecHitMatcher::matchRecHitsToSimTrack(const RPCRecHitCollection& rechits)
{
  
  auto det_ids = simhit_matcher_->detIdsRPC();
  for (auto id: det_ids)
  {
    RPCDetId p_id(id);

    auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    if (verboseRPCRecHit_)
    {
      cout<<"hit_strips_fat ";
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }

    auto rechits_in_det = rechits.get(RPCDetId(id));
    for (auto d = rechits_in_det.first; d != rechits_in_det.second; ++d)
    {
      if (verboseRPCRecHit_) cout<<"recHit "<<p_id<<" "<<*d<<endl;
      // check that the rechit is within BX range
      if (d->BunchX() < minBXRPC_ || d->BunchX() > maxBXRPC_) continue;
      // check that it matches a strip that was hit by SimHits from our track

      int firstStrip = d->firstClusterStrip();
      int cls = d->clusterSize();
      bool stripFound = false;

      for(int i = firstStrip; i < (firstStrip + cls); i++){

	if (hit_strips.find(i) != hit_strips.end()) stripFound = true;
        //std::cout<<i<<" "<<firstStrip<<" "<<cls<<" "<<stripFound<<std::endl;
	
      }

      if (!stripFound) continue;
      if (verboseRPCRecHit_) cout<<"oki"<<endl;

      auto myrechit = make_digi(id, d->firstClusterStrip(), d->BunchX(), RPC_STRIP, d->clusterSize());
      detid_to_recHits_[id].push_back(myrechit);
      chamber_to_recHits_[ p_id.chamberId().rawId() ].push_back(myrechit);

    }
  }
  
}


std::set<unsigned int>
RPCRecHitMatcher::detIds() const
{
  std::set<unsigned int> result;
  for (auto& p: detid_to_recHits_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
RPCRecHitMatcher::chamberIds() const
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_recHits_) result.insert(p.first);
  return result;
}

const RPCRecHitMatcher::RecHitContainer&
RPCRecHitMatcher::recHitsInDetId(unsigned int detid) const
{
  if (detid_to_recHits_.find(detid) == detid_to_recHits_.end()) return no_recHits_;
  return detid_to_recHits_.at(detid);
}

const RPCRecHitMatcher::RecHitContainer&
RPCRecHitMatcher::recHitsInChamber(unsigned int detid) const
{
  if (chamber_to_recHits_.find(detid) == chamber_to_recHits_.end()) return no_recHits_;
  return chamber_to_recHits_.at(detid);
}

std::set<int>
RPCRecHitMatcher::stripNumbersInDetId(unsigned int detid) const
{
  set<int> result;
  /*
  auto recHits = recHitsInDetId(detid);
  for (auto& d: recHits)
  {
    result.insert( digi_channel(d) );
  }
  */
  return result;
}

std::set<int>
RPCRecHitMatcher::partitionNumbers() const
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

GlobalPoint
RPCRecHitMatcher::recHitPosition(const RecHit& rechit) const
{
  unsigned int id = digi_id(rechit);
  int strip = digi_channel(rechit);
  DigiType t = digi_type(rechit);

  GlobalPoint gp;
  if ( t == RPC_STRIP )
  {
    RPCDetId idd(id);
    LocalPoint lp = getRPCGeometry()->roll(idd)->centreOfStrip(strip);
    gp = getRPCGeometry()->idToDet(id)->surface().toGlobal(lp);
  }

  return gp;
}


GlobalPoint
RPCRecHitMatcher::recHitMeanPosition(const RecHitContainer& rechit) const
{
  GlobalPoint point_zero;
  if (rechit.empty()) return point_zero; // point "zero"

  float sumx, sumy, sumz;
  sumx = sumy = sumz = 0.f;
  size_t n = 0;
  for (auto& d: rechit)
  {
    GlobalPoint gp = recHitPosition(d);
    if (gp == point_zero) continue;

    sumx += gp.x();
    sumy += gp.y();
    sumz += gp.z();
    ++n;
  }
  if (n == 0) return GlobalPoint();
  return GlobalPoint(sumx/n, sumy/n, sumz/n);
}

