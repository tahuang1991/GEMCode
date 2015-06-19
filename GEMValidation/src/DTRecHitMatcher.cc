#include "GEMCode/GEMValidation/interface/DTRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;


DTRecHitMatcher::DTRecHitMatcher(SimHitMatcher& sh)
  : BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
  , simhit_matcher_(&sh)
{
  auto dtSegment2D = conf().getParameter<edm::ParameterSet>("dtRecSegment2D");
  dtRecSegment2DInput_ = dtSegment2D.getParameter<std::vector<edm::InputTag>>("validInputTags");
  maxBXDTRecSegment2D_ = dtSegment2D.getParameter<int>("maxBX");
  minBXDTRecSegment2D_ = dtSegment2D.getParameter<int>("minBX");
  verboseDTRecSegment2D_ = dtSegment2D.getParameter<int>("verbose");
  runDTRecSegment2D_ = dtSegment2D.getParameter<bool>("run");

  auto dtSegment4D = conf().getParameter<edm::ParameterSet>("dtRecSegment4D");
  dtRecSegment4DInput_ = dtSegment4D.getParameter<std::vector<edm::InputTag>>("validInputTags");
  maxBXDTRecSegment4D_ = dtSegment4D.getParameter<int>("maxBX");
  minBXDTRecSegment4D_ = dtSegment4D.getParameter<int>("minBX");
  verboseDTRecSegment4D_ = dtSegment4D.getParameter<int>("verbose");
  runDTRecSegment4D_ = dtSegment4D.getParameter<bool>("run");

  if (hasDTGeometry_) {
    edm::Handle<DTRecSegment2DCollection> dt_2DSegments;
    if (gemvalidation::getByLabel(dtRecSegment2DInput_, dt_2DSegments, event())) if (verboseDTRecSegment2D_) matchDTRecSegment2DsToSimTrack(*dt_2DSegments.product());

    edm::Handle<DTRecSegment4DCollection> dt_4DSegments;
    if (gemvalidation::getByLabel(dtRecSegment4DInput_, dt_4DSegments, event())) if (verboseDTRecSegment4D_) matchDTRecSegment4DsToSimTrack(*dt_4DSegments.product());
  }
}


void
DTRecHitMatcher::matchDTRecSegment2DsToSimTrack(const DTRecSegment2DCollection& dtRecSegment2Ds)
{
}


void
DTRecHitMatcher::matchDTRecSegment4DsToSimTrack(const DTRecSegment4DCollection& dtRecSegment4Ds)
{
  // fetch all chamberIds with simhits
  auto chamber_ids = simhit_matcher_->chamberIdsDT();
  
  for (auto id: chamber_ids) {
    DTChamberId p_id(id);
    
    //    std::set<int> hitWiresInDTLayerId(unsigned int, int margin_n_wires = 0) const;  // DT
    
    
    // auto hit_strips = simhit_matcher_->hitStripsInDetId(id, matchDeltaStrip_);
    // if (verboseDTRecSegment4D_) {
    //   cout<<"hit_strips_fat ";
    //   copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
    //   cout<<endl;
    // 	}
    
    // get the segments
    auto segments_in_det = dtRecSegment4Ds.get(p_id);
    
    for (auto d = segments_in_det.first; d != segments_in_det.second; ++d) {
      if (verboseDTRecSegment4D_) cout<<"segment "<<p_id<<" "<<*d<<endl;
      // check that the rechit is within BX range
      //	bunch crossing is calculated from the 2D segments
      //	if (d->BunchX() < minBXGEM_ || d->BunchX() > maxBXGEM_) continue;
      // check that it matches a wire that was hit by SimHits from our track
      
      // int firstStrip = d->firstClusterStrip();
      // int cls = d->clusterSize();
      // bool stripFound = false;
      
      // for(int i = firstStrip; i < (firstStrip + cls); i++){
      
      //   if (hit_strips.find(i) != hit_strips.end()) stripFound = true;
      //   //std::cout<<i<<" "<<firstStrip<<" "<<cls<<" "<<stripFound<<std::endl;
      
      // }
      
      // if (!stripFound) continue;
      // if (verboseGEMRecHit_) cout<<"oki"<<endl;
      
      // auto myrechit = make_digi(id, d->firstClusterStrip(), d->BunchX(), GEM_STRIP, d->clusterSize());
      // detid_to_recHits_[id].push_back(myrechit);
      // chamber_to_recHits_[ p_id.chamberId().rawId() ].push_back(myrechit);
      // superchamber_to_recHits_[ superch_id() ].push_back(myrechit);
    }
  }
}


// std::set<unsigned int>
// DTRecHitMatcher::superLayerIdsDTRecSegment2D() const
// {
//   std::set<unsigned int> result;
//   for (auto& p: chamber_to_dtRecSegment4D_) result.insert(p.first);
//   return result;
// }


const DTRecHitMatcher::DTRecSegment4DContainer&
DTRecHitMatcher::dtRecSegment4DInChamber(unsigned int detid) const
{
  if (chamber_to_dtRecSegment4D_.find(detid) == chamber_to_dtRecSegment4D_.end()) return no_dtRecSegment4Ds_;
  return chamber_to_dtRecSegment4D_.at(detid);
}


int 
DTRecHitMatcher::nDTRecSegment4DInChamber(unsigned int detid) const
{
  return dtRecSegment4DInChamber(detid).size();
}

