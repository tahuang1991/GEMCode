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
    if (gemvalidation::getByLabel(dtRecSegment2DInput_, dt_2DSegments, event())) if (runDTRecSegment2D_) matchDTRecSegment2DsToSimTrack(*dt_2DSegments.product());

    edm::Handle<DTRecSegment4DCollection> dt_4DSegments;
    if (gemvalidation::getByLabel(dtRecSegment4DInput_, dt_4DSegments, event())) if (runDTRecSegment4D_) matchDTRecSegment4DsToSimTrack(*dt_4DSegments.product());
  }
}


void
DTRecHitMatcher::matchDTRecSegment2DsToSimTrack(const DTRecSegment2DCollection& dtRecSegment2Ds)
{
}


void
DTRecHitMatcher::matchDTRecSegment4DsToSimTrack(const DTRecSegment4DCollection& dtRecSegment4Ds)
{
  cout << "Matching simtrack to segments" << endl;
  // fetch all chamberIds with simhits
  auto chamber_ids = simhit_matcher_->chamberIdsDT();
  
  for (auto id: chamber_ids) {
    DTChamberId p_id(id);
    
    // print all the wires in the DTChamber    
    auto hit_wires(simhit_matcher_->hitWiresInDTChamberId(id));
    if (verboseDTRecSegment4D_) {
      cout<<"hit wires dt from simhit"<<endl;
      for (auto wire: hit_wires) cout << "\t"<<DTWireId(wire) << endl;
      cout<<endl;
    }
    
    // get the segments
    auto segments_in_det = dtRecSegment4Ds.get(p_id);
    
    for (auto d = segments_in_det.first; d != segments_in_det.second; ++d) {
      if (verboseDTRecSegment4D_) cout<<"segment "<<p_id<<" "<<*d<<endl;
      
      const float time(d->hasPhi()? d->phiSegment()->t0() : d->zSegment()->t0());
      if (verboseDTRecSegment4D_) cout << "time " << time << endl;
      // check that the rechit is within BX range
      // bunch crossing is calculated from the 2D segments
      //      if (d->BunchX() < minBXGEM_ || d->BunchX() > maxBXGEM_) continue;
      // check that it matches a wire that was hit by SimHits from our track
      
      // access the rechits in the 4D segment
      vector<DTRecHit1D> recHits;
      if (d->hasPhi()) {
	vector<DTRecHit1D> phiHits = d->phiSegment()->specificRecHits();
	recHits.insert(recHits.end(), phiHits.begin(), phiHits.end());
      }
      if (d->hasZed()) {
	vector<DTRecHit1D> zedHits = d->zSegment()->specificRecHits();
	recHits.insert(recHits.end(), zedHits.begin(), zedHits.end());
      }

      int wiresFound = 0;
      if (verboseDTRecSegment4D_)cout<< recHits.size() << " hit wires dt from segment "<<endl;
      for (auto rh: recHits) {
	auto rhid(rh.wireId());
	if (verboseDTRecSegment4D_)cout << "\t"<<rh << " " << rhid << endl;
	// is this "rechit" wire also a "simhit wire"?
	if (hit_wires.find(rhid.rawId()) != hit_wires.end()) ++wiresFound;

      }
      if (verboseDTRecSegment4D_)cout << "Found " << wiresFound << " rechit wires out of " << hit_wires.size() << " simhit wires" << endl;
      if (wiresFound==0) continue;

      chamber_to_dtRecSegment4D_[ p_id.rawId() ].push_back(*d);
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

