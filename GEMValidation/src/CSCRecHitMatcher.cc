#include "GEMCode/GEMValidation/interface/CSCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;


CSCRecHitMatcher::CSCRecHitMatcher(SimHitMatcher& sh)
  : BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
  , simhit_matcher_(&sh)
{
  auto cscRecHit2D = conf().getParameter<edm::ParameterSet>("cscRecHit");
  cscRecHit2DInput_ = cscRecHit2D.getParameter<std::vector<edm::InputTag>>("validInputTags");
  maxBXCSCRecHit2D_ = cscRecHit2D.getParameter<int>("maxBX");
  minBXCSCRecHit2D_ = cscRecHit2D.getParameter<int>("minBX");
  verboseCSCRecHit2D_ = cscRecHit2D.getParameter<int>("verbose");
  runCSCRecHit2D_ = cscRecHit2D.getParameter<bool>("run");

  auto cscSegment2D = conf().getParameter<edm::ParameterSet>("cscSegment");
  cscSegmentInput_ = cscSegment2D.getParameter<std::vector<edm::InputTag>>("validInputTags");
  maxBXCSCSegment_ = cscSegment2D.getParameter<int>("maxBX");
  minBXCSCSegment_ = cscSegment2D.getParameter<int>("minBX");
  verboseCSCSegment_ = cscSegment2D.getParameter<int>("verbose");
  runCSCSegment_ = cscSegment2D.getParameter<bool>("run");

  if (hasCSCGeometry_) {
    edm::Handle<CSCRecHit2DCollection> csc_rechits;
    if (gemvalidation::getByLabel(cscRecHit2DInput_, csc_rechits, event())) if (runCSCRecHit2D_) matchCSCRecHit2DsToSimTrack(*csc_rechits.product());

    edm::Handle<CSCSegmentCollection> csc_2DSegments;
    if (gemvalidation::getByLabel(cscSegmentInput_, csc_2DSegments, event())) if (runCSCSegment_) matchCSCSegmentsToSimTrack(*csc_2DSegments.product());

  }
}


void 
CSCRecHitMatcher::matchCSCRecHit2DsToSimTrack(const CSCRecHit2DCollection& rechits)
{
  cout << "Matching simtrack to CSC rechits" << endl;
  // fetch all chamberIds with simhits
  auto layer_ids = simhit_matcher_->detIdsCSC();
  
  for (auto id: layer_ids) {
    CSCDetId p_id(id);
    
    // print all the wires in the CSCChamber    
    auto hit_wg(simhit_matcher_->hitWiregroupsInDetId(id));
    if (verboseCSCRecHit2D_) {
      cout<<"hit wg csc from simhit"<<endl;
      copy(hit_wg.begin(), hit_wg.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }
    
    // // get the segments
    // auto rechits_in_det = rechits.get(p_id);    
    // for (auto d = rechits_in_det.first; d != rechits_in_det.second; ++d) {
    //   if (verboseCSCRecHit2D_) cout<<"rechit "<<p_id<<" "<<*d<<endl;
      
    //   int wiresFound = 0;
    //   if (verboseCSCRecHit2D_) { 
    // 	cout<< rechitsIds.size() << " hit wires csc from rechit "<<endl;
    // 	cout << "\t"<<rh << " " << CSCWireId(rh) << endl;
    //   }
    //   // is this "rechit" wire also a "simhit wire"?
    //   if (hit_wires.find(rh) != hit_wires.end()) ++wiresFound;
    //   if (verboseCSCRecHit2D_)cout << "Found " << wiresFound << " rechit wires out of " << hit_wires.size() << " simhit wires" << endl;
    //   if (wiresFound==0) continue;
      
    //   layer_to_cscRecHit2D_[id].push_back(*d);
    //   chamber_to_cscRecHit2D_[p_id.chamberId().rawId()].push_back(*d);
    // }
  }
}


void
CSCRecHitMatcher::matchCSCSegmentsToSimTrack(const CSCSegmentCollection& cscRecSegments)
{
  /*
  cout << "Matching simtrack to segments" << endl;
  // fetch all chamberIds with simhits
  auto chamber_ids = simhit_matcher_->chamberIdsCSC();
  
  for (auto id: chamber_ids) {
    CSCChamberId p_id(id);
    
    // print all the wires in the CSCChamber    
    auto hit_wires(simhit_matcher_->hitWiresInCSCChamberId(id));
    if (verboseCSCRecSegment4D_) {
      cout<<"hit wires csc from simhit"<<endl;
      for (auto wire: hit_wires) cout << "\t"<<CSCWireId(wire) << endl;
      cout<<endl;
    }
    
    // get the segments
    auto segments_in_det = cscRecSegment4Ds.get(p_id);
    
    for (auto d = segments_in_det.first; d != segments_in_det.second; ++d) {
      if (verboseCSCRecSegment4D_) cout<<"segment "<<p_id<<" "<<*d<<endl;
      
      const float time(d->hasPhi()? d->phiSegment()->t0() : d->zSegment()->t0());
      if (verboseCSCRecSegment4D_) cout << "time " << time << endl;
      // check that the rechit is within BX range
      // bunch crossing is calculated from the 2D segments
      //      if (d->BunchX() < minBXGEM_ || d->BunchX() > maxBXGEM_) continue;
      // check that it matches a wire that was hit by SimHits from our track
      
      // access the rechits in the 4D segment
      vector<CSCRecHit1D> recHits;
      if (d->hasPhi()) {
	vector<CSCRecHit1D> phiHits = d->phiSegment()->specificRecHits();
	recHits.insert(recHits.end(), phiHits.begin(), phiHits.end());
      }
      if (d->hasZed()) {
	vector<CSCRecHit1D> zedHits = d->zSegment()->specificRecHits();
	recHits.insert(recHits.end(), zedHits.begin(), zedHits.end());
      }

      int wiresFound = 0;
      if (verboseCSCRecSegment4D_)cout<< recHits.size() << " hit wires csc from segment "<<endl;
      for (auto rh: recHits) {
	auto rhid(rh.wireId());
	if (verboseCSCRecSegment4D_)cout << "\t"<<rh << " " << rhid << endl;
	// is this "rechit" wire also a "simhit wire"?
	if (hit_wires.find(rhid.rawId()) != hit_wires.end()) ++wiresFound;

      }
      if (verboseCSCRecSegment4D_)cout << "Found " << wiresFound << " rechit wires out of " << hit_wires.size() << " simhit wires" << endl;
      if (wiresFound==0) continue;

      chamber_to_cscRecSegment4D_[ p_id.rawId() ].push_back(*d);
    }
  }
  */
}


std::set<unsigned int> 
CSCRecHitMatcher::layerIdsCSCRecHit2D() const 
{
  std::set<unsigned int> result;
  for (auto& p: layer_to_cscRecHit2D_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
CSCRecHitMatcher::chamberIdsCSCRecHit2D() const 
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_cscRecHit2D_) result.insert(p.first);
  return result;
}


std::set<unsigned int> 
CSCRecHitMatcher::chamberIdsCSCSegment() const
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_cscSegment_) result.insert(p.first);
  return result;
}


const CSCRecHitMatcher::CSCRecHit2DContainer& 
CSCRecHitMatcher::cscRecHit2DsInLayer(unsigned int detid) const
{
  if (layer_to_cscRecHit2D_.find(detid) == layer_to_cscRecHit2D_.end()) return no_cscRecHit2Ds_;
  return layer_to_cscRecHit2D_.at(detid);
}

 
const CSCRecHitMatcher::CSCRecHit2DContainer& 
CSCRecHitMatcher::cscRecHit2DsInChamber(unsigned int detid) const
{
  if (chamber_to_cscRecHit2D_.find(detid) == chamber_to_cscRecHit2D_.end()) return no_cscRecHit2Ds_;
  return chamber_to_cscRecHit2D_.at(detid);
}


const CSCRecHitMatcher::CSCSegmentContainer& 
CSCRecHitMatcher::cscSegmentsInChamber(unsigned int detid) const
{
  if (chamber_to_cscSegment_.find(detid) == chamber_to_cscSegment_.end()) return no_cscSegments_;
  return chamber_to_cscSegment_.at(detid);
}


int 
CSCRecHitMatcher::nCSCRecHit2DsInLayer(unsigned int detid) const
{
  return cscRecHit2DsInLayer(detid).size();
}


int 
CSCRecHitMatcher::nCSCRecHit2DsInChamber(unsigned int detid) const
{
  return cscRecHit2DsInChamber(detid).size();
}


int 
CSCRecHitMatcher::nCSCSegmentsInChamber(unsigned int detid) const
{
  return cscSegmentsInChamber(detid).size();
}
