#include "GEMCode/GEMValidation/interface/DTRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;


DTRecHitMatcher::DTRecHitMatcher(SimHitMatcher& sh)
  : BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
  , simhit_matcher_(&sh)
{
  auto dtRecHit1DPair = conf().getParameter<edm::ParameterSet>("dtRecHit");
  dtRecHit1DPairInput_ = dtRecHit1DPair.getParameter<vector<edm::InputTag>>("validInputTags");
  maxBXDTRecHit1DPair_ = dtRecHit1DPair.getParameter<int>("maxBX");
  minBXDTRecHit1DPair_ = dtRecHit1DPair.getParameter<int>("minBX");
  verboseDTRecHit1DPair_ = dtRecHit1DPair.getParameter<int>("verbose");
  runDTRecHit1DPair_ = dtRecHit1DPair.getParameter<bool>("run");

  auto dtSegment2D = conf().getParameter<edm::ParameterSet>("dtRecSegment2D");
  dtRecSegment2DInput_ = dtSegment2D.getParameter<vector<edm::InputTag>>("validInputTags");
  maxBXDTRecSegment2D_ = dtSegment2D.getParameter<int>("maxBX");
  minBXDTRecSegment2D_ = dtSegment2D.getParameter<int>("minBX");
  verboseDTRecSegment2D_ = dtSegment2D.getParameter<int>("verbose");
  runDTRecSegment2D_ = dtSegment2D.getParameter<bool>("run");

  auto dtSegment4D = conf().getParameter<edm::ParameterSet>("dtRecSegment4D");
  dtRecSegment4DInput_ = dtSegment4D.getParameter<vector<edm::InputTag>>("validInputTags");
  maxBXDTRecSegment4D_ = dtSegment4D.getParameter<int>("maxBX");
  minBXDTRecSegment4D_ = dtSegment4D.getParameter<int>("minBX");
  verboseDTRecSegment4D_ = dtSegment4D.getParameter<int>("verbose");
  runDTRecSegment4D_ = dtSegment4D.getParameter<bool>("run");

  if (hasDTGeometry_) {
    edm::Handle<DTRecHitCollection> dt_rechits;
    if (gemvalidation::getByLabel(dtRecHit1DPairInput_, dt_rechits, event())) if (runDTRecHit1DPair_) matchDTRecHit1DPairsToSimTrack(*dt_rechits.product());

    edm::Handle<DTRecSegment2DCollection> dt_2DSegments;
    if (gemvalidation::getByLabel(dtRecSegment2DInput_, dt_2DSegments, event())) if (runDTRecSegment2D_) matchDTRecSegment2DsToSimTrack(*dt_2DSegments.product());

    edm::Handle<DTRecSegment4DCollection> dt_4DSegments;
    if (gemvalidation::getByLabel(dtRecSegment4DInput_, dt_4DSegments, event())) if (runDTRecSegment4D_) matchDTRecSegment4DsToSimTrack(*dt_4DSegments.product());
  }
}


void 
DTRecHitMatcher::matchDTRecHit1DPairsToSimTrack(const DTRecHitCollection& rechits)
{
  if (verboseDTRecHit1DPair_) cout << "Matching simtrack to DT rechits" << endl;
  // fetch all chamberIds with simhits
  auto layer_ids = simhit_matcher_->layerIdsDT();
  
  for (auto id: layer_ids) {
    DTLayerId p_id(id);
    
    // print all the wires in the DTChamber    
    auto hit_wires(simhit_matcher_->hitWiresInDTLayerId(id));
    if (verboseDTRecHit1DPair_) {
      cout<<"hit wires dt from simhit"<<endl;
      for (auto wire: hit_wires) cout << "\t"<<DTWireId(wire) << endl;
      cout<<endl;
    }
    
    // get the segments
    auto rechits_in_det = rechits.get(p_id);
    
    for (auto d = rechits_in_det.first; d != rechits_in_det.second; ++d) {
      if (verboseDTRecHit1DPair_) cout<<"rechit "<<p_id<<" "<<*d<<endl;
      
      unsigned int rightId(d->componentRecHit(DTEnums::DTCellSide::Right)->wireId());
      unsigned int leftId(d->componentRecHit(DTEnums::DTCellSide::Left)->wireId());
      set<unsigned int> rechitsIds;
      rechitsIds.insert(rightId);
      rechitsIds.insert(leftId);

      int wiresFound = 0;
      if (verboseDTRecHit1DPair_)cout<< rechitsIds.size() << " hit wires dt from rechit "<<endl;
      for (auto rh: rechitsIds) {
	if (verboseDTRecHit1DPair_)cout << "\t"<<rh << " " << DTWireId(rh) << endl;
	// is this "rechit" wire also a "simhit wire"?
	if (hit_wires.find(rh) != hit_wires.end()) ++wiresFound;

      }
      if (verboseDTRecHit1DPair_)cout << "Found " << wiresFound << " rechit wires out of " << hit_wires.size() << " simhit wires" << endl;
      if (wiresFound==0) continue;

      layer_to_dtRecHit1DPair_[p_id.rawId()].push_back(*d);
      superLayer_to_dtRecHit1DPair_[p_id.superlayerId().rawId()].push_back(*d);
      chamber_to_dtRecHit1DPair_[p_id.chamberId().rawId()].push_back(*d);
    }
  }
  if (verboseDTRecHit1DPair_) {
    cout << "Number of matching DTRecHit1DPairs: " << nDTRecHit1DPairs() << endl;    
    for (auto id: chamberIdsDTRecHit1DPair()) {
      for (auto rh: dtRecHit1DPairInChamber(id)) {
	cout << "\t" << id << " " << rh << endl;	
      }
    }
    cout<<endl;
  }
}


void
DTRecHitMatcher::matchDTRecSegment2DsToSimTrack(const DTRecSegment2DCollection& dtRecSegment2Ds)
{
}


void
DTRecHitMatcher::matchDTRecSegment4DsToSimTrack(const DTRecSegment4DCollection& dtRecSegment4Ds)
{
  if (verboseDTRecSegment4D_) cout << "Matching simtrack to segments" << endl;
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
      if (verboseDTRecSegment4D_) cout << "Found " << wiresFound << " rechit wires out of " << hit_wires.size() << " simhit wires" << endl;
      if (wiresFound==0) continue;
      if (verboseDTRecSegment4D_) cout << "Matching DTRecSegment4D " << *d << endl;

      chamber_to_dtRecSegment4D_[ p_id.rawId() ].push_back(*d);
    }
  }
  if (verboseDTRecSegment4D_) {
    cout << "Number of matching DTRecSegment4D: " << nDTRecSegment4Ds() << endl;
    for (auto id: chamberIdsDTRecSegment4D()) {
      for (auto rh: dtRecSegment4DInChamber(id)) {
	cout << "\t" << id << " " << rh << endl;	
      }
    }
    cout<<endl;
  }
}


set<unsigned int> 
DTRecHitMatcher::layerIdsDTRecHit1DPair() const 
{
  set<unsigned int> result;
  for (auto& p: layer_to_dtRecHit1DPair_) result.insert(p.first);
  return result;
}


set<unsigned int> 
DTRecHitMatcher::superLayerIdsDTRecHit1DPair() const
{
  set<unsigned int> result;
  for (auto& p: superLayer_to_dtRecHit1DPair_) result.insert(p.first);
  return result;
}


set<unsigned int>
DTRecHitMatcher::chamberIdsDTRecHit1DPair() const
{
  set<unsigned int> result;
  for (auto& p: chamber_to_dtRecHit1DPair_) result.insert(p.first);
  return result;
}


set<unsigned int>
DTRecHitMatcher::superLayerIdsDTRecSegment2D() const
{
  set<unsigned int> result;
  for (auto& p: superLayer_to_dtRecSegment2D_) result.insert(p.first);
  return result;
}


set<unsigned int> 
DTRecHitMatcher::chamberIdsDTRecSegment2D() const
{
  set<unsigned int> result;
  for (auto& p: chamber_to_dtRecSegment2D_) result.insert(p.first);
  return result;
}


set<unsigned int> 
DTRecHitMatcher::chamberIdsDTRecSegment4D() const
{
  set<unsigned int> result;
  for (auto& p: chamber_to_dtRecSegment4D_) result.insert(p.first);
  return result;
}

const DTRecHitMatcher::DTRecHit1DPairContainer& 
DTRecHitMatcher::dtRecHit1DPairInLayer(unsigned int detid) const
{
  if (layer_to_dtRecHit1DPair_.find(detid) == layer_to_dtRecHit1DPair_.end()) return no_dtRecHit1DPairs_;
  return layer_to_dtRecHit1DPair_.at(detid);
}


const DTRecHitMatcher::DTRecHit1DPairContainer& 
DTRecHitMatcher::dtRecHit1DPairInSuperLayer(unsigned int detid) const
{
  if (superLayer_to_dtRecHit1DPair_.find(detid) == superLayer_to_dtRecHit1DPair_.end()) return no_dtRecHit1DPairs_;
  return superLayer_to_dtRecHit1DPair_.at(detid);
}


const DTRecHitMatcher::DTRecHit1DPairContainer& 
DTRecHitMatcher::dtRecHit1DPairInChamber(unsigned int detid) const
{
  if (chamber_to_dtRecHit1DPair_.find(detid) == chamber_to_dtRecHit1DPair_.end()) return no_dtRecHit1DPairs_;
  return chamber_to_dtRecHit1DPair_.at(detid);
}


const DTRecHitMatcher::DTRecSegment2DContainer& 
DTRecHitMatcher::dtRecSegment2DInSuperLayer(unsigned int detid) const
{
  if (superLayer_to_dtRecSegment2D_.find(detid) == superLayer_to_dtRecSegment2D_.end()) return no_dtRecSegment2Ds_;
  return superLayer_to_dtRecSegment2D_.at(detid);
}


const DTRecHitMatcher::DTRecSegment2DContainer& 
DTRecHitMatcher::dtRecSegment2DInChamber(unsigned int detid) const
{
  if (chamber_to_dtRecSegment2D_.find(detid) == chamber_to_dtRecSegment2D_.end()) return no_dtRecSegment2Ds_;
  return chamber_to_dtRecSegment2D_.at(detid);
}


const DTRecHitMatcher::DTRecSegment4DContainer&
DTRecHitMatcher::dtRecSegment4DInChamber(unsigned int detid) const
{
  if (chamber_to_dtRecSegment4D_.find(detid) == chamber_to_dtRecSegment4D_.end()) return no_dtRecSegment4Ds_;
  return chamber_to_dtRecSegment4D_.at(detid);
}


const DTRecHitMatcher::DTRecHit1DPairContainer
DTRecHitMatcher::dtRecHit1DPairs() const
{
  DTRecHit1DPairContainer result;
  for (auto id: chamberIdsDTRecHit1DPair()){
    auto rechitsInChamber(dtRecHit1DPairInChamber(id));
    result.insert(result.end(), rechitsInChamber.begin(), rechitsInChamber.end());
  }
  return result;
}


const DTRecHitMatcher::DTRecSegment2DContainer
DTRecHitMatcher::dtRecSegment2Ds() const
{
  DTRecSegment2DContainer result;
  for (auto id: chamberIdsDTRecSegment2D()){
    auto segmentsInChamber(dtRecSegment2DInChamber(id));
    result.insert(result.end(), segmentsInChamber.begin(), segmentsInChamber.end());
  }
  return result;
}


const DTRecHitMatcher::DTRecSegment4DContainer
DTRecHitMatcher::dtRecSegment4Ds() const
{
  DTRecSegment4DContainer result;
  for (auto id: chamberIdsDTRecSegment4D()){
    auto segmentsInChamber(dtRecSegment4DInChamber(id));
    result.insert(result.end(), segmentsInChamber.begin(), segmentsInChamber.end());
  }
  return result;
}


int 
DTRecHitMatcher::nDTRecSegment4DInChamber(unsigned int detid) const
{
  return dtRecSegment4DInChamber(detid).size();
}


int 
DTRecHitMatcher::nDTRecHit1DPairs() const
{
  int n = 0;
  auto ids = chamberIdsDTRecHit1DPair();
  for (auto id: ids) n += dtRecHit1DPairInChamber(id).size();
  return n;  
}


int 
DTRecHitMatcher::nDTRecSegment2Ds() const
{
  int n = 0;
  auto ids = chamberIdsDTRecSegment2D();
  for (auto id: ids) n += dtRecSegment2DInChamber(id).size();
  return n;  
}


int
DTRecHitMatcher::nDTRecSegment4Ds() const
{
  int n = 0;
  auto ids = chamberIdsDTRecSegment4D();
  for (auto id: ids) n += dtRecSegment4DInChamber(id).size();
  return n;  
}

bool 
DTRecHitMatcher::dtRecSegment4DInContainer(const DTRecSegment4D& thisSegment, const DTRecSegment4DContainer& c) const
{
  bool isSame = false;
  for (auto& segment: c) if (areDTRecSegment4DSame(thisSegment,segment)) isSame = true;
  return isSame;
}


bool 
DTRecHitMatcher::isDTRecSegment4DMatched(const DTRecSegment4D& thisSegment) const
{
  return dtRecSegment4DInContainer(thisSegment, dtRecSegment4Ds());
}


bool 
DTRecHitMatcher::areDTRecSegment4DSame(const DTRecSegment4D& dtrh1,const DTRecSegment4D& dtrh2) const
{
  return dtrh1.localDirection()==dtrh2.localDirection() and dtrh1.localPosition()==dtrh2.localPosition();
}

 
