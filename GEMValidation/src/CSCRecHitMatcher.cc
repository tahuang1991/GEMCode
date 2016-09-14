#include "GEMCode/GEMValidation/interface/CSCRecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;


CSCRecHitMatcher::CSCRecHitMatcher(SimHitMatcher& sh, edm::ConsumesCollector && iC)
  : BaseMatcher(sh.trk(), sh.vtx(), sh.conf(), sh.event(), sh.eventSetup())
  , simhit_matcher_(&sh)
{
  auto cscRecHit2D = conf().getParameter<edm::ParameterSet>("cscRecHit");
  cscRecHit2DInput_ = iC.consumes<CSCRecHit2DCollection>(cscRecHit2D.getParameter<edm::InputTag>("validInputTags"));
  maxBXCSCRecHit2D_ = cscRecHit2D.getParameter<int>("maxBX");
  minBXCSCRecHit2D_ = cscRecHit2D.getParameter<int>("minBX");
  verboseCSCRecHit2D_ = cscRecHit2D.getParameter<int>("verbose");
  runCSCRecHit2D_ = cscRecHit2D.getParameter<bool>("run");

  auto cscSegment2D = conf().getParameter<edm::ParameterSet>("cscSegment");
  cscSegmentInput_ = iC.consumes<CSCSegmentCollection>(cscSegment2D.getParameter<edm::InputTag>("validInputTags"));
  maxBXCSCSegment_ = cscSegment2D.getParameter<int>("maxBX");
  minBXCSCSegment_ = cscSegment2D.getParameter<int>("minBX");
  verboseCSCSegment_ = cscSegment2D.getParameter<int>("verbose");
  runCSCSegment_ = cscSegment2D.getParameter<bool>("run");

  if (hasCSCGeometry_) {
    edm::Handle<CSCRecHit2DCollection> csc_rechits;
    if (gemvalidation::getByToken(cscRecHit2DInput_, csc_rechits, event())) if (runCSCRecHit2D_) matchCSCRecHit2DsToSimTrack(*csc_rechits.product());

    edm::Handle<CSCSegmentCollection> csc_2DSegments;
    if (gemvalidation::getByToken(cscSegmentInput_, csc_2DSegments, event())) if (runCSCSegment_) matchCSCSegmentsToSimTrack(*csc_2DSegments.product());
  }
}


void 
CSCRecHitMatcher::matchCSCRecHit2DsToSimTrack(const CSCRecHit2DCollection& rechits)
{
  if (verboseCSCRecHit2D_) cout << "Matching simtrack to CSC rechits" << endl;
  // fetch all chamberIds with simhits
  auto layer_ids = simhit_matcher_->detIdsCSC(0);
  if (verboseCSCRecHit2D_) cout << "Number of matched csc layer_ids " << layer_ids.size() << endl;
  
  for (auto id: layer_ids) {
    CSCDetId p_id(id);
    
    // print all the wires in the CSCChamber    
    auto hit_wg(simhit_matcher_->hitWiregroupsInDetId(id));
    if (verboseCSCRecHit2D_) {
      cout<<"hit wg csc from simhit"<<endl;
      copy(hit_wg.begin(), hit_wg.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }
    
    // print all the strips in the CSCChamber    
    auto hit_strips(simhit_matcher_->hitStripsInDetId(id));
    if (verboseCSCRecHit2D_) {
      cout<<"hit strip csc from simhit"<<endl;
      copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }
    
    // get the rechits
    auto rechits_in_det = rechits.get(p_id);    
    for (auto d = rechits_in_det.first; d != rechits_in_det.second; ++d) {
      if (verboseCSCRecHit2D_) cout<<"rechit "<<p_id<<" "<<*d;
      
      const bool wireMatch(std::find(hit_wg.begin(), hit_wg.end(),d->hitWire())!=hit_wg.end());
      bool stripMatch(false);
      for(size_t iS =0;iS< d->nStrips();++iS) {
	if (std::find(hit_strips.begin(), hit_strips.end(),d->hitWire())!=hit_wg.end()) stripMatch = true;
      }
      // this rechit was matched to a matching simhit
      if (wireMatch and stripMatch) {
	if (verboseCSCRecHit2D_) cout << "\t...was matched!" << endl;
	layer_to_cscRecHit2D_[id].push_back(*d);
	chamber_to_cscRecHit2D_[p_id.chamberId().rawId()].push_back(*d);
      }
    }
  }
}


void
CSCRecHitMatcher::matchCSCSegmentsToSimTrack(const CSCSegmentCollection& cscSegments)
{
  if (verboseCSCSegment_) cout << "Matching simtrack to segments" << endl;
  // fetch all chamberIds with simhits
  auto chamber_ids = simhit_matcher_->chamberIdsCSC(0);
  if (verboseCSCSegment_) cout << "Number of matched csc segments " << chamber_ids.size() << endl;
  for (auto id: chamber_ids) {
    CSCDetId p_id(id);
    
    // print all CSCRecHit2D in the CSCChamber
    auto csc_rechits(cscRecHit2DsInChamber(id));
    if (verboseCSCSegment_) {
      cout<<"hit csc rechits" <<endl;
      for (auto rh: csc_rechits) cout << "\t"<< rh << endl;
      cout<<endl;
    }
    
    // get the segments
    auto segments_in_det = cscSegments.get(p_id);
    for (auto d = segments_in_det.first; d != segments_in_det.second; ++d) {
      if (verboseCSCSegment_) cout<<"segment "<<p_id<<" "<<*d<<endl;
      
      //access the rechits
      auto recHits(d->recHits());
      
      int rechitsFound = 0;
      if (verboseCSCSegment_) cout<< recHits.size() << " csc rechits from segment "<<endl;
      for (auto& rh: recHits) {
	const CSCRecHit2D* cscrh(dynamic_cast<const CSCRecHit2D*>(rh));
       	if (isCSCRecHit2DMatched(*cscrh))
       	  ++rechitsFound;
      }
      if (rechitsFound==0) continue;
      if (verboseCSCSegment_) {
	cout << "Found " << rechitsFound << " rechits out of " << cscRecHit2DsInChamber(id).size() << endl;
	cout << "\t...was matched!" << endl;
      }
      chamber_to_cscSegment_[ p_id.rawId() ].push_back(*d);
    }
  }
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


const CSCRecHitMatcher::CSCRecHit2DContainer 
CSCRecHitMatcher::cscRecHit2Ds() const
{
  CSCRecHitMatcher::CSCRecHit2DContainer result;
  for (auto id: chamberIdsCSCRecHit2D()){
    auto segmentsInChamber(cscRecHit2DsInChamber(id));
    result.insert(result.end(), segmentsInChamber.begin(), segmentsInChamber.end());
  }
  return result;
}


const CSCRecHitMatcher::CSCSegmentContainer 
CSCRecHitMatcher::cscSegments() const
{
  CSCRecHitMatcher::CSCSegmentContainer result;
  for (auto id: chamberIdsCSCSegment()){
    auto segmentsInChamber(cscSegmentsInChamber(id));
    result.insert(result.end(), segmentsInChamber.begin(), segmentsInChamber.end());
  }
  return result;
}

  
bool 
CSCRecHitMatcher::cscRecHit2DInContainer(const CSCRecHit2D& sg, const CSCRecHit2DContainer& c) const
{
  bool isSame = false;
  for (auto& segment: c) if (areCSCRecHit2DsSame(sg,segment)) isSame = true;
  return isSame;
}


bool 
CSCRecHitMatcher::cscSegmentInContainer(const CSCSegment& sg, const CSCSegmentContainer& c) const
{
  bool isSame = false;
  for (auto& segment: c) if (areCSCSegmentsSame(sg,segment)) isSame = true;
  return isSame;
}


bool 
CSCRecHitMatcher::isCSCRecHit2DMatched(const CSCRecHit2D& thisSg) const
{
  return cscRecHit2DInContainer(thisSg, cscRecHit2Ds());
}


bool 
CSCRecHitMatcher::isCSCSegmentMatched(const CSCSegment& thisSg) const
{
  return cscSegmentInContainer(thisSg, cscSegments());
}


int 
CSCRecHitMatcher::nCSCRecHits() const
{
  int n = 0;
  auto ids = chamberIdsCSCRecHit2D();
  for (auto id: ids) n += cscRecHit2DsInChamber(id).size();
  return n;  
}


int 
CSCRecHitMatcher::nCSCSegments() const
{
  int n = 0;
  auto ids = chamberIdsCSCSegment();
  for (auto id: ids) n += cscSegmentsInChamber(id).size();
  return n;  
}


bool 
CSCRecHitMatcher::areCSCRecHit2DsSame(const CSCRecHit2D& l,const CSCRecHit2D& r) const
{
  return l.localPosition() == r.localPosition();
}


bool 
CSCRecHitMatcher::areCSCSegmentsSame(const CSCSegment& l,const CSCSegment& r) const
{
  return (l.localPosition() == r.localPosition() and l.localDirection() == r.localDirection());
}


CSCSegment 
CSCRecHitMatcher::bestCSCSegment(unsigned int id)
{
  CSCSegment emptySegment;
  double chi2overNdf = 99;
  int index=0;
  int foundIndex=-99;
  
  for (auto& seg: chamber_to_cscSegment_[id]){
    double newChi2overNdf(seg.chi2()/seg.degreesOfFreedom());
    if (newChi2overNdf < chi2overNdf) {
      chi2overNdf = newChi2overNdf;
      foundIndex = index;
    }
    ++index;
  }
  if (foundIndex != -99) return chamber_to_cscSegment_[id][foundIndex];
  return emptySegment;
}

