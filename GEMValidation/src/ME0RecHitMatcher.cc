#include "GEMCode/GEMValidation/interface/ME0RecHitMatcher.h"
#include "GEMCode/GEMValidation/interface/ME0DigiMatcher.h"

using namespace std;


ME0RecHitMatcher::ME0RecHitMatcher(ME0DigiMatcher& dg,
                                   edm::EDGetTokenT<ME0RecHitCollection>& me0RecHitInput_,
                                   edm::EDGetTokenT<ME0SegmentCollection>& me0SegmentInput_)
  : BaseMatcher(dg.trk(), dg.vtx(), dg.conf(), dg.event(), dg.eventSetup())
  , digi_matcher_(&dg)
{
  auto me0RecHit = conf().getParameter<edm::ParameterSet>("me0RecHit");
  maxBXME0RecHit_ = me0RecHit.getParameter<int>("maxBX");
  minBXME0RecHit_ = me0RecHit.getParameter<int>("minBX");
  verboseME0RecHit_ = me0RecHit.getParameter<int>("verbose");
  runME0RecHit_ = me0RecHit.getParameter<bool>("run");

  auto me0Segment = conf().getParameter<edm::ParameterSet>("me0Segment");
  maxBXME0Segment_ = me0Segment.getParameter<int>("maxBX");
  minBXME0Segment_ = me0Segment.getParameter<int>("minBX");
  verboseME0Segment_ = me0Segment.getParameter<int>("verbose");
  runME0Segment_ = me0Segment.getParameter<bool>("run");
  minNHitsSegment_ = me0Segment.getParameter<int>("minNHits");

  if (hasME0Geometry_) {
    edm::Handle<ME0RecHitCollection> me0_rechits;
    if (gemvalidation::getByToken(me0RecHitInput_, me0_rechits, event())) if (runME0RecHit_) matchME0RecHitsToSimTrack(*me0_rechits.product());

    edm::Handle<ME0SegmentCollection> me0_segments;
    if (gemvalidation::getByToken(me0SegmentInput_, me0_segments, event())) if (runME0Segment_) matchME0SegmentsToSimTrack(*me0_segments.product());
  }
}


void
ME0RecHitMatcher::matchME0RecHitsToSimTrack(const ME0RecHitCollection& rechits)
{
  if (verboseME0RecHit_) cout << "Matching simtrack to ME0 rechits" << endl;
  // fetch all detIds with digis
  auto layer_ids = digi_matcher_->detIds();
  if (verboseME0RecHit_) cout << "Number of matched me0 layer_ids " << layer_ids.size() << endl;

  for (auto id: layer_ids) {
    ME0DetId p_id(id);

    // get the rechits
    auto rechits_in_det = rechits.get(p_id);
    for (auto rr = rechits_in_det.first; rr != rechits_in_det.second; ++rr) {

      // check that the rechit is within BX range
      if (rr->tof() < minBXME0RecHit_ || rr->tof() > maxBXME0RecHit_) continue;

      if (verboseME0RecHit_) cout<<"rechit "<<p_id<<" "<<*rr << endl;;

      // match the rechit to the digis if the TOF, x and y are the same
      bool match = false;
      ME0DigiPreRecoContainer digis = digi_matcher_->digisInDetId(id);
      for (const auto& digi: digis){
        if (std::abs(digi.x() - rr->localPosition().x())<0.001 and
            std::abs(digi.y() - rr->localPosition().y())<0.001 ) {
          match = true;
        }
      }

      // this rechit was matched to a matching digi
      if (match) {
        if (verboseME0RecHit_) cout << "\t...was matched!" << endl;
        chamber_to_me0RecHit_[id].push_back(*rr);
        superChamber_to_me0RecHit_[p_id.chamberId().rawId()].push_back(*rr);
      }
    }
  }
}


void
ME0RecHitMatcher::matchME0SegmentsToSimTrack(const ME0SegmentCollection& me0Segments)
{
  if (verboseME0Segment_) cout << "Matching simtrack to segments" << endl;
  // fetch all chamberIds with digis
  auto chamber_ids = digi_matcher_->superChamberIds();
  if (verboseME0Segment_) cout << "Number of matched me0 segments " << chamber_ids.size() << endl;
  for (auto id: chamber_ids) {
    ME0DetId p_id(id);

    // print all ME0RecHit in the ME0SuperChamber
    auto me0_rechits(me0RecHitsInSuperChamber(id));
    if (verboseME0Segment_) {
      cout<<"hit me0 rechits" <<endl;
      for (auto rh: me0_rechits) cout << "\t"<< rh << endl;
      cout<<endl;
    }

    // get the segments
    auto segments_in_det = me0Segments.get(p_id);
    for (auto d = segments_in_det.first; d != segments_in_det.second; ++d) {
      if (verboseME0Segment_) cout<<"segment "<<p_id<<" "<<*d<<endl;

      //access the rechits
      auto recHits(d->recHits());

      int rechitsFound = 0;
      if (verboseME0Segment_) cout << "\t has " << recHits.size() << " me0 rechits"<<endl;
      for (auto& rh: recHits) {
        const ME0RecHit* me0rh(dynamic_cast<const ME0RecHit*>(rh));
        if (verboseME0Segment_) cout << "Candidate rechit " << *me0rh << endl;
       	if (isME0RecHitMatched(*me0rh)) {
          if (verboseME0Segment_) cout << "\t...was matched earlier to SimTrack!" << endl;
       	  ++rechitsFound;
        }
      }
      if (rechitsFound<minNHitsSegment_) continue;
      if (verboseME0Segment_) {
        cout << "Found " << rechitsFound << " rechits out of " << me0RecHitsInSuperChamber(id).size() << endl;
        cout << "\t...was matched!" << endl;
      }
      superChamber_to_me0Segment_[ p_id.rawId() ].push_back(*d);
    }
  }
  for (auto& p : superChamber_to_me0Segment_)
    superChamber_to_bestME0Segment_[ p.first] = findbestME0Segment(p.second);
}


std::set<unsigned int>
ME0RecHitMatcher::chamberIdsME0RecHit() const
{
  std::set<unsigned int> result;
  for (auto& p: chamber_to_me0RecHit_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
ME0RecHitMatcher::superChamberIdsME0RecHit() const
{
  std::set<unsigned int> result;
  for (auto& p: superChamber_to_me0RecHit_) result.insert(p.first);
  return result;
}


std::set<unsigned int>
ME0RecHitMatcher::superChamberIdsME0Segment() const
{
  std::set<unsigned int> result;
  for (auto& p: superChamber_to_me0Segment_) result.insert(p.first);
  return result;
}


const ME0RecHitMatcher::ME0RecHitContainer&
ME0RecHitMatcher::me0RecHitsInChamber(unsigned int detid) const
{
  if (chamber_to_me0RecHit_.find(detid) == chamber_to_me0RecHit_.end()) return no_me0RecHits_;
  return chamber_to_me0RecHit_.at(detid);
}


const ME0RecHitMatcher::ME0RecHitContainer&
ME0RecHitMatcher::me0RecHitsInSuperChamber(unsigned int detid) const
{
  if (superChamber_to_me0RecHit_.find(detid) == superChamber_to_me0RecHit_.end()) return no_me0RecHits_;
  return superChamber_to_me0RecHit_.at(detid);
}


const ME0RecHitMatcher::ME0SegmentContainer&
ME0RecHitMatcher::me0SegmentsInSuperChamber(unsigned int detid) const
{
  if (superChamber_to_me0Segment_.find(detid) == superChamber_to_me0Segment_.end()) return no_me0Segments_;
  return superChamber_to_me0Segment_.at(detid);
}


int
ME0RecHitMatcher::nME0RecHitsInChamber(unsigned int detid) const
{
  return me0RecHitsInChamber(detid).size();
}


int
ME0RecHitMatcher::nME0RecHitsInSuperChamber(unsigned int detid) const
{
  return me0RecHitsInSuperChamber(detid).size();
}


int
ME0RecHitMatcher::nME0SegmentsInSuperChamber(unsigned int detid) const
{
  return me0SegmentsInSuperChamber(detid).size();
}


const ME0RecHitMatcher::ME0RecHitContainer
ME0RecHitMatcher::me0RecHits() const
{
  ME0RecHitMatcher::ME0RecHitContainer result;
  for (auto id: superChamberIdsME0RecHit()){
    auto segmentsInSuperChamber(me0RecHitsInSuperChamber(id));
    result.insert(result.end(), segmentsInSuperChamber.begin(), segmentsInSuperChamber.end());
  }
  return result;
}


const ME0RecHitMatcher::ME0SegmentContainer
ME0RecHitMatcher::me0Segments() const
{
  ME0RecHitMatcher::ME0SegmentContainer result;
  for (auto id: superChamberIdsME0Segment()){
    auto segmentsInSuperChamber(me0SegmentsInSuperChamber(id));
    result.insert(result.end(), segmentsInSuperChamber.begin(), segmentsInSuperChamber.end());
  }
  return result;
}


bool
ME0RecHitMatcher::me0RecHitInContainer(const ME0RecHit& rh, const ME0RecHitContainer& c) const
{
  bool isSame = false;
  for (auto& rechit: c) if (areME0RecHitsSame(rh,rechit)) isSame = true;
  return isSame;
}


bool
ME0RecHitMatcher::me0SegmentInContainer(const ME0Segment& sg, const ME0SegmentContainer& c) const
{
  bool isSame = false;
  for (auto& segment: c) if (areME0SegmentsSame(sg,segment)) isSame = true;
  return isSame;
}


bool
ME0RecHitMatcher::isME0RecHitMatched(const ME0RecHit& thisSg) const
{
  return me0RecHitInContainer(thisSg, me0RecHits());
}


bool
ME0RecHitMatcher::isME0SegmentMatched(const ME0Segment& thisSg) const
{
  return me0SegmentInContainer(thisSg, me0Segments());
}


int
ME0RecHitMatcher::nME0RecHits() const
{
  int n = 0;
  auto ids = superChamberIdsME0RecHit();
  for (auto id: ids) n += me0RecHitsInSuperChamber(id).size();
  return n;
}


int
ME0RecHitMatcher::nME0Segments() const
{
  int n = 0;
  auto ids = superChamberIdsME0Segment();
  for (auto id: ids) n += me0SegmentsInSuperChamber(id).size();
  return n;
}


bool
ME0RecHitMatcher::areME0RecHitsSame(const ME0RecHit& l,const ME0RecHit& r) const
{
  return ( l.localPosition() == r.localPosition() and
           l.localPositionError().xx() == r.localPositionError().xx() and
           l.localPositionError().xy() == r.localPositionError().xy() and
           l.localPositionError().yy() == r.localPositionError().yy() and
           l.tof() == r.tof() and
           l.me0Id() == r.me0Id() );
}


bool
ME0RecHitMatcher::areME0SegmentsSame(const ME0Segment& l,const ME0Segment& r) const
{
  return (l.localPosition() == r.localPosition() and l.localDirection() == r.localDirection());
}


ME0Segment
ME0RecHitMatcher::bestME0Segment(unsigned int id)
{
  return superChamber_to_bestME0Segment_[id];
}

ME0Segment
ME0RecHitMatcher::findbestME0Segment(ME0SegmentContainer allSegs) const
{
  ME0Segment bestSegment;
  double chi2overNdf = 99;

  for (auto& seg: allSegs){
    double newChi2overNdf(seg.chi2()/seg.degreesOfFreedom());
    if (newChi2overNdf < chi2overNdf) {
      chi2overNdf = newChi2overNdf;
      bestSegment = seg;
    }
  }
  return bestSegment;
}


GlobalPoint
ME0RecHitMatcher::globalPoint(const ME0Segment& c) const
{
  return getME0Geometry()->idToDet(c.me0DetId())->surface().toGlobal(c.localPosition());
}
