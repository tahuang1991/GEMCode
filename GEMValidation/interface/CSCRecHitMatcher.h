#ifndef GEMCode_GEMValidation_CSCRecHitMatcher_h
#define GEMCode_GEMValidation_CSCRecHitMatcher_h

/**\class DigiMatcher

 Description: Matching of rechits and segments for SimTrack in CSC

 Original Author:  Sven Dildick
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "GEMCode/GEMValidation/interface/GenericDigi.h"
#include "GEMCode/GEMValidation/interface/DigiMatcher.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;

class CSCRecHitMatcher : public BaseMatcher
{
public:
  
  typedef std::vector<CSCRecHit2D> CSCRecHit2DContainer;
  typedef std::vector<CSCSegment> CSCSegmentContainer;
  
  CSCRecHitMatcher(SimHitMatcher& sh, edm::ConsumesCollector & iC);
  
  ~CSCRecHitMatcher() {}

  // layer detIds with CSCRecHit2D
  std::set<unsigned int> layerIdsCSCRecHit2D() const;
  // chamber detIds with CSCRecHit2D
  std::set<unsigned int> chamberIdsCSCRecHit2D() const;
  // chamber detIds with CSCSegment
  std::set<unsigned int> chamberIdsCSCSegment() const;

  //CSC rechits from a particular layer or chamber
  const CSCRecHit2DContainer& cscRecHit2DsInLayer(unsigned int) const;
  const CSCRecHit2DContainer& cscRecHit2DsInChamber(unsigned int) const;
  //CSC segments from a particular chamber
  const CSCSegmentContainer& cscSegmentsInChamber(unsigned int) const;

  const CSCSegmentContainer cscSegments() const;
  const CSCRecHit2DContainer cscRecHit2Ds() const;

  bool cscRecHit2DInContainer(const CSCRecHit2D&, const CSCRecHit2DContainer&) const;
  bool cscSegmentInContainer(const CSCSegment&, const CSCSegmentContainer&) const;
  bool isCSCRecHit2DMatched(const CSCRecHit2D&) const;
  bool isCSCSegmentMatched(const CSCSegment&) const;

  int nCSCRecHits() const;
  int nCSCSegments() const;
  bool areCSCSegmentsSame(const CSCSegment&,const CSCSegment&) const;
  bool areCSCRecHit2DsSame(const CSCRecHit2D&,const CSCRecHit2D&) const;

  int nCSCRecHit2DsInLayer(unsigned int) const;
  int nCSCRecHit2DsInChamber(unsigned int) const;
  int nCSCSegmentsInChamber(unsigned int) const;

  CSCSegment bestCSCSegment(unsigned int);

private:

  const SimHitMatcher* simhit_matcher_;

  void matchCSCRecHit2DsToSimTrack(const CSCRecHit2DCollection&);
  void matchCSCSegmentsToSimTrack(const CSCSegmentCollection&);

  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHit2DInput_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentInput_;

  int verboseCSCRecHit2D_;
  bool runCSCRecHit2D_;
  int maxBXCSCRecHit2D_;
  int minBXCSCRecHit2D_;

  int verboseCSCSegment_;
  bool runCSCSegment_;
  int maxBXCSCSegment_;
  int minBXCSCSegment_;

  std::map<unsigned int, CSCRecHit2DContainer> layer_to_cscRecHit2D_;
  std::map<unsigned int, CSCRecHit2DContainer> chamber_to_cscRecHit2D_;
  std::map<unsigned int, CSCSegmentContainer> chamber_to_cscSegment_;

  CSCRecHit2DContainer no_cscRecHit2Ds_;
  CSCSegmentContainer no_cscSegments_;
};

#endif
