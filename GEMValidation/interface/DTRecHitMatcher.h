#ifndef GEMCode_GEMValidation_DTRecHitMatcher_h
#define GEMCode_GEMValidation_DTRecHitMatcher_h

/**\class DigiMatcher

 Description: Matching of rechits and segments for SimTrack in DT

 Original Author:  Sven Dildick
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;

class DTRecHitMatcher : public BaseMatcher
{
public:
  
  typedef std::vector<DTRecHit1DPair> DTRecHit1DPairContainer;
  typedef std::vector<DTRecSegment2D> DTRecSegment2DContainer;
  typedef std::vector<DTRecSegment4D> DTRecSegment4DContainer;
  
  DTRecHitMatcher(SimHitMatcher& sh, edm::ConsumesCollector & iC);
  
  ~DTRecHitMatcher() {}

  // layer detIds with DTRecHit1DPair
  std::set<unsigned int> layerIdsDTRecHit1DPair() const;
  // superlayer detIds with DTRecHit1DPair
  std::set<unsigned int> superLayerIdsDTRecHit1DPair() const;
  // chamber detIds with DTRecHit1DPair
  std::set<unsigned int> chamberIdsDTRecHit1DPair() const;
  // superlayer detIds with DTRecSegment2D
  std::set<unsigned int> superLayerIdsDTRecSegment2D() const;
  // chamber detIds with DTRecSegment2D
  std::set<unsigned int> chamberIdsDTRecSegment2D() const;
  // chamber detIds with DTRecSegment4D
  std::set<unsigned int> chamberIdsDTRecSegment4D() const;

  //DT rechits from a particular layer, superlayer or chamber
  const DTRecHit1DPairContainer& dtRecHit1DPairInLayer(unsigned int) const;
  const DTRecHit1DPairContainer& dtRecHit1DPairInSuperLayer(unsigned int) const;
  const DTRecHit1DPairContainer& dtRecHit1DPairInChamber(unsigned int) const;
  //DT segments from a particular superlayer or chamber
  const DTRecSegment2DContainer& dtRecSegment2DInSuperLayer(unsigned int) const;
  const DTRecSegment2DContainer& dtRecSegment2DInChamber(unsigned int) const;
  const DTRecSegment4DContainer& dtRecSegment4DInChamber(unsigned int) const;

  const DTRecHit1DPairContainer dtRecHit1DPairs() const;
  const DTRecSegment2DContainer dtRecSegment2Ds() const;
  const DTRecSegment4DContainer dtRecSegment4Ds() const;

  int nDTRecHit1DPairInLayer(unsigned int) const;
  int nDTRecHit1DPairInSuperLayer(unsigned int) const;
  int nDTRecHit1DPairInChamber(unsigned int) const;
  int nDTRecSegment2DInSuperLayer(unsigned int) const;
  int nDTRecSegment2DInChamber(unsigned int) const;
  int nDTRecSegment4DInChamber(unsigned int) const;

  int nDTRecHit1DPairs() const;
  int nDTRecSegment2Ds() const;
  int nDTRecSegment4Ds() const;

  bool dtRecHit1DPairInContainer(const DTRecHit1DPair&, const DTRecHit1DPairContainer&) const;
  bool dtRecSegment2DInContainer(const DTRecSegment2D&, const DTRecSegment2DContainer&) const;
  bool dtRecSegment4DInContainer(const DTRecSegment4D&, const DTRecSegment4DContainer&) const;

  bool isDTRecHit1DPairMatched(const DTRecHit1DPair&) const;
  bool isDTRecSegment2DMatched(const DTRecSegment2D&) const;
  bool isDTRecSegment4DMatched(const DTRecSegment4D&) const;
  bool areDTRecSegment4DSame(const DTRecSegment4D&,const DTRecSegment4D&) const;

private:

  const SimHitMatcher* simhit_matcher_;

  void matchDTRecHit1DPairsToSimTrack(const DTRecHitCollection&);
  void matchDTRecSegment2DsToSimTrack(const DTRecSegment2DCollection&);
  void matchDTRecSegment4DsToSimTrack(const DTRecSegment4DCollection&);

  edm::EDGetTokenT<DTRecHitCollection> dtRecHit1DPairInput_;
  edm::EDGetTokenT<DTRecSegment2DCollection> dtRecSegment2DInput_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtRecSegment4DInput_;

  bool verboseDTRecHit1DPair_;
  bool runDTRecHit1DPair_;
  int maxBXDTRecHit1DPair_;
  int minBXDTRecHit1DPair_;

  bool verboseDTRecSegment2D_;
  bool runDTRecSegment2D_;
  int maxBXDTRecSegment2D_;
  int minBXDTRecSegment2D_;

  bool verboseDTRecSegment4D_;
  bool runDTRecSegment4D_;
  int maxBXDTRecSegment4D_;
  int minBXDTRecSegment4D_;

  std::map<unsigned int, DTRecHit1DPairContainer> layer_to_dtRecHit1DPair_;
  std::map<unsigned int, DTRecHit1DPairContainer> superLayer_to_dtRecHit1DPair_;
  std::map<unsigned int, DTRecHit1DPairContainer> chamber_to_dtRecHit1DPair_;
  std::map<unsigned int, DTRecSegment2DContainer> superLayer_to_dtRecSegment2D_;
  std::map<unsigned int, DTRecSegment2DContainer> chamber_to_dtRecSegment2D_;
  std::map<unsigned int, DTRecSegment4DContainer> chamber_to_dtRecSegment4D_;

  DTRecHit1DPairContainer no_dtRecHit1DPairs_;
  DTRecSegment2DContainer no_dtRecSegment2Ds_;
  DTRecSegment4DContainer no_dtRecSegment4Ds_;
};

#endif
