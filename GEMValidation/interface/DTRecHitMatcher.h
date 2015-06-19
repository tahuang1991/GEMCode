#ifndef GEMCode_GEMValidation_DTRecHitMatcher_h
#define GEMCode_GEMValidation_DTRecHitMatcher_h

/**\class DigiMatcher

 Description: Matching of rechits and segments for SimTrack in DT

 Original Author:  Sven Dildick
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;

class DTRecHitMatcher : public BaseMatcher
{
public:
  
  typedef std::vector<DTRecSegment2D> DTRecSegment2DContainer;
  typedef std::vector<DTRecSegment4D> DTRecSegment4DContainer;
  
  DTRecHitMatcher(SimHitMatcher& sh);
  
  ~DTRecHitMatcher() {}

  // superlayer detIds with DTRecSegment2D
  std::set<unsigned int> superLayerIdsDTRecSegment2D() const;
  // chamber detIds with DTRecSegment2D
  std::set<unsigned int> chamberIdsDTRecSegment2D() const;
  // chamber detIds with DTRecSegment4D
  std::set<unsigned int> chamberIdsDTRecSegment4D() const;

  //DT segments from a particular superlayer or chamber
  const DTRecSegment2DContainer& dtRecSegment2DInSuperLayer(unsigned int) const;
  const DTRecSegment2DContainer& dtRecSegment2DInChamber(unsigned int) const;
  const DTRecSegment4DContainer& dtRecSegment4DInChamber(unsigned int) const;

  int nDTRecSegment2DInSuperLayer(unsigned int) const;
  int nDTRecSegment2DInChamber(unsigned int) const;
  int nDTRecSegment4DInChamber(unsigned int) const;

private:

  const SimHitMatcher* simhit_matcher_;

  void matchDTRecSegment2DsToSimTrack(const DTRecSegment2DCollection&);
  void matchDTRecSegment4DsToSimTrack(const DTRecSegment4DCollection&);

  std::vector<edm::InputTag> dtRecSegment2DInput_;
  std::vector<edm::InputTag> dtRecSegment4DInput_;

  bool verboseDTRecSegment2D_;
  bool runDTRecSegment2D_;
  int maxBXDTRecSegment2D_;
  int minBXDTRecSegment2D_;

  bool verboseDTRecSegment4D_;
  bool runDTRecSegment4D_;
  int maxBXDTRecSegment4D_;
  int minBXDTRecSegment4D_;

  std::map<unsigned int, DTRecSegment2DContainer> superLayer_to_dtRecSegment2D_;
  std::map<unsigned int, DTRecSegment2DContainer> chamber_to_dtRecSegment2D_;
  std::map<unsigned int, DTRecSegment4DContainer> chamber_to_dtRecSegment4D_;

  DTRecSegment4DContainer no_dtRecSegment2Ds_;
  DTRecSegment4DContainer no_dtRecSegment4Ds_;
};

#endif
