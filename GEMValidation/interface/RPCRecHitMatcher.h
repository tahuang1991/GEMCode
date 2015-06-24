#ifndef GEMCode_GEMValidation_RPCRecHitMatcher_h
#define GEMCode_GEMValidation_RPCRecHitMatcher_h

/**\class RecHitMatcher

 Description: Matching of RecHits for SimTrack in GEM

 Original Author:  "Vadim Khotilovich"
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "GEMCode/GEMValidation/interface/GenericDigi.h"
#include "GEMCode/GEMValidation/interface/DigiMatcher.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;
class RPCGeometry;

class RPCRecHitMatcher : BaseMatcher
{
public:

  typedef matching::Digi RecHit;
  typedef matching::DigiContainer RecHitContainer;

  RPCRecHitMatcher(SimHitMatcher& sh);
  
  ~RPCRecHitMatcher();

  // partition RPC detIds with rechits
  std::set<unsigned int> detIds() const;

  // chamber detIds with rechits
  std::set<unsigned int> chamberIds() const;

  // RPC recHits from a particular partition, chamber or superchamber
  const RecHitContainer& recHitsInDetId(unsigned int) const;
  const RecHitContainer& recHitsInChamber(unsigned int) const;

  /// How many recHits in RPC did this simtrack get in total?
  int nRecHits() const;

  std::set<int> stripNumbersInDetId(unsigned int) const;

  // what unique partitions numbers with recHits from this simtrack?
  std::set<int> partitionNumbers() const;

  GlobalPoint recHitPosition(const RecHit& rechit) const;
  GlobalPoint recHitMeanPosition(const RecHitContainer& rechits) const;

private:

  void matchRecHitsToSimTrack(const RPCRecHitCollection& recHits);

  std::vector<edm::InputTag> rpcRecHitInput_;

  const SimHitMatcher* simhit_matcher_;

  int minBXRPC_, maxBXRPC_;

  int matchDeltaStrip_;

  std::map<unsigned int, RecHitContainer> detid_to_recHits_;
  std::map<unsigned int, RecHitContainer> chamber_to_recHits_;

  const RecHitContainer no_recHits_;

  bool verboseRPCRecHit_;
  bool runRPCRecHit_;
};

#endif
