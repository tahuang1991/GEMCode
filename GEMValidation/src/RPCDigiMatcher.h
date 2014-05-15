#ifndef RPCValidation_RPCDigiMatcher_h
#define RPCValidation_RPCDigiMatcher_h

/**\class DigiMatcher

 Description: Matching of Digis for SimTrack in RPC

*/

#include "DigiMatcher.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;

class RPCDigiMatcher : public DigiMatcher
{
public:

  RPCDigiMatcher(SimHitMatcher& sh);
  
  ~RPCDigiMatcher();

  // partition RPC detIds with digis
  std::set<unsigned int> detIds() const;

  // chamber detIds with digis
  std::set<unsigned int> chamberIds() const;

  // RPC digis from a particular partition, chamber or superchamber
  const DigiContainer& digisInDetId(unsigned int) const;
  const DigiContainer& digisInChamber(unsigned int) const;

  // #layers with digis from this simtrack
  int nLayersWithDigisInSuperChamber(unsigned int) const;

  std::set<int> stripNumbersInDetId(unsigned int) const;

  // what unique partitions numbers with digis from this simtrack?
  std::set<int> partitionNumbers() const;

private:

  void init();

  void matchDigisToSimTrack(const RPCDigiCollection& digis);

  edm::InputTag rpcDigiInput_;

  int minBXRPC_, maxBXRPC_;

  int matchDeltaStrip_;

  std::map<unsigned int, DigiContainer> detid_to_digis_;
  std::map<unsigned int, DigiContainer> chamber_to_digis_;

  bool verboseDigi_;
};

#endif
