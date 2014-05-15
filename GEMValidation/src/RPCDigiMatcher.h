#ifndef GEMValidation_RPCDigiMatcher_h
#define GEMValidation_RPCDigiMatcher_h

/**\class DigiMatcher

 Description: Matching of Digis for SimTrack in GEM

 Original Author:  "Vadim Khotilovich"
*/

#include "DigiMatcher.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>

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
  //const DigiContainer& digisInChamber(unsigned int) const;
  //const DigiContainer& digisInSuperChamber(unsigned int) const;


  // GEM co-pads from a particular partition or superchamber
  //const DigiContainer& coPadsInDetId(unsigned int) const;
  //const DigiContainer& coPadsInSuperChamber(unsigned int) const;

  // #layers with digis from this simtrack
  //int nLayersWithDigisInSuperChamber(unsigned int) const;
  //int nLayersWithPadsInSuperChamber(unsigned int) const;

  /// How many pads in RPC did this simtrack get in total?
  int nStrips() const;

  /// How many coincidence pads in GEM did this simtrack get in total?
  //int nCoPads() const;

  std::set<int> stripsInDetId(unsigned int) const;
  //std::set<int> padNumbersInDetId(unsigned int) const;
  //std::set<int> coPadNumbersInDetId(unsigned int) const;

  // what unique partitions numbers with digis from this simtrack?
  std::set<int> partitionNumbers() const;
 // std::set<int> partitionNumbersWithCoPads() const;

private:

  void init();

  void matchDigisToSimTrack(const RPCDigiCollection& digis);
  //void matchPadsToSimTrack(const GEMCSCPadDigiCollection& pads);
 // void matchCoPadsToSimTrack(const GEMCSCPadDigiCollection& co_pads);

  edm::InputTag rpcDigiInput_;
  //edm::InputTag rpcPadDigiInput_;
  //edm::InputTag rpcCoPadDigiInput_;

  int minBXRPC_, maxBXRPC_;

  int matchDeltaStrip_;

  std::map<unsigned int, DigiContainer> detid_to_digis_;
  //std::map<unsigned int, DigiContainer> chamber_to_digis_;
 // std::map<unsigned int, DigiContainer> superchamber_to_digis_;


  //std::map<unsigned int, DigiContainer> detid_to_copads_;
  //std::map<unsigned int, DigiContainer> chamber_to_copads_;
  //std::map<unsigned int, DigiContainer> superchamber_to_copads_;

  bool verboseDigi_;
  //bool verboseCoPad_;
};

#endif
