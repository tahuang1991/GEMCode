#ifndef GEMCode_GEMValidation_DTDigiMatcher_h
#define GEMCode_GEMValidation_DTDigiMatcher_h

/**\class DigiMatcher

 Description: Matching of Digis for SimTrack in GEM

 Original Author:  "Vadim Khotilovich"
*/

#include "DigiMatcher.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"

#include <vector>
#include <map>
#include <set>

typedef std::vector<DTDigi> DTDigiContainer;

class SimHitMatcher;

class DTDigiMatcher : public DigiMatcher
{
public:

  DTDigiMatcher(SimHitMatcher& sh);
  
  ~DTDigiMatcher();

  // partition GEM detIds with digis
  std::set<unsigned int> detIds() const;

  // chamber detIds with digis
  std::set<unsigned int> chamberIds() const;

  // chamber detIds with digis
  std::set<unsigned int> layerIds() const;

  // superchamber detIds with digis
  std::set<unsigned int> superLayerIds() const;

  //DT digis from a particular partition, chamber or superchamber
  const DTDigiContainer& digisInDetId(unsigned int) const;
  const DTDigiContainer& digisInLayer(unsigned int) const;
  const DTDigiContainer& digisInSuperLayer(unsigned int) const;
  const DTDigiContainer& digisInChamber(unsigned int) const;

  // #layers with digis from this simtrack
  int nLayersWithDigisInSuperChamber(unsigned int) const;

  // what unique partitions numbers with digis from this simtrack?
  std::set<int> partitionNumbers() const;

private:

  void init();

  void matchWireDigisToSimTrack(const DTDigiCollection& digis);

  edm::InputTag dtDigiInput_;

  int minBXDT_, maxBXDT_;

  int matchDeltaWire_;

  std::map<unsigned int, DTDigiContainer> detid_to_digis_;
  std::map<unsigned int, DTDigiContainer> layer_to_digis_;
  std::map<unsigned int, DTDigiContainer> superLayer_to_digis_;
  std::map<unsigned int, DTDigiContainer> chamber_to_digis_;

  bool verboseDigi_;

  bool runDTDigi_;
};

#endif
