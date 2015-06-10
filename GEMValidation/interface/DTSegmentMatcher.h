#ifndef GEMCode_GEMValidation_DTSegmentMatcher_h
#define GEMCode_GEMValidation_DTSegmentMatcher_h

/*

#include "GEMCode/GEMValidation/interface/DigiMatcher.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DTSegment/interface/DTSegmentCollection.h"

#include <vector>
#include <map>
#include <set>

typedef std::vector<DTSegment> DTSegmentContainer;

class SimHitMatcher;

class DTSegmentMatcher : public DigiMatcher
{
public:

  DTSegmentMatcher(SimHitMatcher& sh);
  
  ~DTSegmentMatcher();

  // partition GEM detIds with digis
  std::set<unsigned int> detIds() const;
  // chamber detIds with digis
  std::set<unsigned int> layerIds() const;
  // superchamber detIds with digis
  std::set<unsigned int> superLayerIds() const;
  // chamber detIds with digis
  std::set<unsigned int> chamberIds() const;

  //DT digis from a particular partition, chamber or superchamber
  const DTSegmentContainer& digisInDetId(unsigned int) const;
  const DTSegmentContainer& digisInLayer(unsigned int) const;
  const DTSegmentContainer& digisInSuperLayer(unsigned int) const;
  const DTSegmentContainer& digisInChamber(unsigned int) const;

  // #tubes with digis in layer from this simtrack
  int nTubesWithDigisInLayer(unsigned int) const;
  // #layers with digis from this simtrack
  int nLayersWithDigisInSuperLayer(unsigned int) const;
  // #layers with digis from this simtrack
  int nSuperLayersWithDigisInChamber(unsigned int) const;

  // wire numbers from this simtrack in a detId
  std::set<int> wireNumbersInDetId(unsigned int detid) const;

private:

  void matchDigisToSimTrack(const DTSegmentCollection& digis);

  std::vector<edm::InputTag> dtDigiInput_;

  bool verboseDigi_;
  bool runDTSegment_;
  int minBXDT_, maxBXDT_;
  int matchDeltaWire_;

  std::map<unsigned int, DTSegmentContainer> detid_to_digis_;
  std::map<unsigned int, DTSegmentContainer> layer_to_digis_;
  std::map<unsigned int, DTSegmentContainer> superLayer_to_digis_;
  std::map<unsigned int, DTSegmentContainer> chamber_to_digis_;

  DTSegmentContainer no_dt_digis_;
};

*/

#endif
