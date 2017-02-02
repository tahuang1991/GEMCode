#ifndef GEMCode_GEMValidation_GEMDigiMatcher_h
#define GEMCode_GEMValidation_GEMDigiMatcher_h

/**\class DigiMatcher

 Description: Matching of Digis for SimTrack in GEM

 Original Author:  "Vadim Khotilovich"
*/
#include "GEMCode/GEMValidation/interface/DigiMatcher.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCoPadDigiCollection.h"

#include <vector>
#include <map>
#include <set>
#include <unordered_set>

typedef std::vector<GEMDigi> GEMDigiContainer;
typedef std::vector<GEMPadDigi> GEMPadDigiContainer;
typedef std::vector<GEMCoPadDigi> GEMCoPadDigiContainer;

class SimHitMatcher;

class GEMDigiMatcher : public DigiMatcher
{
public:

  GEMDigiMatcher(SimHitMatcher& sh,
                 edm::EDGetTokenT<GEMDigiCollection>& gemDigiInput_, 
                 edm::EDGetTokenT<GEMPadDigiCollection>& gemPadDigiInput_, 
                 edm::EDGetTokenT<GEMCoPadDigiCollection>& gemCoPadDigiInput_);
  
  ~GEMDigiMatcher();

  // partition GEM detIds with digis
  std::set<unsigned int> detIdsDigi(int gem_type = GEM_ALL) const;
  std::set<unsigned int> detIdsPad(int gem_type = GEM_ALL) const;
  //std::set<unsigned int> detIdsCoPad(int gem_type = GEM_ALL) const;

  // chamber detIds with digis
  std::set<unsigned int> chamberIdsDigi(int gem_type = GEM_ALL) const;
  std::set<unsigned int> chamberIdsPad(int gem_type = GEM_ALL) const;

  // superchamber detIds with digis
  std::set<unsigned int> superChamberIdsDigi(int gem_type = GEM_ALL) const;
  std::set<unsigned int> superChamberIdsPad(int gem_type = GEM_ALL) const;
  std::set<unsigned int> superChamberIdsCoPad(int gem_type = GEM_ALL) const;

  // GEM digis from a particular partition, chamber or superchamber
  const DigiContainer& digisInDetId(unsigned int) const;
  const DigiContainer& digisInChamber(unsigned int) const;
  const DigiContainer& digisInSuperChamber(unsigned int) const;

  // GEM pads from a particular partition, chamber or superchamber
  const DigiContainer& padsInDetId(unsigned int) const;
  const DigiContainer& padsInChamber(unsigned int) const;
  const DigiContainer& padsInSuperChamber(unsigned int) const;
  //const DigiContainer& copadsInDetId(unsigned int) const;

  // GEM co-pads from a particular partition or superchamber
  const DigiContainer& coPadsInSuperChamber(unsigned int) const;

  // GEM digis from a particular partition, chamber or superchamber
  const GEMDigiContainer& gemDigisInDetId(unsigned int) const;
  const GEMDigiContainer& gemDigisInChamber(unsigned int) const;
  const GEMDigiContainer& gemDigisInSuperChamber(unsigned int) const;

  // GEM pads from a particular partition, chamber or superchamber
  const GEMPadDigiContainer& gemPadsInDetId(unsigned int) const;
  const GEMPadDigiContainer& gemPadsInChamber(unsigned int) const;
  const GEMPadDigiContainer& gemPadsInSuperChamber(unsigned int) const;

  // GEM co-pads from a particular partition or superchamber
  const GEMCoPadDigiContainer& gemCoPadsInSuperChamber(unsigned int) const;

  const std::map<unsigned int, GEMPadDigiContainer> allGempadsMatch2SimMuon() const { return detid_to_gempads_; }
  const std::map<unsigned int, GEMPadDigiContainer> allGempadsMatch2SimMuon_2strip() const { return detid_to_gempads_2strip_; }
  // #layers with digis from this simtrack
  int nLayersWithDigisInSuperChamber(unsigned int) const;
  int nLayersWithPadsInSuperChamber(unsigned int) const;

  /// How many pads in GEM did this simtrack get in total?
  int nPads() const;

  /// How many coincidence pads in GEM did this simtrack get in total?
  int nCoPads() const;

  std::set<int> stripNumbersInDetId(unsigned int) const;
  std::set<int> padNumbersInDetId(unsigned int) const;

  int extrapolateHsfromGEMPad(unsigned int , int) const;
  int extrapolateHsfromGEMStrip(unsigned int , int) const;
  // what unique partitions numbers with digis from this simtrack?
  std::set<int> partitionNumbers() const;
  std::set<int> partitionNumbersWithCoPads() const;

  std::vector<GlobalPoint> positionPad1InDetId(unsigned int) const;
  std::vector<GlobalPoint> positionPad2InDetId(unsigned int) const;
  std::vector<GlobalPoint> positionPad4InDetId(unsigned int) const;
  std::vector<GlobalPoint> positionPad8InDetId(unsigned int) const;

  GlobalPoint getGlobalPointDigi(unsigned int rawId, const GEMDigi& d) const;
  GlobalPoint getGlobalPointPad(unsigned int rawId, const GEMPadDigi& tp) const;

private:

  void matchDigisToSimTrack(const GEMDigiCollection&);
  void matchPadsToSimTrack(const GEMPadDigiCollection&);
  void matchCoPadsToSimTrack(const GEMCoPadDigiCollection&);

  std::set<unsigned int> selectDetIds(const Id2DigiContainer &, int) const;
  

  int minBXGEMDigi_, maxBXGEMDigi_;
  int minBXGEMPad_, maxBXGEMPad_;
  int minBXGEMCoPad_, maxBXGEMCoPad_;

  int matchDeltaStrip_;

  std::map<unsigned int, DigiContainer> detid_to_digis_;
  std::map<unsigned int, DigiContainer> chamber_to_digis_;
  std::map<unsigned int, DigiContainer> superchamber_to_digis_;

  std::map<unsigned int, DigiContainer> detid_to_pads_;
  std::map<unsigned int, DigiContainer> chamber_to_pads_;
  std::map<unsigned int, DigiContainer> superchamber_to_pads_;

  std::map<unsigned int, DigiContainer> superchamber_to_copads_;

  std::map<unsigned int, GEMDigiContainer> detid_to_gemdigis_;
  std::map<unsigned int, GEMDigiContainer> chamber_to_gemdigis_;
  std::map<unsigned int, GEMDigiContainer> superchamber_to_gemdigis_;

  std::map<unsigned int, GEMPadDigiContainer> detid_to_gempads_;
  std::map<unsigned int, GEMPadDigiContainer> detid_to_gempads_2strip_;
  std::map<unsigned int, GEMPadDigiContainer> chamber_to_gempads_;
  std::map<unsigned int, GEMPadDigiContainer> superchamber_to_gempads_;

  //std::map<unsigned int, DigiContainer> detid_to_copads_;
  std::map<unsigned int, GEMCoPadDigiContainer> superchamber_to_gemcopads_;

  bool verboseDigi_;
  bool verbosePad_;
  bool verboseCoPad_;

  bool runGEMDigi_;
  bool runGEMPad_;
  bool runGEMCoPad_;

  GEMDigiContainer no_gem_digis_;
  GEMPadDigiContainer no_gem_pads_;
  GEMCoPadDigiContainer no_gem_copads_;

};

#endif
