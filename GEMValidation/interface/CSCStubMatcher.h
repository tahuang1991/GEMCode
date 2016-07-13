#ifndef GEMCode_GEMValidation_CSCStubMatcher_h
#define GEMCode_GEMValidation_CSCStubMatcher_h

/**\class CSCStubMatcher

 Description: Matching of CSC L1 trigger stubs to SimTrack

 Original Author:  "Vadim Khotilovich"
 $Id: $
*/

#include "GEMCode/GEMValidation/interface/CSCDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/GEMDigiMatcher.h"
#include "GEMCode/GEMValidation/interface/RPCDigiMatcher.h"

#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include <vector>
#include <map>
#include <set>

typedef std::vector<CSCALCTDigi> CSCALCTDigiContainer;
typedef std::vector<CSCCLCTDigi> CSCCLCTDigiContainer;
typedef std::vector<CSCCorrelatedLCTDigi> CSCCorrelatedLCTDigiContainer;

class SimHitMatcher;

class CSCStubMatcher : public DigiMatcher
{
public:
  
  typedef std::map<int, std::vector<std::pair<unsigned int, const GEMCSCPadDigi*> > > GEMPads;
  typedef std::pair<unsigned int, const GEMCSCPadDigi*> GEMPadBX;
  typedef std::vector<GEMPadBX> GEMPadsBX;

  CSCStubMatcher(SimHitMatcher& sh, CSCDigiMatcher& dg, GEMDigiMatcher& gem_dg, RPCDigiMatcher& rpc_dg);
  
  ~CSCStubMatcher();

  /// crossed chamber detIds with not necessarily matching stubs
  std::set<unsigned int> chamberIdsAllCLCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsAllALCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsAllLCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsAllMPLCT(int csc_type = CSC_ALL) const;
 
  /// chamber detIds with matching stubs
  /// by default, only returns those from ME1b; use al chambers if csc_type=0
  std::set<unsigned int> chamberIdsCLCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsALCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsLCT(int csc_type = CSC_ALL) const;
  std::set<unsigned int> chamberIdsMPLCT(int csc_type = CSC_ALL) const;

  /// single matched stubs from a particular chamber
  Digi clctInChamber(unsigned int) const;
  Digi alctInChamber(unsigned int) const;
  Digi lctInChamber(unsigned int) const;
  Digi mplctInChamber(unsigned int) const;

  /// all stubs (not necessarily matching) from a particular crossed chamber
  const DigiContainer& allCLCTsInChamber(unsigned int) const;
  const DigiContainer& allALCTsInChamber(unsigned int) const;
  const DigiContainer& allLCTsInChamber(unsigned int) const;
  const DigiContainer& allMPLCTsInChamber(unsigned int) const;

  /// all matching from a particular crossed chamber
  const DigiContainer& clctsInChamber(unsigned int) const;
  const DigiContainer& alctsInChamber(unsigned int) const;
  const DigiContainer& lctsInChamber(unsigned int) const;
  const DigiContainer& mplctsInChamber(unsigned int) const;

  /// all stubs (not necessarily matching) from a particular crossed chamber
  const CSCCLCTDigiContainer& allCscCLCTsInChamber(unsigned int) const;
  const CSCALCTDigiContainer& allCscALCTsInChamber(unsigned int) const;
  const CSCCorrelatedLCTDigiContainer& allCscLCTsInChamber(unsigned int) const;
  const CSCCorrelatedLCTDigiContainer& allCscMPLCTsInChamber(unsigned int) const;

  /// all matching from a particular crossed chamber
  const CSCCLCTDigiContainer& cscClctsInChamber(unsigned int) const;
  const CSCALCTDigiContainer& cscAlctsInChamber(unsigned int) const;
  const CSCCorrelatedLCTDigiContainer& cscLctsInChamber(unsigned int) const;
  const CSCCorrelatedLCTDigiContainer& cscMplctsInChamber(unsigned int) const;

  const DigiContainer lctsInStation(int) const;
  /// How many CSC chambers with matching stubs of some minimal quality did this SimTrack hit?
  int nChambersWithCLCT(int min_quality = 0) const;
  int nChambersWithALCT(int min_quality = 0) const;
  int nChambersWithLCT(int min_quality = 0) const;
  int nChambersWithMPLCT(int min_quality = 0) const;

  bool checkStubInChamber(CSCDetId id, CSCCorrelatedLCTDigi lct) const;

private:

  void matchCLCTsToSimTrack(const CSCCLCTDigiCollection&);
  void matchALCTsToSimTrack(const CSCALCTDigiCollection&);
  void matchLCTsToSimTrack(const CSCCorrelatedLCTDigiCollection&);
  void matchMPLCTsToSimTrack(const CSCCorrelatedLCTDigiCollection&);

  const CSCDigiMatcher* digi_matcher_;
  const GEMDigiMatcher* gem_digi_matcher_;
  const RPCDigiMatcher* rpc_digi_matcher_;
  const SimHitMatcher* sh_matcher_;

  // matched stubs in crossed chambers
  Id2Digi chamber_to_clct_;
  Id2Digi chamber_to_alct_;
  Id2Digi chamber_to_lct_;
  Id2Digi chamber_to_mplct_;

  // all stubs (not necessarily matching) in crossed chambers with digis
  Id2DigiContainer chamber_to_clcts_all_;
  Id2DigiContainer chamber_to_alcts_all_;
  Id2DigiContainer chamber_to_lcts_all_;
  Id2DigiContainer chamber_to_mplcts_all_;

  // all matching stubs in crossed chambers with digis
  Id2DigiContainer chamber_to_clcts_;
  Id2DigiContainer chamber_to_alcts_;
  Id2DigiContainer chamber_to_lcts_;
  Id2DigiContainer chamber_to_mplcts_;

  // all stubs (not necessarily matching) in crossed chambers with digis
  std::map<unsigned int, CSCCLCTDigiContainer> chamber_to_cscClcts_all_;
  std::map<unsigned int, CSCALCTDigiContainer> chamber_to_cscAlcts_all_;
  std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamber_to_cscLcts_all_;
  std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamber_to_cscMplcts_all_;

  // all matching stubs in crossed chambers with digis
  std::map<unsigned int, CSCCLCTDigiContainer> chamber_to_cscClcts_;
  std::map<unsigned int, CSCALCTDigiContainer> chamber_to_cscAlcts_;
  std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamber_to_cscLcts_;
  std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamber_to_cscMplcts_;

  template<class D>
  std::set<unsigned int> selectDetIds(D &, int) const;

  bool addGhostLCTs_;
  bool addGhostMPLCTs_;
  bool matchAlctGemME11_;
  bool matchAlctGemME21_;
  bool matchClctGemME11_;
  bool matchClctGemME21_;
  bool matchAlctRpc_;
  bool matchClctRpc_;
  bool hsFromSimHitMean_;

  int minNHitsChamber_;
  int minNHitsChamberALCT_;
  int minNHitsChamberCLCT_;
  int minNHitsChamberLCT_;
  int minNHitsChamberMPLCT_;
  
  bool verboseALCT_;
  bool verboseCLCT_;
  bool verboseLCT_;
  bool verboseMPLCT_;

  bool runALCT_;
  bool runCLCT_;
  bool runLCT_;
  bool runMPLCT_;

  int minBXCLCT_, maxBXCLCT_;
  int minBXALCT_, maxBXALCT_;
  int minBXLCT_, maxBXLCT_;
  int minBXMPLCT_, maxBXMPLCT_;

  std::vector<edm::InputTag> clctInputs_;
  std::vector<edm::InputTag> alctInputs_;
  std::vector<edm::InputTag> lctInputs_;
  std::vector<edm::InputTag> mplctInputs_;

  CSCCLCTDigiContainer no_csc_clcts_;
  CSCALCTDigiContainer no_csc_alcts_;
  CSCCorrelatedLCTDigiContainer no_csc_lcts_;
  CSCCorrelatedLCTDigiContainer no_csc_mplcts_;
};


template<class D>
std::set<unsigned int>
CSCStubMatcher::selectDetIds(D &digis, int csc_type) const
{
  std::set<unsigned int> result;
  for (auto& p: digis)
  {
    auto id = p.first;
    if (csc_type > 0)
    {
      CSCDetId detId(id);
      if (gemvalidation::toCSCType(detId.station(), detId.ring()) != csc_type) continue;
    }
    result.insert(p.first);
  }
  return result;
}

#endif
