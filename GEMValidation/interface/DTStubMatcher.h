#ifndef GEMCode_GEMValidation_DTStubMatcher_h
#define GEMCode_GEMValidation_DTStubMatcher_h

/**\class StubMatcher

 Description: Matching of Stubs for SimTrack in DT

 Original Author:  Sven Dildick
*/

#include "GEMCode/GEMValidation/interface/DigiMatcher.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DTDigi/interface/DTLocalTriggerCollection.h"

#include <vector>
#include <map>
#include <set>

class SimHitMatcher;

class DTStubMatcher : public DigiMatcher
{
public:

  typedef std::vector<DTLocalTrigger> DTLocalTriggerContainer;

  DTStubMatcher(SimHitMatcher& sh);
  
  ~DTStubMatcher();

  // chamber detIds with stubs
  std::set<unsigned int> chamberIds(int dt_type = DT_ALL) const;

  //DT stubs from a particular chamber
  const DTLocalTriggerContainer& dtLocalTriggersInChamber(unsigned int) const;

  int nDTLocalTriggersInChamber(unsigned int) const;

private:

  void matchDTLocalTriggersToSimTrack(const DTLocalTriggerCollection&);

  std::set<unsigned int> selectDetIds(const std::map<unsigned int, DTLocalTriggerContainer>&, int) const;

  bool verbose_;
  bool run_;
  int minBX_, maxBX_;

  std::map<unsigned int, DTLocalTriggerContainer> chamber_to_stubs_;

  DTLocalTriggerContainer no_dt_stubs_;
};

#endif
