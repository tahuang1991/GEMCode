#ifndef GEMCode_GEMValidation_Helpers_h
#define GEMCode_GEMValidation_Helpers_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

namespace gemvalidation
{
  template<typename PROD>
  bool
  getByLabel(std::vector<edm::InputTag> const& tags, edm::Handle<PROD>& result, const edm::Event& iEvent)
  {
    int i=0;
    const bool verbose(false);
    bool inputTagIsNotValid(true);
    while(inputTagIsNotValid) {
      iEvent.getByLabel(tags.at(i), result);
      if (result.isValid()) {
        if (verbose) std::cout << tags.at(i) << " is a valid inputTag" << std::endl;
        inputTagIsNotValid = false;
      } else {
      if (verbose) std::cout << tags.at(i) << " is an invalid inputTag" << std::endl;
      }
      ++i;
    }
    return (!inputTagIsNotValid);
  }
}

#endif
