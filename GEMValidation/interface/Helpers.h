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
    const bool verbose(false);
    bool inputTagIsNotValid(true);
    for (unsigned i=0; i<tags.size(); ++i){
      iEvent.getByLabel(tags[i], result);
      if (result.isValid()) {
        if (verbose) std::cout << tags[i] << " is a valid inputTag " << i << std::endl;
        inputTagIsNotValid = false;
        break;
      } else {
        if (verbose) std::cout << tags[i] << " is an invalid inputTag " << i << std::endl;
      }
    }
    return (!inputTagIsNotValid);
  }
}

#endif
