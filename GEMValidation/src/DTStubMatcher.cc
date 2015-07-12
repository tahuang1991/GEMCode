#include "GEMCode/GEMValidation/interface/DTStubMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

using namespace std;

#include "DataFormats/MuonDetId/interface/DTWireId.h"

DTStubMatcher::DTStubMatcher(SimHitMatcher& sh)
: DigiMatcher(sh)
{
  auto dtStub_= conf().getParameter<edm::ParameterSet>("dtStub");
  input_ = dtStub_.getParameter<std::vector<edm::InputTag>>("validInputTags");
  minBX_ = dtStub_.getParameter<int>("minBX");
  maxBX_ = dtStub_.getParameter<int>("maxBX");
  verbose_ = dtStub_.getParameter<int>("verbose");
  run_ = dtStub_.getParameter<bool>("run");

  if (hasDTGeometry_) {
    edm::Handle<DTLocalTriggerCollection> dt_stubs;
    if(gemvalidation::getByLabel(input_, dt_stubs, event())) if (run_) matchDTLocalTriggersToSimTrack(*dt_stubs.product());
  }
}

DTStubMatcher::~DTStubMatcher() {}


void
DTStubMatcher::matchDTLocalTriggersToSimTrack(const DTLocalTriggerCollection& stubs)
{
  auto det_ids = simhit_matcher_->chamberIdsDT();
  for (auto id: det_ids)
  {
    const DTChamberId c_id(id);

    // get stubs in this layer
    auto stubs_in_det = stubs.get(c_id);

    for (auto d = stubs_in_det.first; d != stubs_in_det.second; ++d)
    {
      if (verbose_) cout<<"dt local trigger digi "<<c_id<<" "<<*d<<endl;
      // check that the digi is within BX range
      if (d->bx() < minBX_ || d->bx() > maxBX_) continue;

      if (verbose_) cout<<"oki"<<endl;

      chamber_to_stubs_[c_id].push_back(*d);
    }
  }
}


std::set<unsigned int>
DTStubMatcher::selectDetIds(const std::map<unsigned int, DTStubMatcher::DTLocalTriggerContainer>&stubs, int dt_type) const
{
  std::set<unsigned int> result;
  for (auto& p: stubs)
  {
    auto id = p.first;
    if (dt_type > 0)
    {
      DTChamberId detId(id);
      if (toDTType(detId) != dt_type) continue;
    }
    result.insert(p.first);
  }
  return result;
}


std::set<unsigned int>
DTStubMatcher::chamberIds(int dt_type) const
{
  return selectDetIds(chamber_to_stubs_, dt_type);
}


const DTStubMatcher::DTLocalTriggerContainer&
DTStubMatcher::dtLocalTriggersInChamber(unsigned int detid) const
{
  if (chamber_to_stubs_.find(detid) == chamber_to_stubs_.end()) return no_dt_stubs_;
  return chamber_to_stubs_.at(detid);
}


int
DTStubMatcher::nDTLocalTriggersInChamber(unsigned int detid) const
{
  return dtLocalTriggersInChamber(detid).size();
}


