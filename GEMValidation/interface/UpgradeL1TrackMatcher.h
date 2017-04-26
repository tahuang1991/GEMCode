#ifndef GEMCode_GEMValidation_UpgradeL1TrackMatcher_h
#define GEMCode_GEMValidation_UpgradeL1TrackMatcher_h

/**\class UpgradeL1TrackMatcher

 Description: Matching of tracks to SimTrack

 Original Author:  "Sven Dildick"
*/

#include "GEMCode/GEMValidation/interface/CSCStubMatcher.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"

class UpgradeL1TrackMatcher : public BaseMatcher
{
 public:
  /// constructor
  UpgradeL1TrackMatcher(CSCStubMatcher&,
                        edm::EDGetTokenT<l1t::EMTFTrackCollection> &emtfTrackInputLabel_
                        );
  /// destructor
  ~UpgradeL1TrackMatcher();

 private:

  void clear();
  void init();

  void matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection&);

  const CSCStubMatcher* csc_stub_matcher_;

  int minBXEMTFTrack_, maxBXEMTFTrack_;
  int verboseEMTFTrack_;
  double deltaREMTFTrack_;

  l1t::EMTFTrackCollection tfTracks_;
};

#endif
