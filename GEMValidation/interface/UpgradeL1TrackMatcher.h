#ifndef GEMCode_GEMValidation_UpgradeL1TrackMatcher_h
#define GEMCode_GEMValidation_UpgradeL1TrackMatcher_h

/**\class UpgradeL1TrackMatcher

 Description: Matching of tracks to SimTrack

 Original Author:  "Sven Dildick"
*/

#include "GEMCode/GEMValidation/interface/CSCStubMatcher.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"

class UpgradeL1TrackMatcher : public BaseMatcher
{
 public:
  /// constructor
  UpgradeL1TrackMatcher(CSCStubMatcher&,
                        edm::EDGetTokenT<l1t::EMTFTrackCollection> &,
                        edm::EDGetTokenT< BXVector<l1t::RegionalMuonCand> > &);
  /// destructor
  ~UpgradeL1TrackMatcher();

 private:

  void clear();

  void matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection&);
  void matchGMTToSimTrack(const BXVector<l1t::RegionalMuonCand>&);

  const CSCStubMatcher* csc_stub_matcher_;

  int minBXEMTFTrack_, maxBXEMTFTrack_;
  int verboseEMTFTrack_;
  double deltaREMTFTrack_;

  int minBXGMT_, maxBXGMT_;
  int verboseGMT_;
  double deltaRGMT_;

  l1t::EMTFTrackCollection tfTracks_;
  l1t::RegionalMuonCand gmt_;
};

#endif
