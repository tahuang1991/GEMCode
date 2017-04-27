#ifndef GEMCode_GEMValidation_UpgradeL1TrackMatcher_h
#define GEMCode_GEMValidation_UpgradeL1TrackMatcher_h

/**\class UpgradeL1TrackMatcher

 Description: Matching of tracks to SimTrack

 Original Author:  "Sven Dildick"
*/
#include "GEMCode/GEMValidation/interface/TFTrack.h" 
#include "GEMCode/GEMValidation/interface/TFCand.h" 

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

  std::vector<TFTrack*> tfTracks() const  {return tfTracks_;}
  TFTrack* bestTFTrack() const { return bestTrack; }
  TFCand* bestGMTCand() const { return bestGMT; }
  std::vector<TFCand*> gmts() const { return gmts_; }
  

 private:

  void clear();

  float simPt;
  float simEta;
  float simPhi;
  float simE;
  float simCharge;

  void matchEmtfTrackToSimTrack(const l1t::EMTFTrackCollection&);
  void matchGMTToSimTrack(const BXVector<l1t::RegionalMuonCand>&);

  float mindREMTFTrack = 10;
  TFTrack* bestTrack;

  float mindRGMT = 10;
  TFCand* bestGMT;

  const CSCStubMatcher* csc_stub_matcher_;

  int minBXEMTFTrack_, maxBXEMTFTrack_;
  int verboseEMTFTrack_;
  double deltaREMTFTrack_;

  int minBXGMT_, maxBXGMT_;
  int verboseGMT_;
  double deltaRGMT_;

  //l1t::EMTFTrackCollection tfTracks_;
  std::vector<TFTrack*> tfTracks_;
  std::vector<TFCand*> gmts_;
};

#endif
