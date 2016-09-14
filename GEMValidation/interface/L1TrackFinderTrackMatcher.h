#ifndef GEMCode_GEMValidation_L1TrackFinderTrackMatcher_h
#define GEMCode_GEMValidation_L1TrackFinderTrackMatcher_h

#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"

class L1TrackFinderTrackMatcher : public BaseMatcher
{
 public:
  /// constructor
  L1TrackFinderTrackMatcher(SimHitMatcher& sh, 
                            edm::EDGetTokenT<L1CSCTrackCollection>& cscTfTrackInputLabel_);
  /// destructor
  ~L1TrackFinderTrackMatcher();
  
 private:
  
  void clear();
  void init();
  
  void matchCSCTfTrackToSimTrack(const L1CSCTrackCollection&);

  int verboseCscTfTrack_;
  bool runCscTfTrack_;
  int minBXCscTfTrack_;
  int maxBXCscTfTrack_;

  L1CSCTrackCollection matchedTfTracks_;
};

#endif
