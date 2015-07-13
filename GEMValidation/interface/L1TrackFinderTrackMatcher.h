#ifndef GEMCode_GEMValidation_L1TrackFinderTrackMatcher_h
#define GEMCode_GEMValidation_L1TrackFinderTrackMatcher_h

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"

class L1TrackFinderTrackMatcher// : public BaseMatcher
{
 public:
  /// constructor
  L1TrackFinderTrackMatcher();
  /// destructor
  ~L1TrackFinderTrackMatcher();
  
 private:
  
  void clear();
  void init();
  
  void matchCSCTfTrackToSimTrack(const L1CSCTrackCollection&);
  void matchDTTfTrackToSimTrack(const L1CSCTrackCollection&);
  void matchRPCTfTrackToSimTrack(const L1CSCTrackCollection&);

  std::vector<edm::InputTag> cscTfTrackInputLabel_;
  std::vector<edm::InputTag> dtTfTrackInputLabel_;
  std::vector<edm::InputTag> rpcTfTrackInputLabel_;
};

#endif
