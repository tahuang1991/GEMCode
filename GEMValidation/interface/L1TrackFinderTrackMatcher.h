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
  void matchDTTfTrackToSimTrack(const L1CSCTrackCollection&);
  void matchRPCTfTrackToSimTrack(const L1CSCTrackCollection&);

  std::vector<edm::InputTag> dtTfTrackInputLabel_;
  std::vector<edm::InputTag> rpcTfTrackInputLabel_;

  int verboseCscTfTrack_;
  int verboseDtTfTrack_;
  int verboseRpcTfTrack_;

  bool runCscTfTrack_;
  bool runDtTfTrack_;
  bool runRpcTfTrack_;

  int minBXCscTfTrack_;
  int minBXDtTfTrack_;
  int minBXRpcTfTrack_;

  int maxBXCscTfTrack_;
  int maxBXDtTfTrack_;
  int maxBXRpcTfTrack_;

  L1CSCTrackCollection matchedTfTracks_;
};

#endif
