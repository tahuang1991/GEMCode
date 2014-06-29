#ifndef GEMValidation_TrackMatcher_h
#define GEMValidation_TrackMatcher_h

/**\class TrackMatcher

 Description: Matching of tracks to SimTrack

 Original Author:  "Sven Dildick"
*/

#include "GEMCode/GEMValidation/src/SimHitMatcher.h"
#include "GEMCode/GEMValidation/src/GEMDigiMatcher.h"
#include "GEMCode/GEMValidation/src/CSCDigiMatcher.h"
#include "GEMCode/GEMValidation/src/RPCDigiMatcher.h"
#include "GEMCode/GEMValidation/src/CSCStubMatcher.h"

#include "GEMCode/GEMValidation/src/TFTrack.h" 
#include "GEMCode/GEMValidation/src/TFCand.h" 
#include "GEMCode/GEMValidation/src/GMTRegCand.h" 
#include "GEMCode/GEMValidation/src/GMTCand.h" 
#include "GEMCode/GEMValidation/src/L1Extra.h" 

#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFSectorProcessor.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h"
#include "L1Trigger/CSCTrackFinder/src/CSCTFDTReceiver.h"

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"


class TrackMatcher : public CSCStubMatcher
{
 public:
  /// constructor
  TrackMatcher(SimHitMatcher& sh, CSCDigiMatcher& dg, 
               GEMDigiMatcher& gem_dg, RPCDigiMatcher& rpc_dg, 
               CSCStubMatcher& st);
  /// destructor
  ~TrackMatcher();

  const std::vector<TFTrack*>& tfTracks() const {return tfTracks_;}
  const std::vector<TFCand*>& tfCands() const {return tfCands_;}
  const std::vector<GMTRegCand*>& gmtRegCands() const {return gmtRegCands_;}
  const std::vector<GMTCand*>& gmtCands() const {return gmtCands_;}
  const std::vector<L1Extra*>& l1Extras() const {return l1Extras_;}

  TFTrack* bestTFTrack(bool sortPtFirst=1);
  TFCand* bestTFCand(bool sortPtFirst=1);
  GMTRegCand* bestGMTRegCand(bool sortPtFirst=1);
  GMTCand* bestGMTCand(bool sortPtFirst=1);
  L1Extra* bestL1Extra(bool sortPtFirst=1);
  
 private:

  void init();
  void clear();

  void matchTfTrackToSimTrack(const L1CSCTrackCollection& tracks);
  void matchTfCandToSimTrack(const std::vector<L1MuRegionalCand>& tracks);
  void matchGmtRegCandToSimTrack(const L1MuRegionalCand& tracks);
  void matchGmtCandToSimTrack(const L1MuGMTExtendedCand& tracks);

  csctf::TrackStub buildTrackStub(const CSCCorrelatedLCTDigi& d, CSCDetId id);
  std::pair<float, float> intersectionEtaPhi(CSCDetId id, int wg, int hs);

  const SimHitMatcher* sh_matcher_;
  const GEMDigiMatcher* gem_digi_matcher_;
  const CSCDigiMatcher* csc_digi_matcher_;
  const RPCDigiMatcher* rpc_digi_matcher_;
  const CSCStubMatcher* csc_stub_matcher_;

  edm::InputTag cscTfTrackInputLabel_;
  edm::InputTag cscTfCandInputLabel_;
  edm::InputTag gmtRegCandInputLabel_;
  edm::InputTag gmtCandInputLabel_;
  edm::InputTag l1ExtraInputLabel_;

  int minBXTFTrack_, maxBXTFTrack_;
  int minBXTFCand_, maxBXTFCand_;
  int minBXGMTRegCand_, maxBXGMTRegCand_;
  int minBXGMTCand_, maxBXGMTCand_;
  int minBXL1Extra_, maxBXL1Extra_;

  int verboseTFTrack_;
  int verboseTFCand_;
  int verboseGMTRegCand_;
  int verboseGMTCand_;
  int verboseL1Extra_;

  std::vector<TFTrack*> tfTracks_;
  std::vector<TFCand*> tfCands_;
  std::vector<GMTRegCand*> gmtRegCands_;
  std::vector<GMTCand*> gmtCands_;
  std::vector<L1Extra*> l1Extras_;

  edm::ParameterSet ptLUTset_;
  edm::ParameterSet CSCTFSPset_;
  CSCTFPtLUT* ptLUT_;
  CSCTFSectorProcessor* my_SPs_[2][6];
  CSCSectorReceiverLUT* srLUTs_[5][6][2];
  CSCTFDTReceiver* my_dtrc_;
  unsigned long long  muScalesCacheID_;
  unsigned long long  muPtScaleCacheID_;

  edm::ESHandle<L1MuTriggerScales> muScales_;
  edm::ESHandle<L1MuTriggerPtScale> muPtScale_;

  bool hasMuScales_;
  bool hasMuPtScale_;
};

#endif
