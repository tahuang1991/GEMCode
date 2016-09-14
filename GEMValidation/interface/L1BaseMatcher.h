#ifndef GEMCode_GEMValidation_L1BaseMatcher_h
#define GEMCode_GEMValidation_L1BaseMatcher_h

#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFSectorProcessor.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h"
#include "L1Trigger/CSCTrackFinder/src/CSCTFDTReceiver.h"

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

class L1BaseMatcher : public BaseMatcher
{
 public:
  /// constructor
  L1BaseMatcher(SimHitMatcher& sh, edm::ConsumesCollector & iC);
  /// destructor
  ~L1BaseMatcher();

 protected:
  
  edm::ESHandle<L1MuTriggerScales> muScalesHd_; 
  edm::ESHandle<L1MuTriggerPtScale> muPtScaleHd_; 
  
 private:
  
  void clear();
  void init(); 
  
  edm::ParameterSet ptLUTset_; 
  edm::ParameterSet CSCTFSPset_; 
  edm::ParameterSet srLUTset_; 

  CSCTFPtLUT* ptLUT_; 
  CSCTFSectorProcessor* my_SPs_[2][6]; 
  CSCSectorReceiverLUT* srLUTs_[5][6][2]; 
  CSCTFDTReceiver* dtrc_; 
  
  unsigned long long  muScalesCacheID_; 
  unsigned long long  muPtScaleCacheID_; 
  
  bool hasMuScales_; 
  bool hasMuPtScale_;   
};

#endif
