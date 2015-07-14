#ifndef GEMCode_GEMValidation_L1GlobalMuonTriggerMatcher_h
#define GEMCode_GEMValidation_L1GlobalMuonTriggerMatcher_h

#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

typedef std::vector<L1MuRegionalCand> L1MuRegionalCandCollection;
typedef std::vector<L1MuGMTExtendedCand> L1MuGMTExtendedCandCollection;
typedef std::vector<L1MuGMTCand> L1MuGMTCandCollection;

class L1GlobalMuonTriggerMatcher : public BaseMatcher
{
 public:
  /// constructor
  L1GlobalMuonTriggerMatcher(SimHitMatcher& sh);
  /// destructor
  ~L1GlobalMuonTriggerMatcher();
  
  bool gmtCandInContainer(const L1MuGMTCand&, const L1MuGMTCandCollection&) const;
  bool isGmtCandMatched(const L1MuGMTCand&) const;

 private:
  
  void clear();
  void init(); 
  
  void matchRegionalCandCSCToSimTrack(const L1MuRegionalCandCollection&); 
  void matchRegionalCandDTToSimTrack(const L1MuRegionalCandCollection&); 
  void matchRegionalCandRPCbToSimTrack(const L1MuRegionalCandCollection&); 
  void matchRegionalCandRPCfToSimTrack(const L1MuRegionalCandCollection&); 
  void matchGMTCandToSimTrack(const L1MuGMTCandCollection&); 
  void matchL1ExtraMuonParticleToSimTrack(const l1extra::L1MuonParticleCollection&); 

  std::vector<edm::InputTag> gmtRegCandCSCInputLabel_;
  std::vector<edm::InputTag> gmtRegCandDTInputLabel_;
  std::vector<edm::InputTag> gmtRegCandRPCfInputLabel_;
  std::vector<edm::InputTag> gmtRegCandRPCbInputLabel_;
  std::vector<edm::InputTag> gmtCandInputLabel_;
  std::vector<edm::InputTag> l1ExtraMuonInputLabel_;

  int verboseGmtRegCandCSC_;
  int verboseGmtRegCandDT_;
  int verboseGmtRegCandRPCf_;
  int verboseGmtRegCandRPCb_;
  int verboseGmtCand_;
  int verboseL1ExtraMuon_;

  bool runGmtRegCandCSC_;
  bool runGmtRegCandDT_;
  bool runGmtRegCandRPCf_;
  bool runGmtRegCandRPCb_;
  bool runGmtCand_;
  bool runL1ExtraMuon_;

  int minBXGmtRegCandCSC_;
  int minBXGmtRegCandDT_;
  int minBXGmtRegCandRPCf_;
  int minBXGmtRegCandRPCb_;
  int minBXGmtCand_;
  int minBXL1ExtraMuon_;

  int maxBXGmtRegCandCSC_;
  int maxBXGmtRegCandDT_;
  int maxBXGmtRegCandRPCf_;
  int maxBXGmtRegCandRPCb_;
  int maxBXGmtCand_;
  int maxBXL1ExtraMuon_;

  float deltaRGmtRegCandCSC_;
  float deltaRGmtRegCandDT_;
  float deltaRGmtRegCandRPCf_;
  float deltaRGmtRegCandRPCb_;
  float deltaRGmtCand_;
  float deltaRL1ExtraMuon_;

  const SimHitMatcher* simhit_matcher_;

  L1MuRegionalCandCollection    matchedL1GmtCSCCands_;
  L1MuRegionalCandCollection    matchedL1GmtRPCfCands_;
  L1MuRegionalCandCollection    matchedL1GmtRPCbCands_;
  L1MuRegionalCandCollection    matchedL1GmtDTCands_;
  L1MuGMTCandCollection matchedL1GmtCands_;
  l1extra::L1MuonParticleCollection matchedL1MuonParticles_;
};

#endif
