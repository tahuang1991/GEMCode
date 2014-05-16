#ifndef SimMuL1_GEMCSCTriggerRateTree_h
#define SimMuL1_GEMCSCTriggerRateTree_h

// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

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

#include "GEMCode/SimMuL1/interface/MuGeometryHelpers.h"
#include "GEMCode/SimMuL1/interface/MatchCSCMuL1.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

struct MyALCT
{
  Int_t event, endcap, station, ring, chamber, bx;
};

struct MyCLCT
{
  Int_t event, endcap, station, ring, chamber, bx;
};

struct MyLCT
{
  Int_t event, endcap, station, ring, chamber, bx;
  Int_t quality, pattern, strip, wiregroup;
  Int_t hasGEM;
};

struct MyMPLCT
{
  Int_t event, endcap, station, ring, chamber, bx;
  Int_t quality, strip, wiregroup;
  Int_t hasGEM;
  Float_t etalut, philut;
};

struct MyTFTrack
{
  Int_t event, bx;
  Float_t pt, eta, phi;
  Int_t quality;
  Int_t hasME1a, hasME1b, hasME12, hasME13;
  Int_t hasME21, hasME22;
  Int_t hasME31, hasME32;
  Int_t hasME41, hasME42;
  Int_t hasGE11, hasGE21S, hasGE21L, hasME0;
  Int_t hasRE12, hasRE13;
  Int_t hasRE22, hasRE23;
  Int_t hasRE31, hasRE32, hasRE33;
  Int_t hasRE41, hasRE42, hasRE43;
};

struct MyTFCand
{
  Int_t event, bx;
  Float_t pt, eta, phi;
  Int_t quality;
  Int_t nStubs, nDetIds; 
  Int_t hasME1a, hasME1b, hasME12, hasME13;
  Int_t hasME21, hasME22;
  Int_t hasME31, hasME32;
  Int_t hasME41, hasME42;
  Int_t hasGE11, hasGE21S, hasGE21L, hasME0;
  Int_t hasRE12, hasRE13;
  Int_t hasRE22, hasRE23;
  Int_t hasRE31, hasRE32, hasRE33;
  Int_t hasRE41, hasRE42, hasRE43;
};

struct MyGMTRegCand
{
  Int_t event, bx;
  Float_t pt, eta, phi;
  Int_t quality;
  Int_t nStubs, nDetIds; 
  Int_t hasME1a, hasME1b, hasME12, hasME13;
  Int_t hasME21, hasME22;
  Int_t hasME31, hasME32;
  Int_t hasME41, hasME42;
  Int_t hasGE11, hasGE21S, hasGE21L, hasME0;
  Int_t hasRE12, hasRE13;
  Int_t hasRE22, hasRE23;
  Int_t hasRE31, hasRE32, hasRE33;
  Int_t hasRE41, hasRE42, hasRE43;
  Int_t isCSC,isDT,isRPCf,isRPCb;
};

struct MyGMT
{
  Int_t event, bx;
  Float_t pt, eta, phi;
  Int_t quality;
  Int_t nStubs, nDetIds; 
  Int_t isGoodSingleMuon;
  Int_t isGoodDoubleMuon;
  Int_t hasME1a, hasME1b, hasME12, hasME13;
  Int_t hasME21, hasME22;
  Int_t hasME31, hasME32;
  Int_t hasME41, hasME42;
  Int_t hasGE11, hasGE21S, hasGE21L, hasME0;
  Int_t hasRE12, hasRE13;
  Int_t hasRE22, hasRE23;
  Int_t hasRE31, hasRE32, hasRE33;
  Int_t hasRE41, hasRE42, hasRE43;
  Int_t isCSC,isDT,isRPCf,isRPCb;
};

class GEMCSCTriggerRateTree : public edm::EDAnalyzer 
{
 public:
  
  explicit GEMCSCTriggerRateTree(const edm::ParameterSet&);
  
  ~GEMCSCTriggerRateTree();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void beginJob();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void endJob();

  enum trig_cscs {MAX_STATIONS = 4, CSC_TYPES = 10};
  enum tfTrack{CSCTF=0, DTTF, RPCPAC};
  enum tfCand{CSCCand=0, DTCand, RPCbCand, RPCfCand};
  enum gmtRegCand{CSCGMT=0, DTGMT, RPCbGMT, RPCfGMT};
  
 private:
  
  // functions
  int getCSCType(CSCDetId& id);
  int isME11(int t);
  int getCSCSpecsType(CSCDetId& id);
  int cscTriggerSubsector(CSCDetId& id);

  // From Ingo:
  // calculates the weight of the event to reproduce a min bias
  //spectrum, from G. Wrochna's note CMSS 1997/096
  void setupTFModeHisto(TH1D* h);
  void runCSCTFSP(const CSCCorrelatedLCTDigiCollection*, const L1MuDTChambPhContainer*);

  std::pair<float, float> intersectionEtaPhi(CSCDetId id, int wg, int hs);
  csctf::TrackStub buildTrackStub(const CSCCorrelatedLCTDigi &d, CSCDetId id);

  void bookALCTTree();
  void bookCLCTTree();
  void bookLCTTree();
  void bookMPCLCTTree();
  void bookTFTrackTree();
  void bookTFCandTree();
  void bookGMTRegCandTree();
  void bookGMTCandTree();

  void intializeTree();

  void analyzeALCTRate(const edm::Event&);
  void analyzeCLCTRate(const edm::Event&);
  void analyzeLCTRate(const edm::Event&);
  void analyzeMPCLCTRate(const edm::Event&);

  void analyzeTFTrackRate(const edm::Event&);
  void analyzeTFTrackRate(const edm::Event&, enum tfTrack type);

  void analyzeTFCandRate(const edm::Event&);
  void analyzeTFCandRate(const edm::Event&, enum tfCand type);

  void analyzeGMTRegCandRate(const edm::Event&);
  void analyzeGMTRegCandRate(const edm::Event&, enum gmtRegCand type);

  void analyzeGMTCandRate(const edm::Event&);

  // parameters
  edm::ParameterSet simTrackMatching;
  edm::ParameterSet CSCTFSPset;
  edm::ParameterSet ptLUTset;

  CSCTFPtLUT* ptLUT;
  CSCTFSectorProcessor* my_SPs[2][6];
  CSCSectorReceiverLUT* srLUTs_[5][6][2];
  CSCTFDTReceiver* my_dtrc;
  unsigned long long  muScalesCacheID_;
  unsigned long long  muPtScaleCacheID_;

  edm::ESHandle< L1MuTriggerScales > muScales;
  edm::ESHandle< L1MuTriggerPtScale > muPtScale;

  // stubs
  edm::InputTag inputALCT_;
  int verboseALCT_;
  int minBXALCT_;
  int maxBXALCT_;

  edm::InputTag inputCLCT_;
  int verboseCLCT_;
  int minBXCLCT_;
  int maxBXCLCT_;

  edm::InputTag inputLCT_;
  int verboseLCT_;
  int minBXLCT_;
  int maxBXLCT_;

  edm::InputTag inputMPLCT_;
  int verboseMPLCT_;
  int minBXMPLCT_;
  int maxBXMPLCT_;
  
  edm::InputTag inputCSCTFTrack_;
  int verboseCSCTFTrack_;
  bool runCSCTFTrack_;
  int minBXCSCTFTrack_;
  int maxBXCSCTFTrack_;
  
  edm::InputTag inputCSCTFCand_;
  int verboseCSCTFCand_;
  bool runCSCTFCand_;
  int minBXCSCTFCand_;
  int maxBXCSCTFCand_;

  edm::InputTag inputRPCfTFCand_;
  int verboseRPCfTFCand_;
  int minBXRPCfTFCand_;
  int maxBXRPCfTFCand_;

  edm::InputTag inputRPCbTFCand_;
  int verboseRPCbTFCand_;
  int minBXRPCbTFCand_;
  int maxBXRPCbTFCand_;

  edm::InputTag inputGMTRegCand_;
  int verboseGMTRegCand_;
  bool runGMTRegCand_;
  int minBXGMTRegCand_;
  int maxBXGMTRegCand_;

  edm::InputTag inputGMTCand_;
  int verboseGMTCand_;
  bool runGMTCand_;
  int minBXGMTCand_;
  int maxBXGMTCand_;

  edm::InputTag inputL1Extra_;
  int verboseL1Extra_;
  int minBXL1Extra_;
  int maxBXL1Extra_;

  bool centralBxOnlyGMT_;
  bool doSelectEtaForGMTRates_;

  const CSCGeometry* cscGeometry;

  TTree* alct_tree_;
  TTree* clct_tree_;
  TTree* lct_tree_;
  TTree* mplct_tree_;
  TTree* tftrack_tree_;
  TTree* tfcand_tree_;
  TTree* gmtregcand_tree_;
  TTree* gmtcand_tree_;

  MyALCT alct_;
  MyCLCT clct_;
  MyLCT lct_;
  MyMPLCT mplct_;
  MyTFTrack tftrack_;
  MyTFCand tfcand_;
  MyGMTRegCand gmtregcand_;
  MyGMT gmtcand_;

  std::vector<MatchCSCMuL1::TFTRACK> rtTFTracks_;
  std::vector<MatchCSCMuL1::TFCAND> rtTFCands_;
  std::vector<MatchCSCMuL1::GMTREGCAND> rtGmtRegCscCands_;
  std::vector<MatchCSCMuL1::GMTREGCAND> rtGmtRegRpcfCands_;
  std::vector<MatchCSCMuL1::GMTREGCAND> rtGmtRegRpcbCands_;
  std::vector<MatchCSCMuL1::GMTREGCAND> rtGmtRegDtCands_;
  std::vector<MatchCSCMuL1::GMTCAND> rtGmtCands_;
  std::vector<MatchCSCMuL1::L1EXTRA> rtL1Extras_;

  edm::Handle<L1MuGMTReadoutCollection> hl1GmtCands_;
  edm::Handle<L1MuDTChambPhContainer> dttrig_;
  std::vector<L1MuGMTExtendedCand> l1GmtCands_;
  std::vector<L1MuGMTExtendedCand> l1GmtfCands_;
  std::vector<L1MuRegionalCand>    l1GmtCSCCands_;
  std::vector<L1MuRegionalCand>    l1GmtRPCfCands_;
  std::vector<L1MuRegionalCand>    l1GmtRPCbCands_;
  std::vector<L1MuRegionalCand>    l1GmtDTCands_;
  std::map<int, std::vector<L1MuRegionalCand> > l1GmtCSCCandsInBXs_;

  TH1D * h_events;
  int n_events;
};

#endif
