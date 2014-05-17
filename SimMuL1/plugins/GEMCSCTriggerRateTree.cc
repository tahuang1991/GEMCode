#include "GEMCode/SimMuL1/plugins/GEMCSCTriggerRateTree.h"

using namespace std;

GEMCSCTriggerRateTree::GEMCSCTriggerRateTree(const edm::ParameterSet& iConfig):
  simTrackMatching(iConfig.getParameter<edm::ParameterSet>("simTrackMatching")),
  CSCTFSPset(iConfig.getParameter<edm::ParameterSet>("sectorProcessor")),
  ptLUTset(CSCTFSPset.getParameter<edm::ParameterSet>("PTLUT")),
  ptLUT(0),
  centralBxOnlyGMT_(iConfig.getUntrackedParameter< bool >("centralBxOnlyGMT",false)),
  doSelectEtaForGMTRates_(iConfig.getUntrackedParameter< bool >("doSelectEtaForGMTRates",false))
{
  // stubs
  auto cscALCT = simTrackMatching.getParameter<edm::ParameterSet>("cscALCT");
  inputALCT_ = cscALCT.getParameter<edm::InputTag>("input");
  verboseALCT_ = cscALCT.getParameter<int>("verbose");
  minBXALCT_ = cscALCT.getParameter<int>("minBX");
  maxBXALCT_ = cscALCT.getParameter<int>("maxBX");

  auto cscCLCT = simTrackMatching.getParameter<edm::ParameterSet>("cscCLCT");
  inputCLCT_ = cscCLCT.getParameter<edm::InputTag>("input");
  verboseCLCT_ = cscCLCT.getParameter<int>("verbose");
  minBXCLCT_ = cscCLCT.getParameter<int>("minBX");
  maxBXCLCT_ = cscCLCT.getParameter<int>("maxBX");

  auto cscLCT = simTrackMatching.getParameter<edm::ParameterSet>("cscLCT");
  inputLCT_ = cscLCT.getParameter<edm::InputTag>("input");
  verboseLCT_ = cscLCT.getParameter<int>("verbose");
  minBXLCT_ = cscLCT.getParameter<int>("minBX");
  maxBXLCT_ = cscLCT.getParameter<int>("maxBX");

  auto cscMPLCT = simTrackMatching.getParameter<edm::ParameterSet>("cscMPLCT");
  inputMPLCT_ = cscMPLCT.getParameter<edm::InputTag>("input");
  verboseMPLCT_ = cscMPLCT.getParameter<int>("verbose");
  minBXMPLCT_ = cscMPLCT.getParameter<int>("minBX");
  maxBXMPLCT_ = cscMPLCT.getParameter<int>("maxBX");
  
  // TFTrack
  auto cscTfTrack = simTrackMatching.getParameter<edm::ParameterSet>("cscTfTrack");
  inputCSCTFTrack_ = cscTfTrack.getParameter<edm::InputTag>("input");
  runCSCTFTrack_ = cscTfTrack.getParameter<bool>("run");
  verboseCSCTFTrack_ = cscTfTrack.getParameter<int>("verbose");
  minBXCSCTFTrack_ = cscTfTrack.getParameter<int>("minBX");
  maxBXCSCTFTrack_ = cscTfTrack.getParameter<int>("maxBX");
  
  // TFCand
  auto cscTfCand = simTrackMatching.getParameter<edm::ParameterSet>("cscTfCand");
  inputCSCTFCand_ = cscTfCand.getParameter<edm::InputTag>("input");
  runCSCTFCand_ = cscTfCand.getParameter<bool>("run");
  verboseCSCTFCand_ = cscTfCand.getParameter<int>("verbose");
  minBXCSCTFCand_ = cscTfCand.getParameter<int>("minBX");
  maxBXCSCTFCand_ = cscTfCand.getParameter<int>("maxBX");

  auto rpcfTfCand = simTrackMatching.getParameter<edm::ParameterSet>("rpcfTfCand");
  inputRPCfTFCand_ = rpcfTfCand.getParameter<edm::InputTag>("input");
  verboseRPCfTFCand_ = rpcfTfCand.getParameter<int>("verbose");
  minBXRPCfTFCand_ = rpcfTfCand.getParameter<int>("minBX");
  maxBXRPCfTFCand_ = rpcfTfCand.getParameter<int>("maxBX");

  auto rpcbTfCand = simTrackMatching.getParameter<edm::ParameterSet>("rpcbTfCand");
  inputRPCbTFCand_ = rpcbTfCand.getParameter<edm::InputTag>("input");
  verboseRPCbTFCand_ = rpcbTfCand.getParameter<int>("verbose");
  minBXRPCbTFCand_ = rpcbTfCand.getParameter<int>("minBX");
  maxBXRPCbTFCand_ = rpcbTfCand.getParameter<int>("maxBX");

  // GMT
  auto gmtRegCand = simTrackMatching.getParameter<edm::ParameterSet>("gmtRegCand");
  verboseGMTRegCand_ = gmtRegCand.getParameter<int>("verbose");
  runGMTRegCand_ = gmtRegCand.getParameter<bool>("run");
  minBXGMTRegCand_ = gmtRegCand.getParameter<int>("minBX");
  maxBXGMTRegCand_ = gmtRegCand.getParameter<int>("maxBX");

  auto gmtCand = simTrackMatching.getParameter<edm::ParameterSet>("gmtCand");
  inputGMTCand_ = gmtCand.getParameter<edm::InputTag>("input"); 
  verboseGMTCand_ = gmtCand.getParameter<int>("verbose");
  runGMTCand_ = gmtCand.getParameter<bool>("run");
  minBXGMTCand_ = gmtCand.getParameter<int>("minBX");
  maxBXGMTCand_ = gmtCand.getParameter<int>("maxBX");

  auto l1Extra = simTrackMatching.getParameter<edm::ParameterSet>("l1Extra");
  inputL1Extra_ = l1Extra.getParameter<edm::InputTag>("input"); 
  verboseL1Extra_ = l1Extra.getParameter<int>("verbose");
  minBXL1Extra_ = l1Extra.getParameter<int>("minBX");
  maxBXL1Extra_ = l1Extra.getParameter<int>("maxBX");

  edm::ParameterSet srLUTset = CSCTFSPset.getParameter<edm::ParameterSet>("SRLUT");

  for(int e=0; e<2; e++) 
    for (int s=0; s<6; s++) 
      my_SPs[e][s] = nullptr;
  
  bool TMB07 = true;
  for(int endcap = 1; endcap<=2; endcap++)
  {
    for(int sector=1; sector<=6; sector++)
    {
      for(int station=1,fpga=0; station<=4 && fpga<5; station++)
      {
	if(station==1) for(int subSector=0; subSector<2; subSector++)
	  srLUTs_[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, subSector+1, station, srLUTset, TMB07);
	else
	  srLUTs_[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, 0, station, srLUTset, TMB07);
      }
    }
  }

  my_dtrc = new CSCTFDTReceiver();

  // cache flags for event setup records
  muScalesCacheID_ = 0ULL ;
  muPtScaleCacheID_ = 0ULL ;

  bookALCTTree();
  bookCLCTTree();
  bookLCTTree();
  bookMPCLCTTree();
  if (runCSCTFTrack_) bookTFTrackTree();
  if (runCSCTFTrack_ and runCSCTFCand_) bookTFCandTree();
  if (runCSCTFTrack_ and runCSCTFCand_ and runGMTRegCand_) bookGMTRegCandTree();
  if (runCSCTFTrack_ and runCSCTFCand_ and runGMTRegCand_ and runGMTCand_) bookGMTCandTree();

  n_events = 0;
}

// ================================================================================================
GEMCSCTriggerRateTree::~GEMCSCTriggerRateTree()
{
  if(ptLUT) delete ptLUT;
  ptLUT = nullptr;

  for(int e=0; e<2; e++) for (int s=0; s<6; s++){
      if  (my_SPs[e][s]) delete my_SPs[e][s];
      my_SPs[e][s] = nullptr;

      for(int fpga=0; fpga<5; fpga++)
	{
	  if (srLUTs_[fpga][s][e]) delete srLUTs_[fpga][s][e];
	  srLUTs_[fpga][s][e] = nullptr;
	}
    }
  
  if(my_dtrc) delete my_dtrc;
  my_dtrc = nullptr;
}

// ================================================================================================
void
GEMCSCTriggerRateTree::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
  edm::ESHandle< CSCGeometry > cscGeom;
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
  cscGeometry = &*cscGeom;
  CSCTriggerGeometry::setGeometry(cscGeometry);

  if (iSetup.get< L1MuTriggerScalesRcd >().cacheIdentifier() != muScalesCacheID_ or
      iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier() != muPtScaleCacheID_ )
    {
      iSetup.get< L1MuTriggerScalesRcd >().get( muScales );
      iSetup.get< L1MuTriggerPtScaleRcd >().get( muPtScale );
      if (ptLUT) delete ptLUT;  
      ptLUT = new CSCTFPtLUT(ptLUTset, muScales.product(), muPtScale.product());
      
      for(int e=0; e<2; e++) for (int s=0; s<6; s++){
	if  (my_SPs[e][s]) delete my_SPs[e][s];
	my_SPs[e][s] = new CSCTFSectorProcessor(e+1, s+1, CSCTFSPset, true, muScales.product(), muPtScale.product());
	my_SPs[e][s]->initialize(iSetup);
      }
      muScalesCacheID_  = iSetup.get< L1MuTriggerScalesRcd >().cacheIdentifier();
      muPtScaleCacheID_ = iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier();
    }
}

// ================================================================================================
void 
GEMCSCTriggerRateTree::beginJob()
{
  edm::Service<TFileService> fs;

  h_events = fs->make<TH1D>("h_events","h_events",1,0,1);
}

void 
GEMCSCTriggerRateTree::endJob()
{
  h_events->SetBinContent(1,n_events);
}

void 
GEMCSCTriggerRateTree::intializeTree()
{
  alct_.event = -999;
  alct_.endcap = -999;  
  alct_.station = -999;   
  alct_.ring = -999;   
  alct_.chamber = -999;   
  alct_.bx = -999;  

  clct_.event = -999;
  clct_.endcap = -999;  
  clct_.station = -999;   
  clct_.ring = -999;   
  clct_.chamber = -999;   
  clct_.bx = -999;  

  lct_.event = -999;
  lct_.endcap = -999;  
  lct_.station = -999;   
  lct_.ring = -999;   
  lct_.chamber = -999;   
  lct_.bx = -999;  
  lct_.quality = -999;
  lct_.strip = -999;
  lct_.wiregroup = -999;
  lct_.hasGEM = 0;

  mplct_.event = -999;
  mplct_.endcap = -999;  
  mplct_.station = -999;   
  mplct_.ring = -999;   
  mplct_.chamber = -999;   
  mplct_.bx = -999;  
  mplct_.quality = -999;
  mplct_.strip = -999;
  mplct_.wiregroup = -999;
  mplct_.hasGEM = 0;
  mplct_.etalut = -999;
  mplct_.philut = -999;

  tftrack_.event = -999;
  tftrack_.bx = -999;
  tftrack_.pt = -999;
  tftrack_.eta = -999;
  tftrack_.phi = -999;
  tftrack_.quality = -999 ;
  tftrack_.hasME1a = 0;
  tftrack_.hasME1b = 0; 
  tftrack_.hasME12 = 0; 
  tftrack_.hasME13 = 0;
  tftrack_.hasME21 = 0; 
  tftrack_.hasME22 = 0;
  tftrack_.hasME31 = 0; 
  tftrack_.hasME32 = 0;
  tftrack_.hasME41 = 0; 
  tftrack_.hasME42 = 0;
  tftrack_.hasGE11 = 0; 
  tftrack_.hasGE21S = 0; 
  tftrack_.hasGE21L = 0; 
  tftrack_.hasME0 = 0;
  tftrack_.hasRE12 = 0; 
  tftrack_.hasRE13 = 0;
  tftrack_.hasRE22 = 0;
  tftrack_.hasRE23 = 0; 
  tftrack_.hasRE31 = 0; 
  tftrack_.hasRE32 = 0;
  tftrack_.hasRE33 = 0;
  tftrack_.hasRE41 = 0; 
  tftrack_.hasRE42 = 0;
  tftrack_.hasRE43 = 0; 

  tfcand_.event = -999;
  tfcand_.bx = -999;
  tfcand_.pt = -999;
  tfcand_.eta = -999;
  tfcand_.phi = -999;
  tfcand_.nStubs = 0;
  tfcand_.nDetIds = 0;
  tfcand_.quality = -999;
  tfcand_.hasME1a = 0;
  tfcand_.hasME1b = 0; 
  tfcand_.hasME12 = 0; 
  tfcand_.hasME13 = 0;
  tfcand_.hasME21 = 0; 
  tfcand_.hasME22 = 0;
  tfcand_.hasME31 = 0; 
  tfcand_.hasME32 = 0;
  tfcand_.hasME41 = 0; 
  tfcand_.hasME42 = 0;
  tfcand_.hasGE11 = 0; 
  tfcand_.hasGE21S = 0; 
  tfcand_.hasGE21L = 0; 
  tfcand_.hasME0 = 0;
  tfcand_.hasRE12 = 0; 
  tfcand_.hasRE13 = 0;
  tfcand_.hasRE22 = 0;
  tfcand_.hasRE23 = 0; 
  tfcand_.hasRE31 = 0; 
  tfcand_.hasRE32 = 0;
  tfcand_.hasRE33 = 0;
  tfcand_.hasRE41 = 0; 
  tfcand_.hasRE42 = 0;
  tfcand_.hasRE43 = 0; 

  gmtregcand_.event = -999;
  gmtregcand_.bx = -999;
  gmtregcand_.pt = -999;
  gmtregcand_.eta = -999;
  gmtregcand_.phi = -999;
  gmtregcand_.quality = -999;
  gmtregcand_.nStubs = 0;
  gmtregcand_.nDetIds = 0;
  gmtregcand_.hasME1a = 0;
  gmtregcand_.hasME1b = 0; 
  gmtregcand_.hasME12 = 0; 
  gmtregcand_.hasME13 = 0;
  gmtregcand_.hasME21 = 0; 
  gmtregcand_.hasME22 = 0;
  gmtregcand_.hasME31 = 0; 
  gmtregcand_.hasME32 = 0;
  gmtregcand_.hasME41 = 0; 
  gmtregcand_.hasME42 = 0;
  gmtregcand_.hasGE11 = 0; 
  gmtregcand_.hasGE21S = 0; 
  gmtregcand_.hasGE21L = 0; 
  gmtregcand_.hasME0 = 0;
  gmtregcand_.hasRE12 = 0; 
  gmtregcand_.hasRE13 = 0;
  gmtregcand_.hasRE22 = 0;
  gmtregcand_.hasRE23 = 0; 
  gmtregcand_.hasRE31 = 0; 
  gmtregcand_.hasRE32 = 0;
  gmtregcand_.hasRE33 = 0;
  gmtregcand_.hasRE41 = 0; 
  gmtregcand_.hasRE42 = 0;
  gmtregcand_.hasRE43 = 0; 

  gmtcand_.event = -999;
  gmtcand_.bx = -999;
  gmtcand_.pt = -999;
  gmtcand_.eta = -999;
  gmtcand_.phi = -999;
  gmtcand_.quality = -999;
  gmtcand_.nStubs = 0;
  gmtcand_.nDetIds = 0;
  gmtcand_.hasME1a = 0;
  gmtcand_.hasME1b = 0; 
  gmtcand_.hasME12 = 0; 
  gmtcand_.hasME13 = 0;
  gmtcand_.hasME21 = 0; 
  gmtcand_.hasME22 = 0;
  gmtcand_.hasME31 = 0; 
  gmtcand_.hasME32 = 0;
  gmtcand_.hasME41 = 0; 
  gmtcand_.hasME42 = 0;
  gmtcand_.hasGE11 = 0; 
  gmtcand_.hasGE21S = 0; 
  gmtcand_.hasGE21L = 0; 
  gmtcand_.hasME0 = 0;
  gmtcand_.hasRE12 = 0; 
  gmtcand_.hasRE13 = 0;
  gmtcand_.hasRE22 = 0;
  gmtcand_.hasRE23 = 0; 
  gmtcand_.hasRE31 = 0; 
  gmtcand_.hasRE32 = 0;
  gmtcand_.hasRE33 = 0;
  gmtcand_.hasRE41 = 0; 
  gmtcand_.hasRE42 = 0;
  gmtcand_.hasRE43 = 0; 
}

// ================================================================================================
void 
GEMCSCTriggerRateTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++n_events;

  // need to reset here
  intializeTree();

  // clear the vectors
  rtTFTracks_.clear();
  rtTFCands_.clear();
  rtGmtRegCscCands_.clear();
  rtGmtRegRpcfCands_.clear();
  rtGmtRegRpcbCands_.clear();
  rtGmtRegDtCands_.clear();
  rtGmtCands_.clear();
  rtL1Extras_.clear();

  l1GmtCands_.clear();
  l1GmtfCands_.clear();
  l1GmtCSCCands_.clear();
  l1GmtRPCfCands_.clear();
  l1GmtRPCbCands_.clear();
  l1GmtDTCands_.clear();
  l1GmtCSCCandsInBXs_.clear();

//   // DT primitives for input to TF for debugging
//   iEvent.getByLabel("simDtTriggerPrimitiveDigis", dttrig_);
//   const L1MuDTChambPhContainer* dttrigs = dttrig_.product();

  // GMT readout collection
  iEvent.getByLabel(inputGMTCand_, hl1GmtCands_);

  if (centralBxOnlyGMT_){
    // Get GMT candidates from central bunch crossing only
    l1GmtCands_ = hl1GmtCands_->getRecord().getGMTCands() ;
    l1GmtfCands_ = hl1GmtCands_->getRecord().getGMTFwdCands() ;
    l1GmtCSCCands_ = hl1GmtCands_->getRecord().getCSCCands() ;
    l1GmtRPCfCands_ = hl1GmtCands_->getRecord().getFwdRPCCands() ;
    l1GmtRPCbCands_ = hl1GmtCands_->getRecord().getBrlRPCCands() ;
    l1GmtDTCands_ = hl1GmtCands_->getRecord().getDTBXCands() ;
    l1GmtCSCCandsInBXs_[hl1GmtCands_->getRecord().getBxInEvent()] = l1GmtCSCCands_;
  }
  else {
    // Get GMT candidates from all bunch crossings
    std::vector<L1MuGMTReadoutRecord> gmt_records = hl1GmtCands_->getRecords();
    for ( std::vector< L1MuGMTReadoutRecord >::const_iterator rItr=gmt_records.begin(); rItr!=gmt_records.end() ; ++rItr ){
      if (rItr->getBxInEvent() < minBXGMTCand_ or rItr->getBxInEvent() > maxBXGMTCand_) continue;
      
      std::vector<L1MuGMTExtendedCand> GMTCands = rItr->getGMTCands();
      for ( std::vector<L1MuGMTExtendedCand>::const_iterator  cItr = GMTCands.begin() ; cItr != GMTCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtCands_.push_back(*cItr);
      
      std::vector<L1MuGMTExtendedCand> GMTfCands = rItr->getGMTFwdCands();
      for ( std::vector<L1MuGMTExtendedCand>::const_iterator  cItr = GMTfCands.begin() ; cItr != GMTfCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtfCands_.push_back(*cItr);
      
      std::vector<L1MuRegionalCand> CSCCands = rItr->getCSCCands();
      l1GmtCSCCandsInBXs_[rItr->getBxInEvent()] = CSCCands;
      for ( std::vector<L1MuRegionalCand>::const_iterator  cItr = CSCCands.begin() ; cItr != CSCCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtCSCCands_.push_back(*cItr);
      
      std::vector<L1MuRegionalCand> RPCfCands = rItr->getFwdRPCCands();
      for ( std::vector<L1MuRegionalCand>::const_iterator  cItr = RPCfCands.begin() ; cItr != RPCfCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtRPCfCands_.push_back(*cItr);
      
      std::vector<L1MuRegionalCand> RPCbCands = rItr->getBrlRPCCands();
      for ( std::vector<L1MuRegionalCand>::const_iterator  cItr = RPCbCands.begin() ; cItr != RPCbCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtRPCbCands_.push_back(*cItr);
      
      std::vector<L1MuRegionalCand> DTCands = rItr->getDTBXCands();
      for ( std::vector<L1MuRegionalCand>::const_iterator  cItr = DTCands.begin() ; cItr != DTCands.end() ; ++cItr )
	if (!cItr->empty()) l1GmtDTCands_.push_back(*cItr);
    }
  }

  analyzeALCTRate(iEvent);
  analyzeCLCTRate(iEvent);
  analyzeLCTRate(iEvent);
  analyzeMPCLCTRate(iEvent);
  if (runCSCTFTrack_) analyzeTFTrackRate(iEvent);
  if (runCSCTFTrack_ and runCSCTFCand_) analyzeTFCandRate(iEvent);
  if (runCSCTFTrack_ and runCSCTFCand_ and runGMTRegCand_) analyzeGMTRegCandRate(iEvent);
  if (runCSCTFTrack_ and runCSCTFCand_ and runGMTRegCand_ and runGMTCand_) analyzeGMTCandRate(iEvent);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookALCTTree()
{
  edm::Service< TFileService > fs;
  alct_tree_ = fs->make<TTree>("ALCT", "ALCT");
  alct_tree_->Branch("event",&alct_.event);
  alct_tree_->Branch("bx",&alct_.bx);
  alct_tree_->Branch("endcap",&alct_.endcap);
  alct_tree_->Branch("station",&alct_.station);
  alct_tree_->Branch("ring",&alct_.ring);
  alct_tree_->Branch("chamber",&alct_.chamber);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookCLCTTree()
{
  edm::Service< TFileService > fs;
  clct_tree_ = fs->make<TTree>("CLCT", "CLCT");
  clct_tree_->Branch("event",&clct_.event);
  clct_tree_->Branch("bx",&clct_.bx);
  clct_tree_->Branch("endcap",&clct_.endcap);
  clct_tree_->Branch("station",&clct_.station);
  clct_tree_->Branch("ring",&clct_.ring);
  clct_tree_->Branch("chamber",&clct_.chamber);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookLCTTree()
{
  edm::Service< TFileService > fs;
  lct_tree_ = fs->make<TTree>("LCT", "LCT");
  lct_tree_->Branch("event",&lct_.event);
  lct_tree_->Branch("bx",&lct_.bx);
  lct_tree_->Branch("endcap",&lct_.endcap);
  lct_tree_->Branch("station",&lct_.station);
  lct_tree_->Branch("ring",&lct_.ring);
  lct_tree_->Branch("chamber",&lct_.chamber);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookMPCLCTTree()
{
  edm::Service< TFileService > fs;
  mplct_tree_ = fs->make<TTree>("MPLCT", "MPLCT");
  mplct_tree_->Branch("event",&mplct_.event);
  mplct_tree_->Branch("bx",&mplct_.bx);
  mplct_tree_->Branch("endcap",&mplct_.endcap);
  mplct_tree_->Branch("station",&mplct_.station);
  mplct_tree_->Branch("ring",&mplct_.ring);
  mplct_tree_->Branch("chamber",&mplct_.chamber);
  mplct_tree_->Branch("etalut",&mplct_.etalut);
  mplct_tree_->Branch("philut",&mplct_.philut);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookTFTrackTree()
{
  edm::Service< TFileService > fs;
  tftrack_tree_ = fs->make<TTree>("TFTrack", "TFTrack");
  tftrack_tree_->Branch("event",&tftrack_.event);
  tftrack_tree_->Branch("bx",&tftrack_.bx);
  tftrack_tree_->Branch("pt",&tftrack_.pt);
  tftrack_tree_->Branch("eta",&tftrack_.eta);
  tftrack_tree_->Branch("phi",&tftrack_.phi);
  tftrack_tree_->Branch("hasME1a",&tftrack_.hasME1a);
  tftrack_tree_->Branch("hasME1b",&tftrack_.hasME1b);
  tftrack_tree_->Branch("hasME12",&tftrack_.hasME12);
  tftrack_tree_->Branch("hasME13",&tftrack_.hasME13);
  tftrack_tree_->Branch("hasME21",&tftrack_.hasME21);
  tftrack_tree_->Branch("hasME22",&tftrack_.hasME22);
  tftrack_tree_->Branch("hasME31",&tftrack_.hasME31);
  tftrack_tree_->Branch("hasME32",&tftrack_.hasME32);
  tftrack_tree_->Branch("hasME41",&tftrack_.hasME41);
  tftrack_tree_->Branch("hasME42",&tftrack_.hasME42);
  tftrack_tree_->Branch("hasRE12",&tftrack_.hasRE12);
  tftrack_tree_->Branch("hasRE13",&tftrack_.hasRE13);
  tftrack_tree_->Branch("hasRE22",&tftrack_.hasRE22);
  tftrack_tree_->Branch("hasRE23",&tftrack_.hasRE23);
  tftrack_tree_->Branch("hasRE31",&tftrack_.hasRE31);
  tftrack_tree_->Branch("hasRE32",&tftrack_.hasRE32);
  tftrack_tree_->Branch("hasRE33",&tftrack_.hasRE33);
  tftrack_tree_->Branch("hasRE41",&tftrack_.hasRE41);
  tftrack_tree_->Branch("hasRE42",&tftrack_.hasRE42);
  tftrack_tree_->Branch("hasRE43",&tftrack_.hasRE43);
  tftrack_tree_->Branch("hasGE11",&tftrack_.hasGE11);
  tftrack_tree_->Branch("hasGE21S",&tftrack_.hasGE21S);
  tftrack_tree_->Branch("hasGE21L",&tftrack_.hasGE21L);
  tftrack_tree_->Branch("hasME0",&tftrack_.hasME0);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookTFCandTree()
{
  edm::Service< TFileService > fs;
  tfcand_tree_ = fs->make<TTree>("TFCand", "TFCand");
  tfcand_tree_->Branch("event",&tfcand_.event);
  tfcand_tree_->Branch("bx",&tfcand_.bx);
  tfcand_tree_->Branch("pt",&tfcand_.pt);
  tfcand_tree_->Branch("eta",&tfcand_.eta);
  tfcand_tree_->Branch("phi",&tfcand_.phi);
  tfcand_tree_->Branch("nDetIds",&tfcand_.nDetIds);
  tfcand_tree_->Branch("nStubs",&tfcand_.nStubs);
  tfcand_tree_->Branch("hasME1a",&tfcand_.hasME1a);
  tfcand_tree_->Branch("hasME1b",&tfcand_.hasME1b);
  tfcand_tree_->Branch("hasME12",&tfcand_.hasME12);
  tfcand_tree_->Branch("hasME13",&tfcand_.hasME13);
  tfcand_tree_->Branch("hasME21",&tfcand_.hasME21);
  tfcand_tree_->Branch("hasME22",&tfcand_.hasME22);
  tfcand_tree_->Branch("hasME31",&tfcand_.hasME31);
  tfcand_tree_->Branch("hasME32",&tfcand_.hasME32);
  tfcand_tree_->Branch("hasME41",&tfcand_.hasME41);
  tfcand_tree_->Branch("hasME42",&tfcand_.hasME42);
  tfcand_tree_->Branch("hasRE12",&tfcand_.hasRE12);
  tfcand_tree_->Branch("hasRE13",&tfcand_.hasRE13);
  tfcand_tree_->Branch("hasRE22",&tfcand_.hasRE22);
  tfcand_tree_->Branch("hasRE23",&tfcand_.hasRE23);
  tfcand_tree_->Branch("hasRE31",&tfcand_.hasRE31);
  tfcand_tree_->Branch("hasRE32",&tfcand_.hasRE32);
  tfcand_tree_->Branch("hasRE33",&tfcand_.hasRE33);
  tfcand_tree_->Branch("hasRE41",&tfcand_.hasRE41);
  tfcand_tree_->Branch("hasRE42",&tfcand_.hasRE42);
  tfcand_tree_->Branch("hasRE43",&tfcand_.hasRE43);
  tfcand_tree_->Branch("hasGE11",&tfcand_.hasGE11);
  tfcand_tree_->Branch("hasGE21S",&tfcand_.hasGE21S);
  tfcand_tree_->Branch("hasGE21L",&tfcand_.hasGE21L);
  tfcand_tree_->Branch("hasME0",&tfcand_.hasME0);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookGMTRegCandTree()
{
  edm::Service< TFileService > fs;
  gmtregcand_tree_ = fs->make<TTree>("GMTRegCand", "GMTRegCand");
  gmtregcand_tree_->Branch("event",&gmtregcand_.event);
  gmtregcand_tree_->Branch("bx",&gmtregcand_.bx);
  gmtregcand_tree_->Branch("pt",&gmtregcand_.pt);
  gmtregcand_tree_->Branch("eta",&gmtregcand_.eta);
  gmtregcand_tree_->Branch("phi",&gmtregcand_.phi);
  gmtregcand_tree_->Branch("quality",&gmtregcand_.quality);
  gmtregcand_tree_->Branch("nStubs",&gmtregcand_.nStubs);
  gmtregcand_tree_->Branch("nDetIds",&gmtregcand_.nDetIds);
  gmtregcand_tree_->Branch("hasME1a",&gmtregcand_.hasME1a);
  gmtregcand_tree_->Branch("hasME1b",&gmtregcand_.hasME1b);
  gmtregcand_tree_->Branch("hasME12",&gmtregcand_.hasME12);
  gmtregcand_tree_->Branch("hasME13",&gmtregcand_.hasME13);
  gmtregcand_tree_->Branch("hasME21",&gmtregcand_.hasME21);
  gmtregcand_tree_->Branch("hasME22",&gmtregcand_.hasME22);
  gmtregcand_tree_->Branch("hasME31",&gmtregcand_.hasME31);
  gmtregcand_tree_->Branch("hasME32",&gmtregcand_.hasME32);
  gmtregcand_tree_->Branch("hasME41",&gmtregcand_.hasME41);
  gmtregcand_tree_->Branch("hasME42",&gmtregcand_.hasME42);
  gmtregcand_tree_->Branch("hasRE12",&gmtregcand_.hasRE12);
  gmtregcand_tree_->Branch("hasRE13",&gmtregcand_.hasRE13);
  gmtregcand_tree_->Branch("hasRE22",&gmtregcand_.hasRE22);
  gmtregcand_tree_->Branch("hasRE23",&gmtregcand_.hasRE23);
  gmtregcand_tree_->Branch("hasRE31",&gmtregcand_.hasRE31);
  gmtregcand_tree_->Branch("hasRE32",&gmtregcand_.hasRE32);
  gmtregcand_tree_->Branch("hasRE33",&gmtregcand_.hasRE33);
  gmtregcand_tree_->Branch("hasRE41",&gmtregcand_.hasRE41);
  gmtregcand_tree_->Branch("hasRE42",&gmtregcand_.hasRE42);
  gmtregcand_tree_->Branch("hasRE43",&gmtregcand_.hasRE43);
  gmtregcand_tree_->Branch("hasGE11",&gmtregcand_.hasGE11);
  gmtregcand_tree_->Branch("hasGE21S",&gmtregcand_.hasGE21S);
  gmtregcand_tree_->Branch("hasGE21L",&gmtregcand_.hasGE21L);
  gmtregcand_tree_->Branch("hasME0",&gmtregcand_.hasME0);
  gmtregcand_tree_->Branch("isCSC",&gmtregcand_.isCSC);
  gmtregcand_tree_->Branch("isDT",&gmtregcand_.isDT);
  gmtregcand_tree_->Branch("isRPCb",&gmtregcand_.isRPCb);
  gmtregcand_tree_->Branch("isRPCf",&gmtregcand_.isRPCf);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::bookGMTCandTree()
{
  edm::Service< TFileService > fs;
  gmtcand_tree_ = fs->make<TTree>("GMTCand", "GMTCand");
  gmtcand_tree_->Branch("event",&gmtcand_.event);
  gmtcand_tree_->Branch("bx",&gmtcand_.bx);
  gmtcand_tree_->Branch("pt",&gmtcand_.pt);
  gmtcand_tree_->Branch("eta",&gmtcand_.eta);
  gmtcand_tree_->Branch("phi",&gmtcand_.phi);
  gmtcand_tree_->Branch("quality",&gmtcand_.quality);
  gmtcand_tree_->Branch("nStubs",&gmtcand_.nStubs);
  gmtcand_tree_->Branch("nDetIds",&gmtcand_.nDetIds);
  gmtcand_tree_->Branch("hasME1a",&gmtcand_.hasME1a);
  gmtcand_tree_->Branch("hasME1b",&gmtcand_.hasME1b);
  gmtcand_tree_->Branch("hasME12",&gmtcand_.hasME12);
  gmtcand_tree_->Branch("hasME13",&gmtcand_.hasME13);
  gmtcand_tree_->Branch("hasME21",&gmtcand_.hasME21);
  gmtcand_tree_->Branch("hasME22",&gmtcand_.hasME22);
  gmtcand_tree_->Branch("hasME31",&gmtcand_.hasME31);
  gmtcand_tree_->Branch("hasME32",&gmtcand_.hasME32);
  gmtcand_tree_->Branch("hasME41",&gmtcand_.hasME41);
  gmtcand_tree_->Branch("hasME42",&gmtcand_.hasME42);
  gmtcand_tree_->Branch("hasRE12",&gmtcand_.hasRE12);
  gmtcand_tree_->Branch("hasRE13",&gmtcand_.hasRE13);
  gmtcand_tree_->Branch("hasRE22",&gmtcand_.hasRE22);
  gmtcand_tree_->Branch("hasRE23",&gmtcand_.hasRE23);
  gmtcand_tree_->Branch("hasRE31",&gmtcand_.hasRE31);
  gmtcand_tree_->Branch("hasRE32",&gmtcand_.hasRE32);
  gmtcand_tree_->Branch("hasRE33",&gmtcand_.hasRE33);
  gmtcand_tree_->Branch("hasRE41",&gmtcand_.hasRE41);
  gmtcand_tree_->Branch("hasRE42",&gmtcand_.hasRE42);
  gmtcand_tree_->Branch("hasRE43",&gmtcand_.hasRE43);
  gmtcand_tree_->Branch("hasGE11",&gmtcand_.hasGE11);
  gmtcand_tree_->Branch("hasGE21S",&gmtcand_.hasGE21S);
  gmtcand_tree_->Branch("hasGE21L",&gmtcand_.hasGE21L);
  gmtcand_tree_->Branch("hasME0",&gmtcand_.hasME0);
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeALCTRate(const edm::Event& iEvent)
{
  edm::Handle< CSCALCTDigiCollection > halcts;
  iEvent.getByLabel(inputALCT_,  halcts);
  const CSCALCTDigiCollection* alcts = halcts.product();
  
  for (CSCALCTDigiCollection::DigiRangeIterator  adetUnitIt = alcts->begin(); adetUnitIt != alcts->end(); ++adetUnitIt)
  {
    CSCDetId detId((*adetUnitIt).first);
    auto range = (*adetUnitIt).second;
    for (CSCALCTDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt)
    {
      if (!(*digiIt).isValid()) continue;
      const int bx((*digiIt).getBX());
      if (bx < minBXALCT_ or bx > maxBXALCT_) continue;
      alct_.event = iEvent.id().event();
      alct_.endcap = detId.zendcap();
      alct_.station = detId.station();
      alct_.ring = detId.ring();
      alct_.chamber = detId.chamber();
      alct_.bx = bx - 6;
      alct_tree_->Fill();

      // debug
      if (verboseALCT_){
        cout << "------------------------------------------------------------------------------" << endl;         
        cout << "Event " << alct_.event << ", detId " << detId << ", ALCT " << digiIt-range.first << endl;
        cout << "endcap " << alct_.endcap << ", station " << alct_.station << ", ring " << alct_.ring << ", chamber " << alct_.chamber << endl;
        cout << *digiIt << endl;
      }
    }
  }
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeCLCTRate(const edm::Event& iEvent)
{
  edm::Handle< CSCCLCTDigiCollection > hclcts;
  iEvent.getByLabel(inputCLCT_,  hclcts);
  const CSCCLCTDigiCollection* clcts = hclcts.product();

  for (CSCCLCTDigiCollection::DigiRangeIterator  adetUnitIt = clcts->begin(); adetUnitIt != clcts->end(); ++adetUnitIt)
  {
    CSCDetId detId((*adetUnitIt).first);
    auto range = (*adetUnitIt).second;
    for (CSCCLCTDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt)
    {
      if (!(*digiIt).isValid()) continue;
      const int bx((*digiIt).getBX());
      if (bx < minBXCLCT_ or bx > maxBXCLCT_) continue;
      clct_.event = iEvent.id().event();
      clct_.endcap = detId.zendcap();
      clct_.station = detId.station();
      clct_.ring = detId.ring();
      clct_.chamber = detId.chamber();
      clct_.bx = bx - 6;
      clct_tree_->Fill();

      // debug
      if (verboseCLCT_){
        cout << "------------------------------------------------------------------------------" << endl;         
        cout << "Event " << clct_.event << ", detId " << detId << ", CLCT " << digiIt-range.first << endl;
        cout << "endcap " << clct_.endcap << ", station " << clct_.station << ", ring " << clct_.ring << ", chamber " << clct_.chamber << endl;
        cout << *digiIt << endl;
      }
    }
  }
}


// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeLCTRate(const edm::Event& iEvent)
{
  edm::Handle< CSCCorrelatedLCTDigiCollection > lcts_tmb;
  iEvent.getByLabel(inputLCT_,  lcts_tmb);
  const CSCCorrelatedLCTDigiCollection* lcts = lcts_tmb.product();

  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator detUnitIt = lcts->begin(); detUnitIt != lcts->end(); detUnitIt++) 
  {
    const CSCDetId& detId = (*detUnitIt).first;
    auto range = (*detUnitIt).second;
    for (CSCCorrelatedLCTDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) 
    {
      if (!(*digiIt).isValid()) continue;
      const int bx((*digiIt).getBX());
      if (bx < minBXLCT_ or bx > maxBXLCT_) continue;
      lct_.event = iEvent.id().event();
      lct_.endcap = detId.zendcap();
      lct_.station = detId.station();
      lct_.ring = detId.ring();
      lct_.chamber = detId.chamber();
      lct_.bx = bx - 6;
      lct_tree_->Fill();

      // debug
      if (verboseLCT_){
        cout << "------------------------------------------------------------------------------" << endl;         
        cout << "Event " << lct_.event << ", detId " << detId << ", LCT " << digiIt-range.first << endl;
        cout << "endcap " << lct_.endcap << ", station " << lct_.station << ", ring " << lct_.ring << ", chamber " << lct_.chamber << endl;
        cout << *digiIt << endl;
      }
    }
  }
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeMPCLCTRate(const edm::Event& iEvent)
{
  edm::Handle< CSCCorrelatedLCTDigiCollection > lcts_mpc;
  iEvent.getByLabel(inputMPLCT_,  lcts_mpc);
  const CSCCorrelatedLCTDigiCollection* mplcts = lcts_mpc.product();

  for (auto detUnitIt = mplcts->begin(); detUnitIt != mplcts->end(); detUnitIt++) 
  {
    const CSCDetId& detId = (*detUnitIt).first;
    auto range = (*detUnitIt).second;
    for (auto digiIt = range.first; digiIt != range.second; digiIt++) 
    {
      if (!(*digiIt).isValid()) continue;
      const int bx((*digiIt).getBX());
      if (bx < minBXMPLCT_ or bx > maxBXMPLCT_) continue;
      mplct_.event = iEvent.id().event();
      mplct_.endcap = detId.zendcap();
      mplct_.station = detId.station();
      mplct_.ring = detId.ring();
      mplct_.chamber = detId.chamber();
      mplct_.bx = bx - 6;
      const csctf::TrackStub stub(buildTrackStub((*digiIt), detId));
      mplct_.etalut = stub.etaValue();
      mplct_.philut = stub.phiValue();
      mplct_tree_->Fill();

      // debug
      if (verboseMPLCT_){
        cout << "------------------------------------------------------------------------------" << endl;         
        cout << "Event " << mplct_.event << ", detId " << detId << ", MPLCT " << digiIt-range.first << endl;
        cout << "endcap " << mplct_.endcap << ", station " << mplct_.station << ", ring " << mplct_.ring << ", chamber " << mplct_.chamber << endl;
        cout << *digiIt << endl;
        cout << "eta " << mplct_.etalut << ", phi " << mplct_.philut << endl;
      }
    }
  }
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeTFTrackRate(const edm::Event& iEvent, enum tfTrack type)
{
  // CSCTF
  // DTTF
  // RPC
  // RPCf
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeTFTrackRate(const edm::Event& iEvent)
{
  edm::Handle< L1CSCTrackCollection > hl1Tracks;
  iEvent.getByLabel(inputCSCTFTrack_,hl1Tracks);
  const L1CSCTrackCollection* l1Tracks = hl1Tracks.product();

  for (auto trk = l1Tracks->begin(); trk != l1Tracks->end(); trk++) {
    if (trk->first.bx() < minBXCSCTFTrack_ or trk->first.bx() > maxBXCSCTFTrack_) continue;
    const bool endcapOnly(true);
    if (endcapOnly and abs(trk->first.endcap())!=1) continue;
    
    MatchCSCMuL1::TFTRACK myTFTrk;
    myTFTrk.init( &(trk->first) , ptLUT, muScales, muPtScale);

    tftrack_.event = iEvent.id().event();
    tftrack_.bx = trk->first.bx();
    tftrack_.pt = myTFTrk.pt;
    tftrack_.eta = myTFTrk.eta;
    tftrack_.phi = myTFTrk.phi;
    if (tftrack_.pt < 0.001) continue;

    for (auto detUnitIt = trk->second.begin(); detUnitIt != trk->second.end(); detUnitIt++) {
      const CSCDetId& id = (*detUnitIt).first;
      auto range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        if (!((*digiIt).isValid())) continue;
        myTFTrk.trgdigis.push_back(&*digiIt);
        myTFTrk.trgids.push_back(id);
        myTFTrk.trgetaphis.push_back(intersectionEtaPhi(id,(*digiIt).getKeyWG(),(*digiIt).getStrip()));
        myTFTrk.trgstubs.push_back(buildTrackStub((*digiIt),id));

        // stub analysis
        if (id.station()==1 and id.ring()==4) tftrack_.hasME1a = 1;
        if (id.station()==1 and id.ring()==1) tftrack_.hasME1b = 1;
        if (id.station()==1 and id.ring()==2) tftrack_.hasME12 = 1;
        if (id.station()==1 and id.ring()==3) tftrack_.hasME13 = 1; 
        if (id.station()==2 and id.ring()==1) tftrack_.hasME21 = 1;
        if (id.station()==2 and id.ring()==2) tftrack_.hasME22 = 1;
        if (id.station()==3 and id.ring()==1) tftrack_.hasME31 = 1;
        if (id.station()==3 and id.ring()==2) tftrack_.hasME32 = 1;
        if (id.station()==4 and id.ring()==1) tftrack_.hasME41 = 1;
        if (id.station()==4 and id.ring()==2) tftrack_.hasME42 = 1;
      }
    }

    // nstubs
    const int nTrgDigis(myTFTrk.trgdigis.size());
    const int nTrigIds(myTFTrk.trgids.size());
    const int nTrgEtaPhis(myTFTrk.trgetaphis.size());
    const int nTrgStubs(myTFTrk.trgstubs.size());
    
    myTFTrk.nTFStubs = myTFTrk.nStubs(1,1,1,1,1);

    if (min({nTrgDigis, nTrigIds, nTrgEtaPhis, nTrgStubs}) != 
	max({nTrgDigis, nTrigIds, nTrgEtaPhis, nTrgStubs})){
      cout << "ALARM!!! Weird CSC TFTrack! " << endl;
    }

    if (verboseCSCTFTrack_){
      cout << "------------------------------------------------------------------------------" << endl
                << "TFTrack " << trk - l1Tracks->begin() << " information" << endl
                << "bx " << tftrack_.bx << ", pt " << tftrack_.pt << ", eta " << tftrack_.eta << ", phi " << tftrack_.phi << endl
		<< "Summary of endcap hits: nTrgDigis " << nTrgDigis << ", nTrigIds " << nTrigIds 
		<< ", nTrgEtaPhis " << nTrgEtaPhis << ", nTrgStubs " << nTrgStubs << endl
		<< ", nStubs " << myTFTrk.nTFStubs << endl
                << "Station 0: " << endl << "\tME0 " << tftrack_.hasME0 << endl
                << "Station 1: " << endl << "\tME1a " << tftrack_.hasME1a << "\tME1b " << tftrack_.hasME1b << "\tME12 " << tftrack_.hasME12 
                << "\tME13 " << tftrack_.hasME13 << "\tGE11 " << tftrack_.hasGE11 << "\tRE12 " << tftrack_.hasRE12 << "\tRE13 " << tftrack_.hasRE13 << endl
                << "Station 2: " << endl << "\tME21 " << tftrack_.hasME21 << "\tME22 " << tftrack_.hasME22 << "\tGE21S "<< tftrack_.hasGE21S 
                << "\tGE21L "<< tftrack_.hasGE21L << "\tRE22 " << tftrack_.hasRE22 << "\tRE23 " << tftrack_.hasRE23 << endl
                << "Station 3: " << endl << "\tME31 " << tftrack_.hasME31 << "\tME32 " << tftrack_.hasME32 << "\tRE31 " << tftrack_.hasRE31 
                << "\tRE32 " << tftrack_.hasRE32 << "\tRE33 " << tftrack_.hasRE33 << endl
                << "Station 4: " << endl << "\tME41 " << tftrack_.hasME41 << "\tME42 " << tftrack_.hasME42 << "\tRE41 " << tftrack_.hasRE41 
                << "\tRE42 " << tftrack_.hasRE42 << "\tRE43 " << tftrack_.hasRE43 << endl;
    }

    // add the TFTrack to the list
    rtTFTracks_.push_back(myTFTrk);

    tftrack_tree_->Fill();
  }
}

// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeTFCandRate(const edm::Event& iEvent)
{
  edm::Handle< std::vector< L1MuRegionalCand > > hl1TfCands;
  iEvent.getByLabel(inputCSCTFCand_, hl1TfCands);
  const std::vector< L1MuRegionalCand > *l1TfCands = hl1TfCands.product();

  for (auto trk = l1TfCands->begin(); trk != l1TfCands->end(); trk++){
    if ( trk->bx() < minBXCSCTFCand_ or trk->bx() > maxBXCSCTFCand_ ) continue;
    //    const int sign_eta(((trk->eta_packed() & 0x20) == 0) ? 1.:-1);
    MatchCSCMuL1::TFCAND myTFCand;
    myTFCand.init( &*trk , ptLUT, muScales, muPtScale);
    // associate the TFTracks to this TFCand
    for (size_t tt = 0; tt<rtTFTracks_.size(); tt++){
      if (trk->bx()         != rtTFTracks_[tt].l1trk->bx() or
          trk->phi_packed() != rtTFTracks_[tt].phi_packed or
          trk->pt_packed()  != rtTFTracks_[tt].pt_packed or
          trk->eta_packed() != rtTFTracks_[tt].eta_packed) continue;
      myTFCand.tftrack = &(rtTFTracks_[tt]);
      // ids now hold *trigger segments IDs*
      myTFCand.ids = myTFCand.tftrack->trgids;
      myTFCand.nTFStubs =myTFCand.tftrack->trgstubs.size(); 
    }

    // analysis
    if (myTFCand.tftrack == nullptr) continue;
    tfcand_.event = iEvent.id().event();
    tfcand_.bx = trk->bx();
    tfcand_.pt = myTFCand.tftrack->pt;
    tfcand_.eta = myTFCand.tftrack->eta;
    tfcand_.phi = myTFCand.tftrack->phi;
    tfcand_.nStubs = myTFCand.nTFStubs;
    tfcand_.nDetIds = myTFCand.ids.size();
    if (tfcand_.pt < 0.001) continue;
    if (tfcand_.nDetIds <= 0 or tfcand_.nStubs <= 0) continue;

    // stub analysis
    for (auto id : myTFCand.tftrack->trgids){
      if (id.station()==1 and id.ring()==4) tfcand_.hasME1a = 1;
      if (id.station()==1 and id.ring()==1) tfcand_.hasME1b = 1;
      if (id.station()==1 and id.ring()==2) tfcand_.hasME12 = 1;
      if (id.station()==1 and id.ring()==3) tfcand_.hasME13 = 1; 
      if (id.station()==2 and id.ring()==1) tfcand_.hasME21 = 1;
      if (id.station()==2 and id.ring()==2) tfcand_.hasME22 = 1;
      if (id.station()==3 and id.ring()==1) tfcand_.hasME31 = 1;
      if (id.station()==3 and id.ring()==2) tfcand_.hasME32 = 1;
      if (id.station()==4 and id.ring()==1) tfcand_.hasME41 = 1;
      if (id.station()==4 and id.ring()==2) tfcand_.hasME42 = 1;
    }

    if (verboseCSCTFCand_){
      cout << "------------------------------------------------------------------------------" << endl
                << "TFCand " << trk - l1TfCands->begin() << " information" << endl
                << "bx " << tfcand_.bx << ", pt " << tfcand_.pt << ", eta " << tfcand_.eta << ", phi " << tfcand_.phi << endl
                << "Summary of endcap hits: " << tfcand_.nStubs << " stubs in " << tfcand_.nDetIds << " detIds " << endl
                << "Station 0: " << endl << "\tME0 " << tfcand_.hasME0 << endl
                << "Station 1: " << endl << "\tME1a " << tfcand_.hasME1a << "\tME1b " << tfcand_.hasME1b << "\tME12 " << tfcand_.hasME12 
                << "\tME13 " << tfcand_.hasME13 << "\tGE11 " << tfcand_.hasGE11 << "\tRE12 " << tfcand_.hasRE12 << "\tRE13 " << tfcand_.hasRE13 << endl
                << "Station 2: " << endl << "\tME21 " << tfcand_.hasME21 << "\tME22 " << tfcand_.hasME22 << "\tGE21S "<< tfcand_.hasGE21S 
                << "\tGE21L "<< tfcand_.hasGE21L << "\tRE22 " << tfcand_.hasRE22 << "\tRE23 " << tfcand_.hasRE23 << endl
                << "Station 3: " << endl << "\tME31 " << tfcand_.hasME31 << "\tME32 " << tfcand_.hasME32 << "\tRE31 " << tfcand_.hasRE31 
                << "\tRE32 " << tfcand_.hasRE32 << "\tRE33 " << tfcand_.hasRE33 << endl
                << "Station 4: " << endl << "\tME41 " << tfcand_.hasME41 << "\tME42 " << tfcand_.hasME42 << "\tRE41 " << tfcand_.hasRE41 
                 << "\tRE42 " << tfcand_.hasRE42 << "\tRE43 " << tfcand_.hasRE43 << endl;
    }
    
    rtTFCands_.push_back(myTFCand);

    tfcand_tree_->Fill();
  }
}


// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeGMTRegCandRate(const edm::Event& iEvent)
{
  analyzeGMTRegCandRate(iEvent, gmtRegCand::CSCGMT);
  analyzeGMTRegCandRate(iEvent, gmtRegCand::DTGMT);
  analyzeGMTRegCandRate(iEvent, gmtRegCand::RPCbGMT);
  analyzeGMTRegCandRate(iEvent, gmtRegCand::RPCfGMT);
}


// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeGMTRegCandRate(const edm::Event& iEvent, enum gmtRegCand type)
{
  std::vector<L1MuRegionalCand> collection;
  switch (type){
  case CSCGMT:
    collection = l1GmtCSCCands_;
    break;
  case DTGMT:
    collection = l1GmtDTCands_;
    break;
  case RPCbGMT:
    collection = l1GmtRPCbCands_;
    break;
  case RPCfGMT:
    collection = l1GmtRPCfCands_;
    break;
  }
  
  for (auto trk = collection.begin(); trk != collection.end(); trk++) {
    if (trk->bx() < minBXGMTRegCand_ or trk->bx() > maxBXGMTRegCand_) continue;
    //double sign_eta = ( (trk->eta_packed() & 0x20) == 0) ? 1.:-1;
    MatchCSCMuL1::GMTREGCAND myGMTREGCand;
    myGMTREGCand.init( &*trk , muScales, muPtScale);
    myGMTREGCand.tfcand = nullptr;

    if (verboseGMTRegCand_){
      cout << "track " << trk - collection.begin() << endl 
	   << "\tnTFCands " << rtTFCands_.size() << endl;
    }
    for (unsigned i=0; i< rtTFCands_.size(); i++) {
      // check if this TFCand has stubs or ids
      // check if GMTRegCand is compatible 
      if (trk->bx()          != rtTFCands_[i].l1cand->bx()         or
          trk->phi_packed()  != rtTFCands_[i].l1cand->phi_packed() or
          trk->eta_packed()  != rtTFCands_[i].l1cand->eta_packed()) {
	continue;
      }
      myGMTREGCand.tfcand = &(rtTFCands_[i]);
      myGMTREGCand.ids = rtTFCands_[i].ids;
      myGMTREGCand.nTFStubs = rtTFCands_[i].nTFStubs;
      break;
    }
    
    
    // analysis
    if (myGMTREGCand.tfcand == nullptr) continue;
    gmtregcand_.event = iEvent.id().event();
    gmtregcand_.bx = trk->bx();
    gmtregcand_.pt = myGMTREGCand.pt;
    gmtregcand_.eta = myGMTREGCand.eta;
    gmtregcand_.phi = myGMTREGCand.phi;
    gmtregcand_.nStubs = myGMTREGCand.nTFStubs;
    gmtregcand_.nDetIds = myGMTREGCand.ids.size();
    gmtregcand_.quality = trk->quality();
    if (gmtregcand_.pt < 0.001) continue;
    if (gmtregcand_.nDetIds <= 0 or gmtregcand_.nStubs <= 0) continue;

    switch(type) {
    case gmtRegCand::CSCGMT:
      rtGmtRegCscCands_.push_back(myGMTREGCand);
      gmtregcand_.isCSC  = 1;
      if (verboseGMTRegCand_) cout << "isCSC"<<endl;
      break;
    case gmtRegCand::DTGMT:
      rtGmtRegDtCands_.push_back(myGMTREGCand);
      gmtregcand_.isDT   = 1;
      if (verboseGMTRegCand_) cout << "isDT"<<endl;
      break;
    case gmtRegCand::RPCbGMT:
      rtGmtRegRpcbCands_.push_back(myGMTREGCand);
      gmtregcand_.isRPCb = 1;
      if (verboseGMTRegCand_) cout << "isRPCb"<<endl;
      break;
    case gmtRegCand::RPCfGMT:
      rtGmtRegRpcfCands_.push_back(myGMTREGCand);
      gmtregcand_.isRPCf = 1;
      if (verboseGMTRegCand_) cout << "isRPCf"<<endl;
      break;
    }

    // stub analysis
//     if (verboseGMTRegCand_) cout << "Nstubs " <<  myGMTREGCand.ids.size() << endl;
    for (auto id : myGMTREGCand.ids){
      if (id.station()==1 and id.ring()==4) gmtregcand_.hasME1a = 1;
      if (id.station()==1 and id.ring()==1) gmtregcand_.hasME1b = 1;
      if (id.station()==1 and id.ring()==2) gmtregcand_.hasME12 = 1;
      if (id.station()==1 and id.ring()==3) gmtregcand_.hasME13 = 1; 
      if (id.station()==2 and id.ring()==1) gmtregcand_.hasME21 = 1;
      if (id.station()==2 and id.ring()==2) gmtregcand_.hasME22 = 1;
      if (id.station()==3 and id.ring()==1) gmtregcand_.hasME31 = 1;
      if (id.station()==3 and id.ring()==2) gmtregcand_.hasME32 = 1;
      if (id.station()==4 and id.ring()==1) gmtregcand_.hasME41 = 1;
      if (id.station()==4 and id.ring()==2) gmtregcand_.hasME42 = 1;
    }

    if (verboseGMTRegCand_){
      cout << "------------------------------------------------------------------------------" << endl
                << "GMTRegCand " << trk - collection.begin() << " information" << endl
                << "bx " << gmtregcand_.bx << ", pt " << gmtregcand_.pt << ", eta " << gmtregcand_.eta << ", phi " << gmtregcand_.phi << endl
                << "Summary of endcap ids: " << gmtregcand_.nStubs << " stubs in " << gmtregcand_.nDetIds << " detIds " << endl
                << "Station 0: " << endl << "\tME0 " << gmtregcand_.hasME0 << endl
                << "Station 1: " << endl << "\tME1a " << gmtregcand_.hasME1a << "\tME1b " << gmtregcand_.hasME1b << "\tME12 " << gmtregcand_.hasME12 
                << "\tME13 " << gmtregcand_.hasME13 << "\tGE11 " << gmtregcand_.hasGE11 << "\tRE12 " << gmtregcand_.hasRE12 << "\tRE13 " << gmtregcand_.hasRE13 << endl
                << "Station 2: " << endl << "\tME21 " << gmtregcand_.hasME21 << "\tME22 " << gmtregcand_.hasME22 << "\tGE21S "<< gmtregcand_.hasGE21S 
                << "\tGE21L "<< gmtregcand_.hasGE21L << "\tRE22 " << gmtregcand_.hasRE22 << "\tRE23 " << gmtregcand_.hasRE23 << endl
                << "Station 3: " << endl << "\tME31 " << gmtregcand_.hasME31 << "\tME32 " << gmtregcand_.hasME32 << "\tRE31 " << gmtregcand_.hasRE31 
                << "\tRE32 " << gmtregcand_.hasRE32 << "\tRE33 " << gmtregcand_.hasRE33 << endl
                << "Station 4: " << endl << "\tME41 " << gmtregcand_.hasME41 << "\tME42 " << gmtregcand_.hasME42 << "\tRE41 " << gmtregcand_.hasRE41 
                << "\tRE42 " << gmtregcand_.hasRE42 << "\tRE43 " << gmtregcand_.hasRE43 << endl;
    }
    
    gmtregcand_tree_->Fill();
  }
}


// ================================================================================================
void  
GEMCSCTriggerRateTree::analyzeGMTCandRate(const edm::Event& iEvent)
{
  const std::vector<L1MuGMTReadoutRecord> gmt_records(hl1GmtCands_->getRecords());

  for (auto rItr=gmt_records.begin(); rItr!=gmt_records.end(); ++rItr) {
    if (rItr->getBxInEvent() < minBXGMTCand_ or rItr->getBxInEvent() > maxBXGMTCand_) continue;
    
    const std::vector<L1MuGMTExtendedCand> GMTCands(rItr->getGMTCands());    
    const std::vector<L1MuRegionalCand> CSCCands(rItr->getCSCCands());
    const std::vector<L1MuRegionalCand> DTCands(rItr->getDTBXCands());
    const std::vector<L1MuRegionalCand> RPCfCands(rItr->getFwdRPCCands());
    const std::vector<L1MuRegionalCand> RPCbCands(rItr->getBrlRPCCands());

    // check for empty collections
    if (GMTCands.size()==0 and CSCCands.size()==0 and DTCands.size()==0 and RPCfCands.size()==0 and RPCbCands.size()==0) break; 
    
    if (verboseGMTCand_){
      cout << "------------------------------------------------------------------------------" << endl;         
      cout << "GMT records for BX " << rItr->getBxInEvent() << endl
	   << "\tnGMTCands " << GMTCands.size()
	   << "\tnCSCCands " << CSCCands.size()
	   << "\tnDTCands " << DTCands.size()
	   << "\tnRPCfCands " << RPCfCands.size()
	   << "\tnRPCbCands " << RPCbCands.size() << endl;
    }
    // check for empty collections
    if (GMTCands.size()!=0 and (CSCCands.size()==0 and DTCands.size()==0 and 
				RPCfCands.size()==0 and RPCbCands.size()==0)){
      if (verboseGMTCand_) cout<<"Warning: non-empty GMT collection with all empty GMTReg collections!!"<<endl;
    }
    
    for (auto trk = GMTCands.begin(); trk != GMTCands.end(); ++trk){
      if(trk->empty()) continue;
      if (trk->bx() < minBXGMTCand_ or trk->bx() > minBXGMTCand_) continue;

      MatchCSCMuL1::GMTCAND myGMTCand;
      myGMTCand.init( &*trk , muScales, muPtScale);
      if (doSelectEtaForGMTRates_ && myGMTCand.eta<0) continue;
      
      myGMTCand.regcand = nullptr;
      myGMTCand.regcand_rpc = nullptr;
      
      //       float gpt = myGMTCand.pt;
      const float geta(fabs(myGMTCand.eta));
      const bool eta_q(geta > 1.2);
      
      MatchCSCMuL1::GMTREGCAND * gmt_csc = nullptr;
      if (trk->isFwd() && ( trk->isMatchedCand() or !trk->isRPC())) {
        const L1MuRegionalCand rcsc(CSCCands[trk->getDTCSCIndex()]);
        unsigned my_i = 999;
	cout << "number of regional candidates " << rtGmtRegCscCands_.size() << endl;
        for (unsigned i=0; i< rtGmtRegCscCands_.size(); i++) {
	  cout << i << ": regional dataword " << rcsc.getDataWord() << " gmt data word " << rtGmtRegCscCands_[i].l1reg->getDataWord() <<endl;
          if (rcsc.getDataWord()!=rtGmtRegCscCands_[i].l1reg->getDataWord()) continue;
          my_i = i;
          break;
        }
        if (my_i<99) gmt_csc = &rtGmtRegCscCands_[my_i];
	else cout<<"DOES NOT EXIST IN rtGmtRegCscCands_! Should not happen!"<<endl;cout<<"DOES NOT EXIST IN rtGmtRegCscCands_! Should not happen!"<<endl;
	myGMTCand.regcand = gmt_csc;
      }
      
      MatchCSCMuL1::GMTREGCAND * gmt_rpcf = nullptr;
      if (trk->isFwd() && (trk->isMatchedCand() or trk->isRPC())) {
        const L1MuRegionalCand rrpcf = RPCfCands[trk->getRPCIndex()];
        unsigned my_i = 999;
        for (unsigned i=0; i< rtGmtRegRpcfCands_.size(); i++) {
          if (rrpcf.getDataWord()!=rtGmtRegRpcfCands_[i].l1reg->getDataWord()) continue;
          my_i = i;
          break;
        }
        if (my_i<99) gmt_rpcf = &rtGmtRegRpcfCands_[my_i];
        else cout<<"DOES NOT EXIST IN rtGmtRegRpcfCands_! Should not happen!"<<endl;
        myGMTCand.regcand_rpc = gmt_rpcf;
      }
      
      MatchCSCMuL1::GMTREGCAND * gmt_rpcb = nullptr;
      if (!(trk->isFwd()) && (trk->isMatchedCand() or trk->isRPC())){
        const L1MuRegionalCand rrpcb = RPCbCands[trk->getRPCIndex()];
        unsigned my_i = 999;
        for (unsigned i=0; i< rtGmtRegRpcbCands_.size(); i++){
          if (rrpcb.getDataWord()!=rtGmtRegRpcbCands_[i].l1reg->getDataWord()) continue;
          my_i = i;
          break;
        }
        if (my_i<99) gmt_rpcb = &rtGmtRegRpcbCands_[my_i];
        else cout<<"DOES NOT EXIST IN rtGmtRegRpcbCands_! Should not happen!"<<endl;
        myGMTCand.regcand_rpc = gmt_rpcb;
      }
      
      MatchCSCMuL1::GMTREGCAND * gmt_dt = nullptr;
      if (!(trk->isFwd()) && (trk->isMatchedCand() or !(trk->isRPC()))){
        const L1MuRegionalCand rdt = DTCands[trk->getDTCSCIndex()];
        unsigned my_i = 999;
        for (unsigned i=0; i< rtGmtRegDtCands_.size(); i++){
          if (rdt.getDataWord()!=rtGmtRegDtCands_[i].l1reg->getDataWord()) continue;
          my_i = i;
          break;
        }
        if (my_i<99) gmt_dt = &rtGmtRegDtCands_[my_i];
        else cout<<"DOES NOT EXIST IN rtGmtRegDtCands_! Should not happen!"<<endl;
        myGMTCand.regcand = gmt_dt;
      }
      // stub analysis 
      rtGmtCands_.push_back(myGMTCand);
      
      // analysis
      if (myGMTCand.regcand == nullptr and myGMTCand.regcand_rpc == nullptr) continue;
      gmtcand_.event = iEvent.id().event();
      gmtcand_.bx = trk->bx();
      gmtcand_.pt = myGMTCand.pt;
      gmtcand_.eta = myGMTCand.eta;
      gmtcand_.phi = myGMTCand.phi;
      gmtcand_.quality = trk->quality();
      gmtcand_.nStubs = gmt_csc->nTFStubs;
      myGMTCand.ids = gmt_csc->ids;
      gmtcand_.nDetIds = myGMTCand.ids.size();
      if (gmtcand_.pt < 0.001) continue;
      // don't process non-CSC tracks
      if (gmtcand_.nDetIds <= 0 or gmtcand_.nStubs <= 0) continue;

      if (trk->useInSingleMuonTrigger()) gmtcand_.isGoodSingleMuon = 1;
      if (trk->useInSingleMuonTrigger()) gmtcand_.isGoodDoubleMuon = 1;
      auto trgids(myGMTCand.ids);
      // stub analysis
      for (auto id : trgids){
        if (id.station()==1 and id.ring()==4) gmtcand_.hasME1a = 1;
        if (id.station()==1 and id.ring()==1) gmtcand_.hasME1b = 1;
        if (id.station()==1 and id.ring()==2) gmtcand_.hasME12 = 1;
        if (id.station()==1 and id.ring()==3) gmtcand_.hasME13 = 1; 
        if (id.station()==2 and id.ring()==1) gmtcand_.hasME21 = 1;
        if (id.station()==2 and id.ring()==2) gmtcand_.hasME22 = 1;
        if (id.station()==3 and id.ring()==1) gmtcand_.hasME31 = 1;
        if (id.station()==3 and id.ring()==2) gmtcand_.hasME32 = 1;
        if (id.station()==4 and id.ring()==1) gmtcand_.hasME41 = 1;
        if (id.station()==4 and id.ring()==2) gmtcand_.hasME42 = 1;
      }

      const bool isCSC(gmt_csc != nullptr);
      const bool isDT(gmt_dt  != nullptr);
      const bool isRPCf(gmt_rpcf != nullptr);
      const bool isRPCb(gmt_rpcb != nullptr);
      if (isCSC)  gmtcand_.isCSC  = 1;
      if (isDT)   gmtcand_.isDT   = 1;
      if (isRPCb) gmtcand_.isRPCb = 1;
      if (isRPCf) gmtcand_.isRPCf = 1;

      const bool hasCSCCand(isCSC and gmt_csc->tfcand != nullptr);
      const bool hasCSCTrack(isCSC and gmt_csc->tfcand->tftrack != nullptr);
      //      const bool hasCSCL1Track(isCSC and gmt_csc->tfcand->tftrack->l1trk != nullptr);
      const bool hasCSCRegCand(isCSC and gmt_csc->l1reg != nullptr);
      const bool hasCSCCandTrack(hasCSCCand and hasCSCTrack);
      const bool isCSC2s(hasCSCCandTrack and gmt_csc->tfcand->tftrack->nStubs()>=2);
      const bool isCSC3s(hasCSCCandTrack and ((!eta_q and isCSC2s) or (eta_q and gmt_csc->tfcand->tftrack->nStubs()>=3)));
      const bool isCSC2q(hasCSCRegCand and gmt_csc->l1reg->quality()>=2);
      const bool isCSC3q(hasCSCRegCand and ( (!eta_q and isCSC2q) or (eta_q and gmt_csc->l1reg->quality()>=3) ) );
      
      myGMTCand.isCSC = isCSC;
      myGMTCand.isDT = isDT;
      myGMTCand.isRPCf = isRPCf;
      myGMTCand.isRPCb = isRPCb;
      myGMTCand.isCSC2s = isCSC2s;
      myGMTCand.isCSC3s = isCSC3s;
      myGMTCand.isCSC2q = isCSC2q;
      myGMTCand.isCSC3q = isCSC3q;

      if (verboseGMTCand_){
        cout << "------------------------------------------------------------------------------" << endl
                  << "GMTCand " << trk - GMTCands.begin() << " information" << endl
                  << "bx " << gmtcand_.bx << ", pt " << gmtcand_.pt << ", eta " << gmtcand_.eta << ", phi " << gmtcand_.phi << endl
                  << "Summary of endcap hits: " << gmtcand_.nStubs << " stubs in " << gmtcand_.nDetIds << " detIds " << endl
                  << "Station 0: " << endl << "\tME0 " << gmtcand_.hasME0 << endl
                  << "Station 1: " << endl << "\tME1a " << gmtcand_.hasME1a << "\tME1b " << gmtcand_.hasME1b << "\tME12 " << gmtcand_.hasME12 
                  << "\tME13 " << gmtcand_.hasME13 << "\tGE11 " << gmtcand_.hasGE11 << "\tRE12 " << gmtcand_.hasRE12 << "\tRE13 " << gmtcand_.hasRE13 << endl
                  << "Station 2: " << endl << "\tME21 " << gmtcand_.hasME21 << "\tME22 " << gmtcand_.hasME22 << "\tGE21S "<< gmtcand_.hasGE21S 
                  << "\tGE21L "<< gmtcand_.hasGE21L << "\tRE22 " << gmtcand_.hasRE22 << "\tRE23 " << gmtcand_.hasRE23 << endl
                  << "Station 3: " << endl << "\tME31 " << gmtcand_.hasME31 << "\tME32 " << gmtcand_.hasME32 << "\tRE31 " << gmtcand_.hasRE31 
                  << "\tRE32 " << gmtcand_.hasRE32 << "\tRE33 " << gmtcand_.hasRE33 << endl
                  << "Station 4: " << endl << "\tME41 " << gmtcand_.hasME41 << "\tME42 " << gmtcand_.hasME42 << "\tRE41 " << gmtcand_.hasRE41 
                  << "\tRE42 " << gmtcand_.hasRE42 << "\tRE43 " << gmtcand_.hasRE43 << endl;
       }      
      
      rtGmtCands_.push_back(myGMTCand);
      gmtcand_tree_->Fill();
    }
  }
}


// ================================================================================================
void 
GEMCSCTriggerRateTree::runCSCTFSP(const CSCCorrelatedLCTDigiCollection* mplcts, const L1MuDTChambPhContainer* dttrig)
   //, L1CSCTrackCollection*, CSCTriggerContainer<csctf::TrackStub>*)
{
// Just run it for the sake of its debug printout, do not return any results

  // Create csctf::TrackStubs collection from MPC LCTs
  CSCTriggerContainer<csctf::TrackStub> stub_list;
  CSCCorrelatedLCTDigiCollection::DigiRangeIterator Citer;
  for(Citer = mplcts->begin(); Citer != mplcts->end(); Citer++)
  {
    CSCCorrelatedLCTDigiCollection::const_iterator Diter = (*Citer).second.first;
    CSCCorrelatedLCTDigiCollection::const_iterator Dend = (*Citer).second.second;
    for(; Diter != Dend; Diter++)
    {
      csctf::TrackStub theStub((*Diter),(*Citer).first);
      stub_list.push_back(theStub);
    }
  }
  
  // Now we append the track stubs the the DT Sector Collector
  // after processing from the DT Receiver.
  CSCTriggerContainer<csctf::TrackStub> dtstubs = my_dtrc->process(dttrig);
  stub_list.push_many(dtstubs);
  
  //for(int e=0; e<2; e++) for (int s=0; s<6; s++) {
  int e=0;
  for (int s=0; s<6; s++) 
  {
    CSCTriggerContainer<csctf::TrackStub> current_e_s = stub_list.get(e+1, s+1);
    if (current_e_s.get().size()>0) 
    {
      cout<<"sector "<<s+1<<":"<<endl<<endl;
      my_SPs[e][s]->run(current_e_s);
    }
  }
}

// ================================================================================================
// Returns chamber type (0-9) according to the station and ring number
int 
GEMCSCTriggerRateTree::getCSCType(CSCDetId &id) 
{
  int type = -999;

  if (id.station() == 1) {
    type = (id.triggerCscId()-1)/3;
    if (id.ring() == 4) {
      type = 3;
    }
  }
  else { // stations 2-4
    type = 3 + id.ring() + 2*(id.station()-2);
  }
  assert(type >= 0 && type < CSC_TYPES); // include ME4/2
  return type;
}

// ================================================================================================
int
GEMCSCTriggerRateTree::isME11(int t)
{
  if (t==0 or t==3) return CSC_TYPES;
  return 0;
}

// ================================================================================================
// Returns chamber type (0-9) according to CSCChamberSpecs type
// 1..10 -> 1/a, 1/b, 1/2, 1/3, 2/1...
int
GEMCSCTriggerRateTree::getCSCSpecsType(CSCDetId &id)
{
  return cscGeometry->chamber(id)->specs()->chamberType();
}

// ================================================================================================
int
GEMCSCTriggerRateTree::cscTriggerSubsector(CSCDetId &id)
{
  if(id.station() != 1) return 0; // only station one has subsectors
  int chamber(id.chamber());
  switch(chamber) // first make things easier to deal with
  {
    case 1:
      chamber = 36;
      break;
    case 2:
      chamber = 35;
      break;
    default:
      chamber -= 2;
  }
  chamber = ((chamber-1)%6) + 1; // renumber all chambers to 1-6
  return ((chamber-1) / 3) + 1; // [1,3] -> 1 , [4,6]->2
}


// ================================================================================================
std::pair<float, float> 
GEMCSCTriggerRateTree::intersectionEtaPhi(CSCDetId id, int wg, int hs)
{
  const CSCDetId layerId(id.endcap(), id.station(), id.ring(), id.chamber(), CSCConstants::KEY_CLCT_LAYER);
  const CSCLayer* csclayer(cscGeometry->layer(layerId));
  const CSCLayerGeometry* layer_geo(csclayer->geometry());
    
  // LCT::getKeyWG() starts from 0
  const float wire(layer_geo->middleWireOfGroup(wg + 1));

  // half-strip to strip
  // note that LCT's HS starts from 0, but in geometry strips start from 1
  const float fractional_strip(0.5 * (hs + 1) - 0.25);
  const LocalPoint csc_intersect(layer_geo->intersectionOfStripAndWire(fractional_strip, wire));
  const GlobalPoint csc_gp(cscGeometry->idToDet(layerId)->surface().toGlobal(csc_intersect));
  return std::make_pair(csc_gp.eta(), csc_gp.phi());
}

// ================================================================================================
csctf::TrackStub 
GEMCSCTriggerRateTree::buildTrackStub(const CSCCorrelatedLCTDigi &d, CSCDetId id)
{
  const unsigned fpga((id.station() == 1) ? CSCTriggerNumbering::triggerSubSectorFromLabels(id) - 1 : id.station());
  const CSCSectorReceiverLUT* srLUT(srLUTs_[fpga][id.triggerSector()-1][id.endcap()-1]);
  const unsigned cscid(CSCTriggerNumbering::triggerCscIdFromLabels(id));
  const unsigned cscid_special(id.station()==1 and id.ring()==4 ? cscid + 9 : cscid);
  const lclphidat lclPhi(srLUT->localPhi(d.getStrip(), d.getPattern(), d.getQuality(), d.getBend()));
  const gblphidat gblPhi(srLUT->globalPhiME(lclPhi.phi_local, d.getKeyWG(), cscid_special));
  const gbletadat gblEta(srLUT->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, d.getKeyWG(), cscid));
  return csctf::TrackStub(d, id, gblPhi.global_phi, gblEta.global_eta);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMCSCTriggerRateTree);
