// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "L1Trigger/CSCTrackFinder/test/src/RunSRLUTs.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"//
#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
//#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"

//#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
//#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"
#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "GEMCode/GEMValidation/interface/GEMCSCdphi_LUT.h"
#include "GEMCode/GEMValidation/interface/TFTrack.h" 

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;
//using namespace L1TMuon;

class L1MuonTTriggerRate : public edm::EDAnalyzer {
public:
  explicit L1MuonTTriggerRate(const edm::ParameterSet&);
  ~L1MuonTTriggerRate();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;
  //  csctf_analysis::RunSRLUTs* runSRLUTs;
  
  double min_pt;
  double max_pt;
  double min_aEta;
  double max_aEta;
  float ntotalEvents = 0.0;
  /*const float ptscale[33] = { 
    -1.,   0.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
    4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,  
    16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0, 
    50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 1.E6 };
  const float xptscale[30] = { 
    2.0,   2.5,   3.0,   3.5,   4.0,
    4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,  
    16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0, 
    50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 200.0};
  unsigned maxPTbins = 33;
*/

private:
  void init();

private:
  TTree* evtree;
//sim
  int event;

  float simpt, simeta, simphi;
  int charge;
  int endcap;
//TFTrack
  float pt, eta, phi;
  unsigned int quality_packed, pt_packed, eta_packed, phi_packed;
  unsigned int chargesign;
  unsigned int rank;
  unsigned int nstubs;
  unsigned int deltaphi12, deltaphi23; 
  bool hasME1,hasME2;
  float GE11dPhi,GE21dPhi;
  bool passGE11,passGE21;

//need to drop this
  enum etabins{eta_all, eta_me1, eta_me2, eta_me3, netabins};
  enum ptbins{pt_all, pt_15, pt_20, nptbins};
  enum stubbins{stub_2, stub_3, stub_4, nstubbins};
  enum MEbins{ME_all, ME_1, GE_1, passGE_1, ME_2, GE_2, passGE_2, GE_12, passGE1_GE2, passGE_12, passGE, nMEbins};
  TH1F* h_L1CSCTrack_pt[netabins][nptbins][nstubbins][nMEbins];
  TH1F* h_L1CSCTrack_eta[nptbins][nstubbins][nMEbins];
  TH1F* h_L1CSCTrack_phi[netabins][nptbins][nstubbins][nMEbins];

  TH1F* h_nStation;

};

L1MuonTTriggerRate::L1MuonTTriggerRate(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  //  runSRLUTs = new csctf_analysis::RunSRLUTs();
  min_pt = iConfig.getParameter<double>("minPt");
  max_pt = iConfig.getParameter<double>("maxPt");
  min_aEta = iConfig.getParameter<double>("minEta");
  max_aEta = iConfig.getParameter<double>("maxEta");

}

L1MuonTTriggerRate::~L1MuonTTriggerRate()
{
}

void L1MuonTTriggerRate::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<L1CSCTrackCollection> l1csctracks;
  iEvent.getByLabel("simCsctfTrackDigis",l1csctracks);
  //edm::Handle<TriggerPrimitiveCollection> trigPrims;
  //iEvent.getByLabel("L1TMuonTriggerPrimitives",trigPrims);
  edm::Handle<CSCCorrelatedLCTDigiCollection> lcts;
  iEvent.getByLabel("simCscTriggerPrimitiveDigis","MPCSORTED", lcts);

  edm::ESHandle<L1MuTriggerScales> muScalesHd_;
  edm::ESHandle<L1MuTriggerPtScale> muPtScaleHd_;
  //const L1MuTriggerScales* muScales_;
  //const L1MuTriggerPtScale* muPtScale_;
  try {
    iSetup.get<L1MuTriggerScalesRcd>().get(muScalesHd_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerScalesRcd>& e) {
    LogDebug("TrackMatcher") << "+++ Info: L1MuTriggerScalesRcd is unavailable. +++\n";
  }

  try {
    iSetup.get<L1MuTriggerPtScaleRcd>().get(muPtScaleHd_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerPtScaleRcd>& e) {
    LogDebug("TrackMatcher") << "+++ Info: L1MuTriggerPtScaleRcd is unavailable. +++\n";
  }

  //CSCTFPtLUT* ptLUT_ = new CSCTFPtLUT(ptLUTset_, muScalesHd_.product(), muPtScaleHd_.product());

  ++ntotalEvents;

  // lcts used in tracks
  L1CSCTrackCollection::const_iterator tmp_trk = l1csctracks->begin();
  for(; tmp_trk != l1csctracks->end(); tmp_trk++){
    //    float GEMdPhi = -99;

    TFTrack *track = new TFTrack(&tmp_trk->first,&tmp_trk->second);
    track->init(muScalesHd_, muPtScaleHd_);

    pt = track->pt();
    eta = track->eta();
    phi = track->phi();
    pt_packed = track->ptPacked();
    eta_packed = track->etaPacked();
    phi_packed = track->phiPacked();
    quality_packed = track->qPacked();
    chargesign = track->chargesign();
    deltaphi12 = track->dPhi12();
    deltaphi23 = track->dPhi23();
    hasME1 = track->hasStubEndcap(1);
    hasME2 = track->hasStubEndcap(2);
    nstubs = track->nStubs();
    

    //std::cout <<" charge sign " << sign << std::endl;
   unsigned int lct1 = 999;
   auto me1b(track->digiInME(1,1));
   auto me1a(track->digiInME(1,4));
   if (me1a != 999) lct1 = me1a;
   if (me1b != 999) lct1 = me1b;

   if (lct1 < (track->getTriggerDigis()).size()) 
   {
       passGE11 = track->passDPhicutTFTrack(1);
       GE11dPhi = ((track->getTriggerDigis()).at(lct1))->getGEMDPhi();
    }
   unsigned int lct2 = 999;
   lct2 = track->digiInME(2,1);
   if (lct2 < (track->getTriggerDigis()).size()) 
   {
       passGE21 = track->passDPhicutTFTrack(2);
       GE21dPhi = ((track->getTriggerDigis()).at(lct2))->getGEMDPhi();
    }
 
    evtree->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void L1MuonTTriggerRate::beginJob()
{
   evtree = fs->make<TTree>("evtree", "evtree");

   evtree->Branch("nEvents", &ntotalEvents);

   evtree->Branch("eta", &eta);
   evtree->Branch("phi", &phi);
   evtree->Branch("pt", &pt);
   evtree->Branch("quality_packed",&quality_packed);
   evtree->Branch("rank",&rank);
   evtree->Branch("pt_packed",&pt_packed);
   evtree->Branch("eta_packed",&eta_packed);
   evtree->Branch("phi_packed",&phi_packed);
   evtree->Branch("chargesign",&chargesign);
   evtree->Branch("deltaphi12",&deltaphi12);
   evtree->Branch("deltaphi23",&deltaphi23);
   evtree->Branch("hasME1",&hasME1);
   evtree->Branch("hasME2",&hasME2);
   evtree->Branch("dphiGE11",&GE11dPhi);
   evtree->Branch("dphiGE21",&GE21dPhi);
   evtree->Branch("passGE11",&passGE11);
   evtree->Branch("passGE21",&passGE21);
   evtree->Branch("nstubs",&nstubs);
   

}

//--------------- init the tree ------------------------------
void L1MuonTTriggerRate::init(){
  
  eta = -99;
  phi = -99;
  pt = 0;
  quality_packed = 999;
  rank = 999;
  pt_packed = 999;
  phi_packed = 999;
  eta_packed =999;
  chargesign = 999;
  deltaphi12 = 999;
  deltaphi23 = 999;
  hasME1 = false;
  hasME2 = false;
  passGE11 = false;
  passGE21 = false;
  GE11dPhi = -999;


}


// ------------ method called once each job just after ending the event loop  ------------
void L1MuonTTriggerRate::endJob() 
{
  cout << "L1MuonTTriggerRate::endJob ntotalEvents " << ntotalEvents<<endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MuonTTriggerRate);
