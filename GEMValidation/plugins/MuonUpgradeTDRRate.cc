// system include files
#include <memory>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TDatime.h"
#include "TTimeStamp.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "L1Trigger/CSCTrackFinder/test/src/TFTrack.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
// #include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
// #include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
// #include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackCand.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackContainer.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "L1Trigger/DTTrackFinder/interface/L1MuDTTrack.h"
#include "L1Trigger/DTTrackFinder/src/L1MuDTTrackSegPhi.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiL1Linkfwd.h"
#include "DataFormats/RPCDigi/interface/RPCDigiL1Link.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>
#include <L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h>
#include "GEMCode/GEMValidation/plugins/CSCStubPatterns.h"

#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCoPadDigiCollection.h"
//
// class declaration
//

const Int_t kMaxL1Mu = 50;
const Int_t kMaxDTTF = 50;
const Int_t kMaxCSCTF = 50;
const Int_t kMaxRPCb = 50;
const Int_t kMaxRPCf = 50;
const Int_t kMaxGEM = 50;
const Int_t kMaxL1Tk = 500;

typedef std::vector<CSCComparatorDigi> CSCComparatorDigiContainer;
typedef std::vector<std::pair<CSCDetId, CSCComparatorDigiContainer> > CSCComparatorDigiContainerIds;

typedef std::pair<CSCDetId, CSCCorrelatedLCTDigi> CSCCorrelatedLCTDigiId;
typedef std::vector<CSCCorrelatedLCTDigi> CSCCorrelatedLCTDigiContainer;
typedef std::pair<CSCDetId, CSCCorrelatedLCTDigiContainer> CSCCorrelatedLCTDigiContainerId;
typedef std::vector<std::pair<CSCDetId, CSCCorrelatedLCTDigiContainer> > CSCCorrelatedLCTDigiContainerIds;

typedef std::pair<GEMDetId, GEMCoPadDigi> GEMCoPadDigiId;
typedef std::vector<GEMCoPadDigi> GEMCoPadDigiContainer;
typedef std::pair<GEMDetId, GEMCoPadDigiContainer> GEMCoPadDigiContainerId;
typedef std::vector<std::pair<GEMDetId, GEMCoPadDigiContainer> > GEMCoPadDigiContainerIds;

typedef std::pair<GEMDetId, GEMPadDigi> GEMPadDigiId;
typedef std::vector<GEMPadDigi> GEMPadDigiContainer;
typedef std::pair<GEMDetId, GEMPadDigiContainer> GEMPadDigiContainerId;
typedef std::vector<std::pair<GEMDetId, GEMPadDigiContainer> > GEMPadDigiContainerIds;

struct MyEvent
{
  Int_t lumi, run, event;

  Float_t beamSpot_x;
  Float_t beamSpot_y;
  Float_t beamSpot_z;

  Int_t nL1Mu;
  Int_t nL1Tk;
  Float_t L1Mu_pt[kMaxL1Mu], L1Mu_eta[kMaxL1Mu], L1Mu_phi[kMaxL1Mu];
  Int_t L1Mu_charge[kMaxL1Mu], L1Mu_bx[kMaxL1Mu];
  Int_t L1Mu_quality[kMaxL1Mu];
  // Float_t L1Mu_L1Tk_dR_corr[kMaxL1Mu];
  // Float_t L1Mu_L1Tk_pt_corr[kMaxL1Mu];
  // Float_t L1Mu_L1Tk_dR_prop[kMaxL1Mu];
  // Float_t L1Mu_L1Tk_dR_prop_true[kMaxL1Mu];
  // Float_t L1Mu_L1Tk_pt_prop[kMaxL1Mu];
  // Float_t L1Tk_pt[kMaxL1Tk], L1Tk_eta[kMaxL1Tk], L1Tk_phi[kMaxL1Tk];
  // Float_t L1Tk_eta_prop[kMaxL1Tk], L1Tk_phi_prop[kMaxL1Tk];
  // Float_t L1Tk_deta_prop[kMaxL1Tk], L1Tk_dphi_prop[kMaxL1Tk], L1Tk_dR_prop[kMaxL1Tk];
  // Float_t L1Tk_eta_corr[kMaxL1Tk], L1Tk_phi_corr[kMaxL1Tk];
  // Float_t L1Tk_deta_corr[kMaxL1Tk], L1Tk_dphi_corr[kMaxL1Tk], L1Tk_dR_corr[kMaxL1Tk];

  // Matching the L1Mu to DTTF
  Int_t nDTTF;
  Int_t L1Mu_DTTF_index[kMaxL1Mu];
  Float_t DTTF_pt[kMaxDTTF], DTTF_eta[kMaxDTTF], DTTF_phi[kMaxDTTF];
  Int_t DTTF_bx[kMaxDTTF], DTTF_nStubs[kMaxDTTF], DTTF_quality[kMaxDTTF];
  Float_t DTTF_phi1[kMaxDTTF], DTTF_phi2[kMaxDTTF], DTTF_phi3[kMaxDTTF], DTTF_phi4[kMaxDTTF];
  Float_t DTTF_phib1[kMaxDTTF], DTTF_phib2[kMaxDTTF], DTTF_phib3[kMaxDTTF], DTTF_phib4[kMaxDTTF];
  Int_t DTTF_quality1[kMaxDTTF], DTTF_quality2[kMaxDTTF], DTTF_quality3[kMaxDTTF], DTTF_quality4[kMaxDTTF];
  Int_t DTTF_bx1[kMaxDTTF], DTTF_bx2[kMaxDTTF], DTTF_bx3[kMaxDTTF], DTTF_bx4[kMaxDTTF];
  Int_t DTTF_wh1[kMaxDTTF], DTTF_wh2[kMaxDTTF], DTTF_wh3[kMaxDTTF], DTTF_wh4[kMaxDTTF];
  Int_t DTTF_se1[kMaxDTTF], DTTF_se2[kMaxDTTF], DTTF_se3[kMaxDTTF], DTTF_se4[kMaxDTTF];
  Int_t DTTF_st1[kMaxDTTF], DTTF_st2[kMaxDTTF], DTTF_st3[kMaxDTTF], DTTF_st4[kMaxDTTF];

  // Matching the L1Mu to CSCTF  
  Int_t nCSCTF;
  Int_t L1Mu_CSCTF_index[kMaxL1Mu];
  Float_t CSCTF_pt[kMaxCSCTF], CSCTF_eta[kMaxCSCTF], CSCTF_phi[kMaxCSCTF];
  Int_t CSCTF_bx[kMaxCSCTF], CSCTF_nStubs[kMaxCSCTF], CSCTF_quality[kMaxCSCTF];

  Int_t CSCTF_st1[kMaxCSCTF], CSCTF_ri1[kMaxCSCTF], CSCTF_ch1[kMaxCSCTF], CSCTF_en1[kMaxCSCTF];
  Int_t CSCTF_trk1[kMaxCSCTF], CSCTF_quality1[kMaxCSCTF], CSCTF_wg1[kMaxCSCTF], CSCTF_hs1[kMaxCSCTF]; 
  Int_t CSCTF_pat1[kMaxCSCTF], CSCTF_bend1[kMaxCSCTF], CSCTF_bx1[kMaxCSCTF], CSCTF_clctpat1[kMaxCSCTF];
  
  Int_t CSCTF_st2[kMaxCSCTF], CSCTF_ri2[kMaxCSCTF], CSCTF_ch2[kMaxCSCTF], CSCTF_en2[kMaxCSCTF];
  Int_t CSCTF_trk2[kMaxCSCTF], CSCTF_quality2[kMaxCSCTF], CSCTF_wg2[kMaxCSCTF], CSCTF_hs2[kMaxCSCTF]; 
  Int_t CSCTF_pat2[kMaxCSCTF], CSCTF_bend2[kMaxCSCTF], CSCTF_bx2[kMaxCSCTF], CSCTF_clctpat2[kMaxCSCTF];

  Int_t CSCTF_st3[kMaxCSCTF], CSCTF_ri3[kMaxCSCTF], CSCTF_ch3[kMaxCSCTF], CSCTF_en3[kMaxCSCTF];
  Int_t CSCTF_trk3[kMaxCSCTF], CSCTF_quality3[kMaxCSCTF], CSCTF_wg3[kMaxCSCTF], CSCTF_hs3[kMaxCSCTF]; 
  Int_t CSCTF_pat3[kMaxCSCTF], CSCTF_bend3[kMaxCSCTF], CSCTF_bx3[kMaxCSCTF], CSCTF_clctpat3[kMaxCSCTF];

  Int_t CSCTF_st4[kMaxCSCTF], CSCTF_ri4[kMaxCSCTF], CSCTF_ch4[kMaxCSCTF], CSCTF_en4[kMaxCSCTF];
  Int_t CSCTF_trk4[kMaxCSCTF], CSCTF_quality4[kMaxCSCTF], CSCTF_wg4[kMaxCSCTF], CSCTF_hs4[kMaxCSCTF]; 
  Int_t CSCTF_pat4[kMaxCSCTF], CSCTF_bend4[kMaxCSCTF], CSCTF_bx4[kMaxCSCTF], CSCTF_clctpat4[kMaxCSCTF];
  
  Int_t CSCTF_id1[kMaxCSCTF], CSCTF_id2[kMaxCSCTF], CSCTF_id3[kMaxCSCTF], CSCTF_id4[kMaxCSCTF];

  Int_t CSCTF_val1[kMaxCSCTF], CSCTF_val2[kMaxCSCTF], CSCTF_val3[kMaxCSCTF], CSCTF_val4[kMaxCSCTF];
  Float_t CSCTF_phi1[kMaxCSCTF], CSCTF_phi2[kMaxCSCTF], CSCTF_phi3[kMaxCSCTF], CSCTF_phi4[kMaxCSCTF];
  Float_t CSCTF_eta1[kMaxCSCTF], CSCTF_eta2[kMaxCSCTF], CSCTF_eta3[kMaxCSCTF], CSCTF_eta4[kMaxCSCTF];
  Float_t CSCTF_phib1[kMaxCSCTF], CSCTF_phib2[kMaxCSCTF], CSCTF_phib3[kMaxCSCTF], CSCTF_phib4[kMaxCSCTF];

  Float_t CSCTF_gemdphi1[kMaxCSCTF], CSCTF_gemdphi2[kMaxCSCTF];
  Float_t CSCTF_R1[kMaxCSCTF], CSCTF_R2[kMaxCSCTF], CSCTF_R3[kMaxCSCTF], CSCTF_R4[kMaxCSCTF];
  Float_t CSCTF_x1[kMaxCSCTF], CSCTF_x2[kMaxCSCTF], CSCTF_x3[kMaxCSCTF], CSCTF_x4[kMaxCSCTF];
  Float_t CSCTF_y1[kMaxCSCTF], CSCTF_y2[kMaxCSCTF], CSCTF_y3[kMaxCSCTF], CSCTF_y4[kMaxCSCTF];
  Float_t CSCTF_z1[kMaxCSCTF], CSCTF_z2[kMaxCSCTF], CSCTF_z3[kMaxCSCTF], CSCTF_z4[kMaxCSCTF];
  
  // fitted positions -- at key layer
  Float_t CSCTF_fit_phi1[kMaxCSCTF], CSCTF_fit_phi2[kMaxCSCTF], CSCTF_fit_phi3[kMaxCSCTF], CSCTF_fit_phi4[kMaxCSCTF];
  Float_t CSCTF_fit_R1[kMaxCSCTF], CSCTF_fit_R2[kMaxCSCTF], CSCTF_fit_R3[kMaxCSCTF], CSCTF_fit_R4[kMaxCSCTF];
  Float_t CSCTF_fit_x1[kMaxCSCTF], CSCTF_fit_x2[kMaxCSCTF], CSCTF_fit_x3[kMaxCSCTF], CSCTF_fit_x4[kMaxCSCTF];
  Float_t CSCTF_fit_y1[kMaxCSCTF], CSCTF_fit_y2[kMaxCSCTF], CSCTF_fit_y3[kMaxCSCTF], CSCTF_fit_y4[kMaxCSCTF];
  Float_t CSCTF_fit_z1[kMaxCSCTF], CSCTF_fit_z2[kMaxCSCTF], CSCTF_fit_z3[kMaxCSCTF], CSCTF_fit_z4[kMaxCSCTF];

  //fitted directions
  Float_t CSCTF_fit_dphi1[kMaxCSCTF], CSCTF_fit_dphi2[kMaxCSCTF], CSCTF_fit_dphi3[kMaxCSCTF], CSCTF_fit_dphi4[kMaxCSCTF];

  // sim positions - at key layer
  Float_t CSCTF_fitline_x1[kMaxCSCTF], CSCTF_fitline_x2[kMaxCSCTF], CSCTF_fitline_x3[kMaxCSCTF], CSCTF_fitline_x4[kMaxCSCTF];
  Float_t CSCTF_fitline_y1[kMaxCSCTF], CSCTF_fitline_y2[kMaxCSCTF], CSCTF_fitline_y3[kMaxCSCTF], CSCTF_fitline_y4[kMaxCSCTF];
  Float_t CSCTF_fitline_z1[kMaxCSCTF], CSCTF_fitline_z2[kMaxCSCTF], CSCTF_fitline_z3[kMaxCSCTF], CSCTF_fitline_z4[kMaxCSCTF];
  

  // Matching the L1Mu to RPCb  
  Int_t nRPCb;
  Int_t L1Mu_RPCb_index[kMaxL1Mu];
  Float_t RPCb_pt[kMaxRPCb], RPCb_eta[kMaxRPCb], RPCb_phi[kMaxRPCb];
  Int_t RPCb_bx[kMaxRPCb], RPCb_nStubs[kMaxRPCb], RPCb_quality[kMaxRPCb];

  Float_t RPCb_phi1[kMaxRPCb], RPCb_phi2[kMaxRPCb], RPCb_phi3[kMaxRPCb];
  Float_t RPCb_phi4[kMaxRPCb], RPCb_phi5[kMaxRPCb], RPCb_phi6[kMaxRPCb];

  Int_t RPCb_bx1[kMaxRPCb], RPCb_strip1[kMaxRPCb];
  Int_t RPCb_re1[kMaxRPCb], RPCb_ri1[kMaxRPCb], RPCb_st1[kMaxRPCb], RPCb_se1[kMaxRPCb];
  Int_t RPCb_la1[kMaxRPCb], RPCb_su1[kMaxRPCb], RPCb_ro1[kMaxRPCb];
 
  Int_t RPCb_bx2[kMaxRPCb], RPCb_strip2[kMaxRPCb];
  Int_t RPCb_re2[kMaxRPCb], RPCb_ri2[kMaxRPCb], RPCb_st2[kMaxRPCb], RPCb_se2[kMaxRPCb];
  Int_t RPCb_la2[kMaxRPCb], RPCb_su2[kMaxRPCb], RPCb_ro2[kMaxRPCb];

  Int_t RPCb_bx3[kMaxRPCb], RPCb_strip3[kMaxRPCb];
  Int_t RPCb_re3[kMaxRPCb], RPCb_ri3[kMaxRPCb], RPCb_st3[kMaxRPCb], RPCb_se3[kMaxRPCb];
  Int_t RPCb_la3[kMaxRPCb], RPCb_su3[kMaxRPCb], RPCb_ro3[kMaxRPCb];

  Int_t RPCb_bx4[kMaxRPCb], RPCb_strip4[kMaxRPCb];
  Int_t RPCb_re4[kMaxRPCb], RPCb_ri4[kMaxRPCb], RPCb_st4[kMaxRPCb], RPCb_se4[kMaxRPCb];
  Int_t RPCb_la4[kMaxRPCb], RPCb_su4[kMaxRPCb], RPCb_ro4[kMaxRPCb];

  Int_t RPCb_bx5[kMaxRPCb], RPCb_strip5[kMaxRPCb];
  Int_t RPCb_re5[kMaxRPCb], RPCb_ri5[kMaxRPCb], RPCb_st5[kMaxRPCb], RPCb_se5[kMaxRPCb];
  Int_t RPCb_la5[kMaxRPCb], RPCb_su5[kMaxRPCb], RPCb_ro5[kMaxRPCb];

  Int_t RPCb_bx6[kMaxRPCb], RPCb_strip6[kMaxRPCb];
  Int_t RPCb_re6[kMaxRPCb], RPCb_ri6[kMaxRPCb], RPCb_st6[kMaxRPCb], RPCb_se6[kMaxRPCb];
  Int_t RPCb_la6[kMaxRPCb], RPCb_su6[kMaxRPCb], RPCb_ro6[kMaxRPCb];


  // Matching the L1Mu to RPCf  
  Int_t nRPCf;
  Int_t L1Mu_RPCf_index[kMaxL1Mu];
  Float_t RPCf_pt[kMaxRPCf], RPCf_eta[kMaxRPCf], RPCf_phi[kMaxRPCf];
  Int_t RPCf_bx[kMaxRPCf], RPCf_nStubs[kMaxRPCf], RPCf_quality[kMaxRPCf];

  Float_t RPCf_phi1[kMaxRPCf], RPCf_phi2[kMaxRPCf], RPCf_phi3[kMaxRPCf];
  Float_t RPCf_phi4[kMaxRPCf], RPCf_phi5[kMaxRPCf], RPCf_phi6[kMaxRPCf];

  Int_t RPCf_bx1[kMaxRPCf], RPCf_strip1[kMaxRPCf];
  Int_t RPCf_re1[kMaxRPCf], RPCf_ri1[kMaxRPCf], RPCf_st1[kMaxRPCf], RPCf_se1[kMaxRPCf];
  Int_t RPCf_la1[kMaxRPCf], RPCf_su1[kMaxRPCf], RPCf_ro1[kMaxRPCf];
 
  Int_t RPCf_bx2[kMaxRPCf], RPCf_strip2[kMaxRPCf];
  Int_t RPCf_re2[kMaxRPCf], RPCf_ri2[kMaxRPCf], RPCf_st2[kMaxRPCf], RPCf_se2[kMaxRPCf];
  Int_t RPCf_la2[kMaxRPCf], RPCf_su2[kMaxRPCf], RPCf_ro2[kMaxRPCf];

  Int_t RPCf_bx3[kMaxRPCf], RPCf_strip3[kMaxRPCf];
  Int_t RPCf_re3[kMaxRPCf], RPCf_ri3[kMaxRPCf], RPCf_st3[kMaxRPCf], RPCf_se3[kMaxRPCf];
  Int_t RPCf_la3[kMaxRPCf], RPCf_su3[kMaxRPCf], RPCf_ro3[kMaxRPCf];

  Int_t RPCf_bx4[kMaxRPCf], RPCf_strip4[kMaxRPCf];
  Int_t RPCf_re4[kMaxRPCf], RPCf_ri4[kMaxRPCf], RPCf_st4[kMaxRPCf], RPCf_se4[kMaxRPCf];
  Int_t RPCf_la4[kMaxRPCf], RPCf_su4[kMaxRPCf], RPCf_ro4[kMaxRPCf];

  Int_t RPCf_bx5[kMaxRPCf], RPCf_strip5[kMaxRPCf];
  Int_t RPCf_re5[kMaxRPCf], RPCf_ri5[kMaxRPCf], RPCf_st5[kMaxRPCf], RPCf_se5[kMaxRPCf];
  Int_t RPCf_la5[kMaxRPCf], RPCf_su5[kMaxRPCf], RPCf_ro5[kMaxRPCf];

  Int_t RPCf_bx6[kMaxRPCf], RPCf_strip6[kMaxRPCf];
  Int_t RPCf_re6[kMaxRPCf], RPCf_ri6[kMaxRPCf], RPCf_st6[kMaxRPCf], RPCf_se6[kMaxRPCf];
  Int_t RPCf_la6[kMaxRPCf], RPCf_su6[kMaxRPCf], RPCf_ro6[kMaxRPCf];

  // Matching the SIM Mu to GEM pad (really no other way to do this)
  Int_t nGEM;
  Float_t GE11_phi_L1[kMaxGEM], GE11_phi_L2[kMaxGEM], GE21_phi_L1[kMaxGEM], GE21_phi_L2[kMaxGEM];
  Int_t GE11_bx_L1[kMaxGEM], GE11_bx_L2[kMaxGEM], GE21_bx_L1[kMaxGEM], GE21_bx_L2[kMaxGEM];
  Int_t GE11_ch_L1[kMaxGEM], GE11_ch_L2[kMaxGEM], GE21_ch_L1[kMaxGEM], GE21_ch_L2[kMaxGEM];
  Float_t GE11_z_L1[kMaxGEM], GE11_z_L2[kMaxGEM], GE21_z_L1[kMaxGEM], GE21_z_L2[kMaxGEM];

  Float_t GE0_phi[kMaxGEM];
  Float_t GE0_phib[kMaxGEM];

  Float_t GE0_sim_phi[kMaxGEM];
  Float_t GE0_sim_phib[kMaxGEM];

  // positions from artificial pads
  Float_t GE21_pad1_phi_L1[kMaxGEM], GE21_pad1_phi_L2[kMaxGEM];
  Float_t GE21_pad2_phi_L1[kMaxGEM], GE21_pad2_phi_L2[kMaxGEM];
  Float_t GE21_pad4_phi_L1[kMaxGEM], GE21_pad4_phi_L2[kMaxGEM];
  Float_t GE21_pad8_phi_L1[kMaxGEM], GE21_pad8_phi_L2[kMaxGEM];
};

double 
dxy(double px, double py, double vx, double vy, double pt)
{
  //Source: https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrackReco/interface/TrackBase.h#119
  return (- vx * py + vy * px ) / pt;
}

double 
phiHeavyCorr(double pt, double eta, double phi, double q)
{
  //  float resEta = eta;
  float etaProp = std::abs(eta);
  if (etaProp< 1.1) etaProp = 1.1;
  float resPhi = phi - 1.464*q*cosh(1.7)/cosh(etaProp)/pt - M_PI/144.;
  if (resPhi > M_PI) resPhi -= 2.*M_PI;
  if (resPhi < -M_PI) resPhi += 2.*M_PI;
  return resPhi;
}

double 
dRWeighted(double eta1, double phi1, double eta2, double phi2, double sigma_eta=2., double sigma_phi=1.)
{
  double dEta = std::abs(eta1 - eta2);
  double dPhi = reco::deltaPhi(phi1, phi2);
  double dR = std::sqrt((dEta*dEta)/(sigma_eta*sigma_eta) + (dPhi*dPhi)/(sigma_phi*sigma_phi));
  return dR;
}

double 
My_dPhi(double phi1, double phi2) {
  double dPhi = phi1 - phi2;
  if (dPhi >  M_PI) dPhi -= 2.*M_PI;
  if (dPhi < -M_PI) dPhi += 2.*M_PI;
  return dPhi;
}

double 
phiL1DTTrack(const L1MuDTTrack& track)
{
  int phi_local = track.phi_packed(); //range: 0 < phi_local < 31
  if ( phi_local > 15 ) phi_local -= 32; //range: -16 < phi_local < 15    
  double dttrk_phi_global = normalizedPhi((phi_local*(M_PI/72.))+((M_PI/6.)*track.spid().sector()));// + 12*i->scNum(); //range: -16 < phi_global < 147 
  // if(dttrk_phi_global < 0) dttrk_phi_global+=2*M_PI; //range: 0 < phi_global < 147
  // if(dttrk_phi_global > 2*M_PI) dttrk_phi_global-=2*M_PI; //range: 0 < phi_global < 143
  return dttrk_phi_global;
}

double 
phiL1CSCTrack(const csc::L1Track& track)
{
  unsigned gbl_phi(track.localPhi() + ((track.sector() - 1)*24) + 6);
  if(gbl_phi > 143) gbl_phi -= 143;
  double phi_packed = gbl_phi & 0xff;
  return phi_packed;
}

int
getHalfStrip(const CSCComparatorDigi& digi) 
{
  return (digi.getStrip() - 1) * 2 + digi.getComparator();
}

float 
getFractionalStrip(const CSCComparatorDigi&d)
{
  return d.getStrip() + d.getComparator()/2. - 3/4.;
}

void fitStraightLineErrors(const std::vector<float>& v, 
                           const std::vector<float>& w, 
                           const std::vector<float>& ev, 
                           const std::vector<float>& ew, 
                           float& alpha, float& beta, 
                           int lumi, int run, int event, int muon, int stub,
                           bool debug)
{
  //std::cout << "size of v: "<<v.size() << std::endl; 
  
  if (v.size()>=3) {
  
  float zmin;
  float zmax;
  if (v.front() < v.back()){
    zmin = v.front();
    zmax = v.back();
  }
  else{
    zmin = v.back();
    zmax = v.front();
  }

  TF1 *fit1 = new TF1("fit1","pol1",zmin,zmax); 
  //where 0 = x-axis_lowest and 48 = x_axis_highest 
  TGraphErrors* gr = new TGraphErrors(v.size(),&(v[0]),&(w[0]),&(ev[0]),&(ew[0]));
  gr->SetMinimum(w[2]-5*0.002);
  gr->SetMaximum(w[2]+5*0.002);
 
  gr->Fit(fit1,"RQ"); 
  
  alpha = fit1->GetParameter(0); //value of 0th parameter
  beta  = fit1->GetParameter(1); //value of 1st parameter

  if (debug){
    TCanvas* c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
    c1->cd();
    c1->SetFillColor(42);
    c1->SetGrid();
    // c1->GetFrame()->SetFillColor(21);
    // c1->GetFrame()->SetBorderSize(12);
    
    TString slumi;  slumi.Form("%d", lumi);
    TString srun;   srun.Form("%d", run);
    TString sevent; sevent.Form("%d", event);
    TString smuon;  smuon.Form("%d", muon);
    TString sstub;  sstub.Form("%d", stub);

    gr->SetTitle("Linear fit to ComparatorDigis for Lumi " + slumi + " Run " + srun + " Event " + sevent + " Muon " + smuon + " Stub " + sstub);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
    
    c1->SaveAs("ComparatoDigiLinearFits/c_debug_fit_L" + slumi + "_R" + srun + "_E" + sevent + "_M" + smuon + "_S" + sstub + ".png");
    delete c1;
  }
  delete fit1;
  delete gr;
  }
  else{
    //std::cout << "ERROR: LCT without at least 3 comparator digis in the chamber!!" << std::endl;
  }
}


void fitStraightLine(const std::vector<float>& v, 
                     const std::vector<float>& w, 
                     float& alpha, float& beta, 
                     bool debug = false)
{
  if (v.size()>=2) {
    
    float zmin;
    float zmax;
    if (v.front() < v.back()){
      zmin = v.front();
      zmax = v.back();
    }
    else{
      zmin = v.back();
      zmax = v.front();
    }
    
    TF1 *fit1 = new TF1("fit1","pol1",zmin,zmax); 
    //where 0 = x-axis_lowest and 48 = x_axis_highest 
    TGraph* gr = new TGraph(v.size(),&(v[0]),&(w[0]));
    gr->Fit(fit1,"RQ"); 
    
    alpha = fit1->GetParameter(0); //value of 0th parameter
    beta  = fit1->GetParameter(1); //value of 1st parameter
    
    delete fit1;
    delete gr;
  }
  else{
    alpha = 0;
    beta = 0;
  }
}


void getPositionsStations(float alpha_x, float beta_x, float alpha_y, float beta_y, 
                          std::vector<float>& xs, std::vector<float>& ys, 
                          int sign_z,
                          float z1 = 600, float z2 = 825, float z3 = 935, float z4 = 1020)
{
  xs.push_back(alpha_x + beta_x * z1*sign_z);
  xs.push_back(alpha_x + beta_x * z2*sign_z);
  xs.push_back(alpha_x + beta_x * z3*sign_z);
  xs.push_back(alpha_x + beta_x * z4*sign_z);
  
  ys.push_back(alpha_y + beta_y * z1*sign_z);
  ys.push_back(alpha_y + beta_y * z2*sign_z);
  ys.push_back(alpha_y + beta_y * z3*sign_z);
  ys.push_back(alpha_y + beta_y * z4*sign_z);
}


bool 
isSimTrackGood(const SimTrack &t)
{
  // SimTrack selection
  //if (t.noVertex()) return false;
  //if (t.noGenpart()) return false;
  // only muons 
  if (std::abs(t.type()) != 13) return false;
  // pt selection
  if (t.momentum().pt() < 0) return false;
  // eta selection
  const float eta(std::abs(t.momentum().eta()));
  if (eta > 3.0) return false; 
  return true;
}

using namespace std;

class DisplacedL1MuFilter : public edm::EDFilter 
{
public:
  explicit DisplacedL1MuFilter(const edm::ParameterSet&);
  ~DisplacedL1MuFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void bookL1MuTree();
  void clearBranches();

  virtual void beginJob() override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  float getGlobalPhi(unsigned int rawid, int stripN);
  double calcCSCSpecificPhi(unsigned int rawId, const CSCCorrelatedLCTDigi& tp) const;
  GlobalPoint getGlobalPointPad(unsigned int rawId, const GEMPadDigi& tp) const;;  
  GlobalPoint getCSCSpecificPoint(unsigned int rawid, const CSCCorrelatedLCTDigi& tp) const;
  GlobalPoint getCSCSpecificPoint2(unsigned int rawId, const CSCCorrelatedLCTDigi& tp) const;
  
  bool isCSCCounterClockwise(const std::unique_ptr<const CSCLayer>& layer) const;
  void fitComparatorsLCT(const CSCComparatorDigiCollection&, const CSCCorrelatedLCTDigi& tp, 
                         CSCDetId chid, int iMuon, float& fit_z, float& fit_phi, float& dphi) const;
  void getStubPositions(int index, std::vector<float>& x, 
                        std::vector<float>& y, std::vector<float>& z) const;

  void printCSCStubProperties(int index) const;

  CSCCorrelatedLCTDigiId 
  pickBestMatchingStub(float x1, float y1,
                       const CSCCorrelatedLCTDigiId& oldCoPad,
                       const CSCCorrelatedLCTDigiId& newCoPad, 
                       int refBx) const;
  
  bool stubInCSCTFTracks(const CSCCorrelatedLCTDigi& stub, const L1CSCTrackCollection& l1Tracks) const;

  GEMPadDigiId
  pickBestMatchingCoPad(float x1, float y1,
                        const GEMPadDigiId& oldCoPad,
                        const GEMPadDigiId& newCoPad, int refBx) const;
  GEMPadDigiId
  pickBestMatchingPad(float x1, float y1,
                      const GEMPadDigiId& oldPad,
                      const GEMPadDigiId& newPad, int refBx) const;
  void fillCSCStubProperties(const CSCDetId& ch_id,
                             const CSCCorrelatedLCTDigi& stub,
                             int index,
                             const GlobalPoint& GP,
                             float z, float phi, float dphi);

  // ----------member data ---------------------------

  enum WhichTrack { None, TrackerTk, MuonTk, GlobalTk };
  enum WhichState { AtVertex, Innermost, Outermost };

  /// Use simplified geometry (cylinders and disks, not individual chambers)
  bool useSimpleGeometry_;

  /// Propagate to MB2 (default) instead of MB1
  bool useMB2_;

  /// Fallback to ME1 if propagation to ME2 fails
  bool fallbackToME1_;

  /// Labels for input collections
  WhichTrack whichTrack_;
  WhichState whichState_;

  /// for cosmics, some things change: the along-opposite is not in-out, nor the innermost/outermost states are in-out really
  bool cosmicPropagation_;

  int min_L1Mu_Quality;
  double max_dR_L1Mu_L1Tk;
  double max_dR_L1Mu_noL1Tk;
  double min_pT_L1Tk;
  double max_pT_L1Tk;
  int verbose;
  bool produceFitPlots_;
  bool processRPCb_;
  bool processRPCf_;
  bool doSimAnalysis_;
  
  const RPCGeometry* rpcGeometry_;
  const CSCGeometry* cscGeometry_;
  const GEMGeometry* gemGeometry_;
 
  edm::InputTag L1Mu_input;
  edm::InputTag L1TkMu_input;
  
  edm::ESHandle<MagneticField> magfield_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<Propagator> propagatorOpposite_;
  edm::ESHandle<Propagator> propagatorAny_;
  edm::ESHandle<MuonDetLayerGeometry> muonGeometry_;
  edm::ESHandle<RPCGeometry> rpc_geom_;
  edm::ESHandle<CSCGeometry> csc_geom_;
  edm::ESHandle<GEMGeometry> gem_geom_;
  
  const  BoundCylinder *barrelCylinder_;
  const  BoundDisk *endcapDiskPos_[3], *endcapDiskNeg_[3];
  double barrelHalfLength_;
  std::pair<float,float> endcapRadii_[3];
  MyEvent event_;
  TTree* event_tree_;

  // trigger scale
  unsigned long long  muScalesCacheID_;
  unsigned long long  muPtScaleCacheID_;

  edm::ESHandle< L1MuTriggerScales > muScales;
  edm::ESHandle< L1MuTriggerPtScale > muPtScale;

  edm::ParameterSet cfg_;
};

DisplacedL1MuFilter::DisplacedL1MuFilter(const edm::ParameterSet& iConfig) : 
  useSimpleGeometry_(iConfig.getParameter<bool>("useSimpleGeometry")),
  useMB2_(iConfig.existsAs<bool>("useStation2") ? iConfig.getParameter<bool>("useStation2") : true),
  fallbackToME1_(iConfig.existsAs<bool>("fallbackToME1") ? iConfig.getParameter<bool>("fallbackToME1") : false),
  whichTrack_(None), whichState_(AtVertex),
  cosmicPropagation_(iConfig.getParameter<bool>("cosmicPropagationHypothesis")),
  cfg_(iConfig.getParameterSet("simTrackMatching"))
{
  //now do what ever initialization is needed
  min_L1Mu_Quality = iConfig.getParameter<int>("min_L1Mu_Quality");
  max_dR_L1Mu_L1Tk = iConfig.getParameter<double>("max_dR_L1Mu_L1Tk");
  max_dR_L1Mu_noL1Tk = iConfig.getParameter<double>("max_dR_L1Mu_noL1Tk");
  min_pT_L1Tk = iConfig.getParameter<double>("min_pT_L1Tk");
  max_pT_L1Tk = iConfig.getParameter<double>("max_pT_L1Tk");
  verbose = iConfig.getParameter<int>("verbose");
  produceFitPlots_ = iConfig.getParameter<bool>("produceFitPlots");
  processRPCb_ = iConfig.getParameter<bool>("processRPCb");
  processRPCf_ = iConfig.getParameter<bool>("processRPCf");
  doSimAnalysis_ = iConfig.getParameter<bool>("doSimAnalysis");
  
  L1Mu_input = iConfig.getParameter<edm::InputTag>("L1Mu_input");
  L1TkMu_input = iConfig.getParameter<edm::InputTag>("L1TkMu_input");
  
  bookL1MuTree();

  muScalesCacheID_ = 0ULL ;
  muPtScaleCacheID_ = 0ULL ;

}

DisplacedL1MuFilter::~DisplacedL1MuFilter()
{
}


bool
DisplacedL1MuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  clearBranches();
  
  // propagator
  iSetup.get<IdealMagneticFieldRecord>().get(magfield_);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);
  iSetup.get<MuonRecoGeometryRecord>().get(muonGeometry_);

  // Get the barrel cylinder
  const DetLayer * dtLay = muonGeometry_->allDTLayers()[useMB2_ ? 1 : 0];
  barrelCylinder_ = dynamic_cast<const BoundCylinder *>(&dtLay->surface());
  barrelHalfLength_ = barrelCylinder_->bounds().length()/2;;
  // std::cout << "L1MuonMatcher: barrel radius = " << barrelCylinder_->radius() << ", half length = " << barrelHalfLength_ << std::endl;

  // Get the endcap disks. Note that ME1 has two disks (ME1/1 and ME2/1-ME3/2-ME4/1), so there's one more index
  for (size_t i = 0; i <= (useMB2_ ? 2 : 1); ++i) {
    endcapDiskPos_[i] = dynamic_cast<const BoundDisk *>(& muonGeometry_->forwardCSCLayers()[i]->surface());
    endcapDiskNeg_[i] = dynamic_cast<const BoundDisk *>(& muonGeometry_->backwardCSCLayers()[i]->surface());
    endcapRadii_[i] = std::make_pair(endcapDiskPos_[i]->innerRadius(), endcapDiskPos_[i]->outerRadius());
    // std::cout << "L1MuonMatcher: endcap " << i << " Z = " << endcapDiskPos_[i]->position().z() << ", radii = " << endcapRadii_[i].first << "," << endcapRadii_[i].second << std::endl;
  }

  iSetup.get<MuonGeometryRecord>().get(rpc_geom_);
  rpcGeometry_ = &*rpc_geom_;

  iSetup.get<MuonGeometryRecord>().get(csc_geom_);
  cscGeometry_ = &*csc_geom_;

  iSetup.get<MuonGeometryRecord>().get(gem_geom_);
  gemGeometry_ = &*gem_geom_;

  typedef std::vector<L1MuGMTCand> GMTs;
  edm::Handle<GMTs> aH;
  iEvent.getByLabel("simGmtDigis", aH);
  const GMTs& l1GmtCands(*aH.product());

  edm::Handle<vector<pair<L1MuDTTrack,vector<L1MuDTTrackSegPhi> > > > L1DTTrackPhiH;
  iEvent.getByLabel("dttfDigis","DTTF", L1DTTrackPhiH);
  const vector<pair<L1MuDTTrack,vector<L1MuDTTrackSegPhi> > >& L1DTTrackPhis(*L1DTTrackPhiH.product());

  // tracks produced by TF
  edm::Handle< L1CSCTrackCollection > hl1Tracks;
  iEvent.getByLabel("simCsctfTrackDigis", hl1Tracks);
  const L1CSCTrackCollection& l1Tracks(*hl1Tracks.product());

  if (iSetup.get< L1MuTriggerScalesRcd >().cacheIdentifier() != muScalesCacheID_ ||
      iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier() != muPtScaleCacheID_ )
    {
      iSetup.get< L1MuTriggerScalesRcd >().get( muScales );
      iSetup.get< L1MuTriggerPtScaleRcd >().get( muPtScale );

      muScalesCacheID_  = iSetup.get< L1MuTriggerScalesRcd >().cacheIdentifier();
      muPtScaleCacheID_ = iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier();
    }

  edm::Handle< std::vector<L1MuRegionalCand> > hL1MuRPCbs;
  iEvent.getByLabel("simRpcTriggerDigis", "RPCb", hL1MuRPCbs);
  const std::vector<L1MuRegionalCand>& l1MuRPCbs(*hL1MuRPCbs.product());
  
  edm::Handle< std::vector<RPCDigiL1Link> > hL1MuRPCbLinks;
  iEvent.getByLabel("simRpcTriggerDigis", "RPCb", hL1MuRPCbLinks);
  const std::vector<RPCDigiL1Link>& l1MuRPCbLinks(*hL1MuRPCbLinks.product());
  
  edm::Handle< std::vector<L1MuRegionalCand> > hL1MuRPCfs;
  iEvent.getByLabel("simRpcTriggerDigis", "RPCf", hL1MuRPCfs);
  const std::vector<L1MuRegionalCand>& l1MuRPCfs(*hL1MuRPCfs.product());
  
  edm::Handle< std::vector<RPCDigiL1Link> > hL1MuRPCfLinks;
  iEvent.getByLabel("simRpcTriggerDigis", "RPCf", hL1MuRPCfLinks);
  const std::vector<RPCDigiL1Link>& l1MuRPCfLinks(*hL1MuRPCfLinks.product());

  // comparator digis
  edm::Handle< CSCComparatorDigiCollection > hCSCComparators;
  iEvent.getByLabel("simMuonCSCDigis", "MuonCSCComparatorDigi", hCSCComparators);
  //const CSCComparatorDigiCollection& CSCComparators(*hCSCComparators.product());
  
  // LCTs
  edm::Handle< CSCCorrelatedLCTDigiCollection > hCSCCorrelatedLCTs;
  iEvent.getByLabel("simCscTriggerPrimitiveDigis", "MPCSORTED", hCSCCorrelatedLCTs);
  const CSCCorrelatedLCTDigiCollection& CSCCorrelatedLCTs(*hCSCCorrelatedLCTs.product());

  // GEM pads and copads
  edm::Handle< GEMPadDigiCollection > hGEMPads;
  iEvent.getByLabel("simMuonGEMPadDigis", hGEMPads);
  //const GEMPadDigiCollection& GEMPads(*hGEMPads.product());

  edm::Handle< GEMPadDigiCollection > hGEMCoPads;
  iEvent.getByLabel("simMuonGEMPadDigis", "Coincidence", hGEMCoPads);
  //const GEMCoPadDigiCollection& GEMCoPads(*hGEMCoPads.product());

  // L1 TrackingTrigger Analysis
  // edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  // iEvent.getByLabel("TTTracksFromPixelDigis", "Level1TTTracks", TTTrackHandle);
  // const std::vector< TTTrack< Ref_PixelDigi_ > >& TTTracks = *TTTrackHandle.product();

  event_.lumi = iEvent.id().luminosityBlock();
  event_.run = iEvent.id().run();
  event_.event = iEvent.id().event();

  event_.nL1Mu = l1GmtCands.size();
  // event_.nL1Tk = TTTracks.size();
  
  event_.beamSpot_x = 0;
  event_.beamSpot_y = 0;
  event_.beamSpot_z = 0;

  /////////////////
  // L1 analysis //
  /////////////////

  // // store the L1Tk variables
  // for (unsigned int j=0; j<TTTracks.size(); ++j) {
  //   auto l1Tk = TTTracks[j];
  //   const double l1Tk_pt = l1Tk.getMomentum().perp();
  //   const double l1Tk_eta = l1Tk.getMomentum().eta();
  //   const double l1Tk_phi = normalizedPhi(l1Tk.getMomentum().phi());
  //   const double l1Tk_charge = l1Tk.getRInv()>0? 1: -1;
  //   const double l1Tk_eta_corr = l1Tk_eta;
  //   const double l1Tk_phi_corr = phiHeavyCorr(l1Tk_pt, l1Tk_eta, l1Tk_phi, l1Tk_charge);

  //   if(verbose and false) {
  //     cout << "l1Tk " << j << endl; 
  //     cout << "l1Tk_pt " << l1Tk_pt << endl;
  //     cout << "l1Tk_eta " << l1Tk_eta << endl;
  //     cout << "l1Tk_phi " << l1Tk_phi << endl;
  //     cout << "l1Tk_phi_corr " << l1Tk_phi_corr << endl;
  //     cout << "l1Tk_charge " << l1Tk_charge << endl;
  //   }
      
  //   double l1Tk_eta_prop = -99;
  //   double l1Tk_phi_prop = -99;
  //   GlobalPoint ex_point(extrapolateGP(l1Tk));
  //   if (!(ex_point == GlobalPoint())) {
  //     l1Tk_eta_prop = ex_point.eta();
  //     l1Tk_phi_prop = ex_point.phi();
  //     if(verbose and false) {
  //       cout << "l1Tk_eta_prop " << l1Tk_eta_prop << endl;
  //       cout << "l1Tk_phi_prop " << l1Tk_phi_prop << endl;
  //     }
  //   }
    
  //   event_.L1Tk_pt[j] = l1Tk_pt;
  //   event_.L1Tk_eta[j] = l1Tk_eta;
  //   event_.L1Tk_phi[j] = l1Tk_phi;
    
  //   event_.L1Tk_eta_prop[j] = l1Tk_eta_prop;
  //   event_.L1Tk_phi_prop[j] = l1Tk_phi_prop;
  //   event_.L1Tk_deta_prop[j] = std::abs(event_.L1Tk_eta_prop[j] - event_.L1Tk_eta[j]);
  //   event_.L1Tk_dphi_prop[j] = My_dPhi(event_.L1Tk_phi_prop[j], event_.L1Tk_phi[j]); 
  //   event_.L1Tk_dR_prop[j] = reco::deltaR(event_.L1Tk_eta_prop[j], event_.L1Tk_phi_prop[j], event_.L1Tk_eta[j], event_.L1Tk_phi[j]);
    
  //   event_.L1Tk_eta_corr[j] = l1Tk_eta_corr;
  //   event_.L1Tk_phi_corr[j] = l1Tk_phi_corr;
  //   event_.L1Tk_deta_corr[j] = std::abs(event_.L1Tk_eta_corr[j] - event_.L1Tk_eta[j]);
  //   event_.L1Tk_dphi_corr[j] = My_dPhi(event_.L1Tk_phi_corr[j], event_.L1Tk_phi[j]); 
  //   event_.L1Tk_dR_corr[j] = reco::deltaR(event_.L1Tk_eta_corr[j], event_.L1Tk_phi_corr[j], event_.L1Tk_eta[j], event_.L1Tk_phi[j]);
  // } // end of loop on TTTracks
  
  

  // store the DTTF variables
  if(verbose) std::cout << std::endl<< "Number of L1DTTrackPhis " <<L1DTTrackPhis.size() << std::endl;
  event_.nDTTF = L1DTTrackPhis.size();
  for (unsigned int j=0; j<L1DTTrackPhis.size(); ++j) { 
    auto track = L1DTTrackPhis[j].first;
    
    event_.DTTF_pt[j] = muPtScale->getPtScale()->getLowEdge(track.pt()) + 1.e-6;
    event_.DTTF_eta[j] = muScales->getRegionalEtaScale(0)->getCenter(track.eta());
    event_.DTTF_phi[j] = phiL1DTTrack(track);
    event_.DTTF_bx[j] = track.bx();
    event_.DTTF_nStubs[j] = L1DTTrackPhis[j].second.size();
    event_.DTTF_quality[j] = track.quality();
    
    if(verbose) {  
      std::cout << std::endl
                << "pt  = " << event_.DTTF_pt[j]
                << ", eta  = " << event_.DTTF_eta[j] 
                << ", phi  = " << event_.DTTF_phi[j]
                << ", bx = " << event_.DTTF_bx[j] 
                << ", quality = " << event_.DTTF_quality[j] 
                << ", nStubs = " << event_.DTTF_nStubs[j]
                << std::endl;
    }
    
    for (auto stub: L1DTTrackPhis[j].second) {
      // std::cout << "\t " << stub << std::endl;
      // std::cout << "\t phiValue = " << stub.phiValue() << ", phibValue = " << stub.phibValue() << std::endl;
      int station = stub.station();
      switch(station) {
      case 1:
        event_.DTTF_phi1[j] = stub.phiValue();
        event_.DTTF_phib1[j] = stub.phibValue();
        event_.DTTF_quality1[j] = stub.quality();
        event_.DTTF_bx1[j] = stub.bx();
        event_.DTTF_wh1[j] = stub.wheel();
        event_.DTTF_se1[j] = stub.sector();
        event_.DTTF_st1[j] = stub.station();
        break;
      case 2:
        event_.DTTF_phi2[j] = stub.phiValue();
        event_.DTTF_phib2[j] = stub.phibValue();
        event_.DTTF_quality2[j] = stub.quality();
        event_.DTTF_bx2[j] = stub.bx();
        event_.DTTF_wh2[j] = stub.wheel();
        event_.DTTF_se2[j] = stub.sector();
        event_.DTTF_st2[j] = stub.station();
        break;
      case 3:
        event_.DTTF_phi3[j] = stub.phiValue();
        event_.DTTF_phib3[j] = stub.phibValue();
        event_.DTTF_quality3[j] = stub.quality();
        event_.DTTF_bx3[j] = stub.bx();
        event_.DTTF_wh3[j] = stub.wheel();
        event_.DTTF_se3[j] = stub.sector();
        event_.DTTF_st3[j] = stub.station();
        break;
      case 4:
        event_.DTTF_phi4[j] = stub.phiValue();
        event_.DTTF_phib4[j] = stub.phibValue();
        event_.DTTF_quality4[j] = stub.quality();
        event_.DTTF_bx4[j] = stub.bx();
        event_.DTTF_wh4[j] = stub.wheel();
        event_.DTTF_se4[j] = stub.sector();
        event_.DTTF_st4[j] = stub.station();
        break;
      };
    }
  }


  // store the CSCTF variables
  event_.nCSCTF = l1Tracks.size();
  if(verbose) std::cout << std::endl<< "Number of L1CSCTracks " <<event_.nCSCTF << std::endl;
  for (int j=0; j<event_.nCSCTF; ++j) {
    auto track = l1Tracks[j].first;
    const int sign(track.endcap()==1 ? 1 : -1);
    unsigned gpt = 0, quality = 0;
    csc::L1Track::decodeRank(track.rank(), gpt, quality);

    // calculate pt, eta and phi (don't forget to store the sign)                                                                                   
    event_.CSCTF_pt[j] = muPtScale->getPtScale()->getLowEdge(gpt & 0x1f) + 1.e-6;
    event_.CSCTF_eta[j] = muScales->getRegionalEtaScale(2)->getCenter(track.eta_packed()) * sign;
    event_.CSCTF_phi[j] = normalizedPhi(muScales->getPhiScale()->getLowEdge(phiL1CSCTrack(track)));
    event_.CSCTF_bx[j] = track.bx();
    event_.CSCTF_quality[j] = quality;
    
    if(verbose) {  
      std::cout << std::endl
                << "L1CSC pt  = " << event_.CSCTF_pt[j]
                << ", eta  = " << event_.CSCTF_eta[j] 
                << ", phi  = " << event_.CSCTF_phi[j]
                << ", bx = " << event_.CSCTF_bx[j]
                << ", quality = " << event_.CSCTF_quality[j] 
                << std::endl<< std::endl;
      std::cout << "Print CSCTF stubs:" << std::endl;
    }
    event_.CSCTF_nStubs[j] = 0;
    auto stubCollection(l1Tracks[j].second);
    for (auto detUnitIt = stubCollection.begin(); detUnitIt != stubCollection.end(); detUnitIt++) {
      const CSCDetId& ch_id = (*detUnitIt).first;
      const auto range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        //if (!(*digiIt).isValid()) continue;
        event_.CSCTF_nStubs[j] += 1;
        auto stub = *digiIt;
        auto gp = getCSCSpecificPoint2(ch_id.rawId(), stub);
        if(verbose)std::cout << "\t" << ch_id << " " << stub 
                             << " key eta " << gp.eta() << " key phi " << gp.phi() << std::endl;
        
        // stub position
        
        float z_pos_L3, bestFitPhi, bestFitDPhi;

        fitComparatorsLCT(*hCSCComparators.product(), stub, ch_id, int(digiIt-range.first), z_pos_L3, bestFitPhi, bestFitDPhi);

        fillCSCStubProperties(ch_id, stub, j, gp, z_pos_L3, bestFitPhi, bestFitDPhi);
      }
    }
    
    /* 
    CSCTF stub recovery per CSCTF track
    The CSC track-finder may drop certain stubs if they don't match the pattern
    First get the station numbers where stubs are not filled
    */
    bool stubMissingSt1 = event_.CSCTF_st1[j] == 99;
    bool stubMissingSt2 = event_.CSCTF_st2[j] == 99;
    bool stubMissingSt3 = event_.CSCTF_st3[j] == 99;
    bool stubMissingSt4 = event_.CSCTF_st4[j] == 99;
    bool atLeast1StubMissing = stubMissingSt1 or stubMissingSt2 or stubMissingSt3 or stubMissingSt4;
    
    bool doStubRecovery = true;
    if (doStubRecovery and atLeast1StubMissing){
      
      std::vector<float> xs;
      std::vector<float> ys;
      std::vector<float> zs;
      //std::cout << "Get stub positions" << std::endl;
      getStubPositions(j, xs, ys, zs);
      
      //std::cout << "fit stub positions with straight line" << std::endl;
      float alpha_x, beta_x, alpha_y, beta_y;
      fitStraightLine(zs, xs, alpha_x, beta_x); 
      fitStraightLine(zs, ys, alpha_y, beta_y); 
      
      //std::cout << "Get positions stations" << std::endl;
      std::vector<float> allxs;
      std::vector<float> allys;
      int sign_z = int(event_.CSCTF_eta[j]/std::abs(event_.CSCTF_eta[j]));
      getPositionsStations(alpha_x, beta_x, alpha_y, beta_y,
                           allxs, allys, sign_z);

      
      int triggerSector = track.sector();
      //std::cout << "trigger sector " << triggerSector << std::endl;
      
      for (int endcap=1; endcap<=2; endcap++){
        //do not consider stubs in the wrong endcap
        int zendcap(endcap!=1 ? -1 : +1 );
        if (zendcap * event_.CSCTF_eta[j] < 0) continue;
        for (int station=1; station<=4; station++){
          
          // ignore station where a L1Mu stub is present!
          if (not stubMissingSt1 and station==1) continue;
          if (not stubMissingSt2 and station==2) continue;
          if (not stubMissingSt3 and station==3) continue;
          if (not stubMissingSt4 and station==4) continue;
          if(verbose) std::cout << "Recovered stubs in station: " << station << std::endl;
          // temp storage of candidate stubs per station and ring
          CSCCorrelatedLCTDigiId bestMatchingStub;
          int iStub = 0;
          for (int ring=1; ring<=3; ring++){
            if (station!=1 and ring==3) continue;
            //std::cout << "Analyzing ring " << ring << std::endl;
            
            for (int chamber=1; chamber<=36; chamber++){
              // do not consider invalid detids
              if ( (station==2 or station==3 or station==4) and 
                   (ring==1) and chamber>18) continue;
              //std::cout << "Analyzing chamber " << chamber << std::endl; 
              // create the detid
              CSCDetId ch_id(endcap, station, ring, chamber);
              //std::cout << "ch_id " <<  ch_id << std::endl;
              // get the stubs in this detid
              auto range = CSCCorrelatedLCTs.get(ch_id);
              for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
                iStub++; 

                auto stub(*digiItr);

                // check that this stub is not already part of a CSC TF track
                if (stubInCSCTFTracks(stub, l1Tracks)) continue;

                // trigger sector must be the same
                if (triggerSector != ch_id.triggerSector()) continue;
                
                // BXs have to match
                int deltaBX = std::abs(stub.getBX() - (6 + event_.CSCTF_bx[j]));
                if (deltaBX > 1) continue;

                std::cout << ch_id << std::endl;
                std::cout<<"Candidate " << stub << std::endl;
                bestMatchingStub = pickBestMatchingStub(allxs[ch_id.station()-1], allys[ch_id.station()-1], 
                                                        bestMatchingStub, std::make_pair(ch_id, stub), 6 + event_.CSCTF_bx[j]);
              }
              // consider the case ME1a
              if (station==1 and ring==1){
                CSCDetId me1a_id(endcap, station, 4, chamber);
                auto range = CSCCorrelatedLCTs.get(me1a_id);
                for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
                  iStub++;
                  auto stub(*digiItr);

                  // check that this stub is not already part of a CSC TF track
                  if (stubInCSCTFTracks(stub, l1Tracks)) continue;
                  
                  // trigger sector must be the same
                  if (triggerSector != me1a_id.triggerSector()) continue;
                  
                  // BXs have to match
                  int deltaBX = std::abs(stub.getBX() - (6 + event_.CSCTF_bx[j]));
                  if (deltaBX > 1) continue; 

                  std::cout << me1a_id << std::endl;
                  std::cout<<"Candidate " << stub << std::endl;
                  bestMatchingStub = pickBestMatchingStub(allxs[me1a_id.station()-1], allys[me1a_id.station()-1], 
                                                          bestMatchingStub, std::make_pair(me1a_id, stub), 6 + event_.CSCTF_bx[j]);
                  
                }
              }
            }
          }
          if (bestMatchingStub.second != CSCCorrelatedLCTDigi()) {
            // stub position
            auto gp = getCSCSpecificPoint2(bestMatchingStub.first.rawId(), bestMatchingStub.second);
            if (reco::deltaR(gp.eta(), normalizedPhi(gp.phi()), 
                             (float)event_.CSCTF_eta[j] , normalizedPhi((float)event_.CSCTF_phi[j])) < 0.2){ 

              if(verbose) {
                std::cout << "\tChoice:" 
                          << bestMatchingStub.first << " " 
                          << bestMatchingStub.second << " " 
                          << "key eta " << gp.eta() << " "
                          << "key phi " << gp.phi() << std::endl;
              }
              
              float z_pos_L3, bestFitPhi, bestFitDPhi;
              fitComparatorsLCT(*hCSCComparators.product(), bestMatchingStub.second, bestMatchingStub.first, j, z_pos_L3, bestFitPhi, bestFitDPhi);
              
              fillCSCStubProperties(bestMatchingStub.first, bestMatchingStub.second, j, 
                                    gp, z_pos_L3, bestFitPhi, bestFitDPhi);
            }
            else{
              if(verbose) if (iStub!=0) std::cout << "\tNone " << std::endl;
            }
          }
          else{
            if(verbose) if (iStub!=0) std::cout << "\tNone " << std::endl;
          }
        }
      } 
    }
    stubMissingSt1 = event_.CSCTF_st1[j] == 99;
    stubMissingSt2 = event_.CSCTF_st2[j] == 99;
    cout << "Stub missing in station 1: " << bool(stubMissingSt1) << endl; 
    cout << "Stub missing in station 2: " << bool(stubMissingSt2) << endl; 
    /*

    std::cout << "Matching copads" << std::endl;
    // Get matching GEM copads

    GEMPadDigiId bestCoPad_GE11;
    GEMPadDigiId bestCoPad_GE21;

    for(auto cItr = hGEMCoPads->begin(); cItr != hGEMCoPads->end(); ++cItr) {
      GEMDetId gem_id = (*cItr).first;
      if (not stubMissingSt1) {
        // get the CSCDetId of station 1
        CSCDetId csc_st1(event_.CSCTF_id1[ j ]);
        std::cout << "CSC " << csc_st1 << endl;
        // chambers need to be compatible
        if (gem_id.station() != 1 or 
            csc_st1.chamber() != gem_id.chamber() or
            (csc_st1.ring() != 4 and csc_st1.ring() != 1)) continue;
        std::cout << "Investigate GE11 chamber " << gem_id << std::endl;
        // get the copads
        auto copad_range = (*cItr).second;
        for (auto digiItr = copad_range.first; digiItr != copad_range.second; ++digiItr){
          auto copad(*digiItr);
          std::cout << "\tCandidate copad " << copad  << std::endl;
          bestCoPad_GE11 = pickBestMatchingCoPad(event_.CSCTF_fit_x1[ j ],
                                                 event_.CSCTF_fit_y1[ j ], 
                                                 bestCoPad_GE11, std::make_pair(gem_id, copad), event_.CSCTF_bx[j]);
        }
      }
      if (not stubMissingSt2) {
        // get the CSCDetId of station 1
        CSCDetId csc_st2(event_.CSCTF_id2[ j ]);
        std::cout << "CSC " << csc_st2 << endl;
        // chambers need to be compatible
        if (gem_id.station() != 3 or 
            csc_st2.chamber() != gem_id.chamber() or
            csc_st2.ring() != 1) continue;
        std::cout << "Investigate GE21 chamber " << gem_id << std::endl;
        // get the copads
        auto copad_range = (*cItr).second;
        for (auto digiItr = copad_range.first; digiItr != copad_range.second; ++digiItr){
          auto copad(*digiItr);
          std::cout << "\tCandidate copad " << copad  << std::endl;
          bestCoPad_GE21 = pickBestMatchingCoPad(event_.CSCTF_fit_x2[ j ],
                                                 event_.CSCTF_fit_y2[ j ], 
                                                 bestCoPad_GE21, std::make_pair(gem_id, copad), event_.CSCTF_bx[j]);
        }
      }
    }
    // found copad is not empty
    if (bestCoPad_GE11.second != GEMPadDigi()){
      std::cout << "\tBest GE11 copad" << bestCoPad_GE11.second << std::endl;
      if (bestCoPad_GE11.first.station()==1) {
        auto gem_gp1 = getGlobalPointPad(bestCoPad_GE11.first, bestCoPad_GE11.second);
        event_.GE11_phi_L1[j] = gem_gp1.phi();
        event_.GE11_bx_L1[j] = bestCoPad_GE11.second.bx();
        event_.GE11_ch_L1[j] = bestCoPad_GE11.first.chamber();
        event_.GE11_z_L1[j] = gem_gp1.z();
      }
    } else {
      std::cout << "No best GE11 copad" << std::endl;
    }
    
    if (bestCoPad_GE21.second != GEMPadDigi()){
      std::cout << "\tBest GE21 copad" << bestCoPad_GE21.second << std::endl;
      if (bestCoPad_GE21.first.station()==2) {
        auto gem_gp1 = getGlobalPointPad(bestCoPad_GE21.first, bestCoPad_GE21.second);
        event_.GE21_phi_L1[j] = gem_gp1.phi();
        event_.GE21_bx_L1[j] = bestCoPad_GE21.second.bx();
        event_.GE21_ch_L1[j] = bestCoPad_GE21.first.chamber();
        event_.GE21_z_L1[j] = gem_gp1.z();
      }
    } else {
      std::cout << "No best GE21 copad" << std::endl;
    } 
    
    

    // check copads were matched
    const bool GE11_copad_matched( event_.GE11_phi_L1[j]!=99 );
    const bool GE21_copad_matched( event_.GE21_phi_L1[j]!=99 );
    */

    std::cout << "Matching pads" << std::endl;
    // Get matching GEM pads
    if (true){
      //( (not GE11_copad_matched) or  (not GE21_copad_matched) )
      GEMPadDigiId bestPad_GE11_L1;
      GEMPadDigiId bestPad_GE11_L2;
      GEMPadDigiId bestPad_GE21_L1;
      GEMPadDigiId bestPad_GE21_L2;
      
      for(auto cItr = hGEMPads->begin(); cItr != hGEMPads->end(); ++cItr) {
        GEMDetId gem_id = (*cItr).first;
        //cout << "Check GEMDetId " << gem_id << endl;
        //float bestGEMDPhi = 99;
        if (not stubMissingSt1) {// and not GE11_copad_matched
          // get the CSCDetId of station 1
          CSCDetId csc_st1(event_.CSCTF_id1[ j ]);
          
          // chambers need to be compatible
          if (gem_id.region() == csc_st1.zendcap() and 
              gem_id.station() == 1 and 
              csc_st1.chamber() == gem_id.chamber() and
              (csc_st1.ring() == 4 or csc_st1.ring() == 1)) {
            std::cout << "Investigate GE11 chamber " << gem_id << std::endl;
            // get the pads
            auto pad_range = (*cItr).second;
            for (auto digiItr = pad_range.first; digiItr != pad_range.second; ++digiItr){
              auto GE11_pad(*digiItr);
              int deltaBX = std::abs(GE11_pad.bx() - event_.CSCTF_bx[j]);
              if (deltaBX <= 1) {
              
                std::cout << "\tCandidate Pad " << GE11_pad  << std::endl;
                if (gem_id.layer()==1){
                  // BXs have to match
                  
                  bestPad_GE11_L1 = pickBestMatchingPad(event_.CSCTF_fit_x1[ j ],
                                                        event_.CSCTF_fit_y1[ j ], 
                                                        bestPad_GE11_L1, std::make_pair(gem_id, GE11_pad), event_.CSCTF_bx[j]);
                }
                if (gem_id.layer()==2){
                  bestPad_GE11_L2 = pickBestMatchingPad(event_.CSCTF_fit_x1[ j ],
                                                        event_.CSCTF_fit_y1[ j ], 
                                                        bestPad_GE11_L2, std::make_pair(gem_id, GE11_pad), event_.CSCTF_bx[j]);
                }
              }
            }
          }
        }
        if (not stubMissingSt2) {// and not GE21_copad_matched
           // get the CSCDetId of station 1
          CSCDetId csc_st2(event_.CSCTF_id2[ j ]);
          //int index = i;
          // chambers need to be compatible
          if (gem_id.region() == csc_st2.zendcap() and
              gem_id.station() == 3 and 
              csc_st2.chamber() == gem_id.chamber() and
              csc_st2.ring() == 1) {
            std::cout << "Investigate GE21 chamber " << gem_id << std::endl;
            // get the pads
            auto pad_range = (*cItr).second;
            for (auto digiItr = pad_range.first; digiItr != pad_range.second; ++digiItr){
              auto pad(*digiItr);
              int deltaBX = std::abs(pad.bx() - event_.CSCTF_bx[j]);
              if (deltaBX <= 1) {
                std::cout << "\tCandidate Pad " << pad  << std::endl;
                if (gem_id.layer()==1){
                  bestPad_GE21_L1 = pickBestMatchingPad(event_.CSCTF_fit_x2[ j ],
                                                        event_.CSCTF_fit_y2[ j ], 
                                                        bestPad_GE21_L1, std::make_pair(gem_id, pad), event_.CSCTF_bx[j]);
                }
                if (gem_id.layer()==2){
                  bestPad_GE21_L2 = pickBestMatchingPad(event_.CSCTF_fit_x2[ j ],
                                                        event_.CSCTF_fit_y2[ j ], 
                                                        bestPad_GE21_L2, std::make_pair(gem_id, pad), event_.CSCTF_bx[j]);
                }                
              }
            }
          }
        }
      }
      // found pad is not empty
      if (bestPad_GE11_L1.second != GEMPadDigi()){
        std::cout << "Best pad GE11 L1" << bestPad_GE11_L1.second << std::endl;
        if (bestPad_GE11_L1.first.station()==1 and bestPad_GE11_L1.first.layer()==1) {
          auto gem_gp1 = getGlobalPointPad(bestPad_GE11_L1.first, bestPad_GE11_L1.second);
          event_.GE11_phi_L1[j] = gem_gp1.phi();
          event_.GE11_bx_L1[j] = bestPad_GE11_L1.second.bx();
          event_.GE11_ch_L1[j] = bestPad_GE11_L1.first.chamber();
          event_.GE11_z_L1[j] = gem_gp1.z();              
        }
      } else {
        std::cout << "No best pad GE11 L1" << std::endl;
      }
      
      if (bestPad_GE11_L2.second != GEMPadDigi()){
        std::cout << "Best pad GE11 L2" << bestPad_GE11_L2.second << std::endl;
        if (bestPad_GE11_L2.first.station()==1 and bestPad_GE11_L2.first.layer()==2) {
          auto gem_gp1 = getGlobalPointPad(bestPad_GE11_L2.first, bestPad_GE11_L2.second);
          event_.GE11_phi_L2[j] = gem_gp1.phi();
          event_.GE11_bx_L2[j] = bestPad_GE11_L2.second.bx();
          event_.GE11_ch_L2[j] = bestPad_GE11_L2.first.chamber();
          event_.GE11_z_L2[j] = gem_gp1.z();              
        }
      } else {
        std::cout << "No best pad GE11 L2" << std::endl;
      }
      // found pad is not empty
      if (bestPad_GE21_L1.second != GEMPadDigi()){
        std::cout << "Best pad GE21 L1" << bestPad_GE21_L1.second << std::endl;
        if (bestPad_GE21_L1.first.station()==3 and bestPad_GE21_L1.first.layer()==1) {
          auto gem_gp1 = getGlobalPointPad(bestPad_GE21_L1.first, bestPad_GE21_L1.second);
          event_.GE21_phi_L1[j] = gem_gp1.phi();
          event_.GE21_bx_L1[j] = bestPad_GE21_L1.second.bx();
          event_.GE21_ch_L1[j] = bestPad_GE21_L1.first.chamber();
          event_.GE21_z_L1[j] = gem_gp1.z();              
        }
      } else {
        std::cout << "No best pad GE21 L1" << std::endl;
      }
      
      if (bestPad_GE21_L2.second != GEMPadDigi()){
        std::cout << "Best pad GE21 L2" << bestPad_GE21_L2.second << std::endl;
        if (bestPad_GE21_L2.first.station()==3 and bestPad_GE21_L2.first.layer()==2) {
          auto gem_gp1 = getGlobalPointPad(bestPad_GE21_L2.first, bestPad_GE21_L2.second);
          event_.GE21_phi_L2[j] = gem_gp1.phi();
          event_.GE21_bx_L2[j] = bestPad_GE21_L2.second.bx();
          event_.GE21_ch_L2[j] = bestPad_GE21_L2.first.chamber();
          event_.GE21_z_L2[j] = gem_gp1.z();              
        }
      } else {
        std::cout << "No best pad GE21 L2" << std::endl;
      }
    } // check if match to pads
  } // loop on csctf tracks

  // Store the RPCb variables
  if(verbose) std::cout << "Number of l1MuRPCbs: " <<l1MuRPCbs.size() << std::endl;
  event_.nRPCb = l1MuRPCbs.size();
  if (processRPCb_) {
    for (unsigned int j=0; j<l1MuRPCbs.size(); ++j) { 
      auto track = l1MuRPCbs[j];
      
      event_.RPCb_pt[j] = muPtScale->getPtScale()->getLowEdge(track.pt_packed());
      event_.RPCb_eta[j] = muScales->getRegionalEtaScale(track.type_idx())->getCenter(track.eta_packed());
      event_.RPCb_phi[j] = normalizedPhi(muScales->getPhiScale()->getLowEdge(track.phi_packed()));
      event_.RPCb_bx[j] = track.bx();
      event_.RPCb_quality[j] = track.quality();
      
      if(verbose) {  
        std::cout << "pt " << event_.RPCb_pt[j]
                  << ", eta " << event_.RPCb_eta[j]
                  << ", phi " << event_.RPCb_phi[j]
                  << ", bx " << event_.RPCb_bx[j]
                  << ", quality " << event_.RPCb_quality[j]
                  << std::endl;
      }
      
      auto link = l1MuRPCbLinks[j];
      event_.RPCb_nStubs[j] = 0;
      for (unsigned int ii=1; ii<=link.nlayer(); ++ii){
        if (link.empty(ii)) continue;
        event_.RPCb_nStubs[j] += 1;
        double phi = getGlobalPhi(link.rawdetId(ii), link.strip(ii));
        auto detId = RPCDetId(link.rawdetId(ii));
        if(verbose) {  
          std::cout << "\t" << ii 
                    << ", RPCDetId " << detId 
                    << ", strip " << link.strip(ii) 
                    << ", bx " << link.bx(ii) 
                    << ", phi " << phi
                    << std::endl;
        } 
        switch(ii) {
        case 1:
          event_.RPCb_bx1[j] = link.bx(ii); 
          event_.RPCb_strip1[j] = link.strip(ii); 
          event_.RPCb_phi1[j] = phi;
          event_.RPCb_re1[j] = detId.region(); 
          event_.RPCb_ri1[j] = detId.ring(); 
          event_.RPCb_st1[j] = detId.station(); 
          event_.RPCb_se1[j] = detId.sector();
          event_.RPCb_la1[j] = detId.layer(); 
          event_.RPCb_su1[j] = detId.subsector(); 
          event_.RPCb_ro1[j] = detId.roll();
          break;
        case 2:
          event_.RPCb_bx2[j] = link.bx(ii); 
          event_.RPCb_strip2[j] = link.strip(ii); 
          event_.RPCb_phi2[j] = phi;
          event_.RPCb_re2[j] = detId.region(); 
          event_.RPCb_ri2[j] = detId.ring(); 
          event_.RPCb_st2[j] = detId.station(); 
          event_.RPCb_se2[j] = detId.sector();
          event_.RPCb_la2[j] = detId.layer(); 
          event_.RPCb_su2[j] = detId.subsector(); 
          event_.RPCb_ro2[j] = detId.roll();
          break;
        case 3:
          event_.RPCb_bx3[j] = link.bx(ii); 
          event_.RPCb_strip3[j] = link.strip(ii); 
          event_.RPCb_phi3[j] = phi;
          event_.RPCb_re3[j] = detId.region(); 
          event_.RPCb_ri3[j] = detId.ring(); 
          event_.RPCb_st3[j] = detId.station(); 
          event_.RPCb_se3[j] = detId.sector();
          event_.RPCb_la3[j] = detId.layer(); 
          event_.RPCb_su3[j] = detId.subsector(); 
          event_.RPCb_ro3[j] = detId.roll();
          break;
        case 4:
          event_.RPCb_bx4[j] = link.bx(ii); 
          event_.RPCb_strip4[j] = link.strip(ii); 
          event_.RPCb_phi4[j] = phi;
          event_.RPCb_re4[j] = detId.region(); 
          event_.RPCb_ri4[j] = detId.ring(); 
          event_.RPCb_st4[j] = detId.station(); 
          event_.RPCb_se4[j] = detId.sector();
          event_.RPCb_la4[j] = detId.layer(); 
          event_.RPCb_su4[j] = detId.subsector(); 
          event_.RPCb_ro4[j] = detId.roll();
          break;
        case 5:
          event_.RPCb_bx5[j] = link.bx(ii); 
          event_.RPCb_strip5[j] = link.strip(ii); 
          event_.RPCb_phi5[j] = phi;
          event_.RPCb_re5[j] = detId.region(); 
          event_.RPCb_ri5[j] = detId.ring(); 
          event_.RPCb_st5[j] = detId.station(); 
          event_.RPCb_se5[j] = detId.sector();
          event_.RPCb_la5[j] = detId.layer(); 
          event_.RPCb_su5[j] = detId.subsector(); 
          event_.RPCb_ro5[j] = detId.roll();
          break;
        case 6:
          event_.RPCb_bx6[j] = link.bx(ii); 
          event_.RPCb_strip6[j] = link.strip(ii); 
          event_.RPCb_phi6[j] = phi;
          event_.RPCb_re6[j] = detId.region(); 
          event_.RPCb_ri6[j] = detId.ring(); 
          event_.RPCb_st6[j] = detId.station(); 
          event_.RPCb_se6[j] = detId.sector();
          event_.RPCb_la6[j] = detId.layer(); 
          event_.RPCb_su6[j] = detId.subsector(); 
          event_.RPCb_ro6[j] = detId.roll();
          break;
        };
      } 
    }
  }


  // Store the RPCf variables
  if(verbose) std::cout << "Number of l1MuRPCfs: " <<l1MuRPCfs.size() << std::endl;
  event_.nRPCf = l1MuRPCfs.size();
  if (processRPCf_) {
    for (unsigned int j=0; j<l1MuRPCfs.size(); ++j) { 
      auto track = l1MuRPCfs[j];
      
      event_.RPCf_pt[j] = muPtScale->getPtScale()->getLowEdge(track.pt_packed());
      event_.RPCf_eta[j] = muScales->getRegionalEtaScale(track.type_idx())->getCenter(track.eta_packed());
      event_.RPCf_phi[j] = normalizedPhi(muScales->getPhiScale()->getLowEdge(track.phi_packed()));
      event_.RPCf_bx[j] = track.bx();
      event_.RPCf_quality[j] = track.quality();
      
      if(verbose) {  
        std::cout << "pt " << event_.RPCf_pt[j]
                  << ", eta " << event_.RPCf_eta[j]
                  << ", phi " << event_.RPCf_phi[j]
                  << ", bx " << event_.RPCf_bx[j]
                  << ", quality " << event_.RPCf_quality[j]
                  << std::endl;
      }
      
      auto link = l1MuRPCfLinks[j];
      event_.RPCf_nStubs[j] = 0;
      for (unsigned int ii=1; ii<=link.nlayer(); ++ii){
        if (link.empty(ii)) continue;
        event_.RPCf_nStubs[j] += 1;
        double phi = getGlobalPhi(link.rawdetId(ii), link.strip(ii));
        auto detId = RPCDetId(link.rawdetId(ii));
        if(verbose) {  
          std::cout << "\t" << ii
                    << ", RPCDetId " << detId 
                    << ", strip " << link.strip(ii) 
                    << ", bx " << link.bx(ii) 
                    << ", phi " << phi
                    << std::endl;
        } 
        switch(ii) {
        case 1:
          event_.RPCf_bx1[j] = link.bx(ii); 
          event_.RPCf_strip1[j] = link.strip(ii); 
          event_.RPCf_phi1[j] = phi;
          event_.RPCf_re1[j] = detId.region(); 
          event_.RPCf_ri1[j] = detId.ring(); 
          event_.RPCf_st1[j] = detId.station(); 
          event_.RPCf_se1[j] = detId.sector();
          event_.RPCf_la1[j] = detId.layer(); 
          event_.RPCf_su1[j] = detId.subsector(); 
          event_.RPCf_ro1[j] = detId.roll();
          break;
        case 2:
          event_.RPCf_bx2[j] = link.bx(ii); 
          event_.RPCf_strip2[j] = link.strip(ii); 
          event_.RPCf_phi2[j] = phi;
          event_.RPCf_re2[j] = detId.region(); 
          event_.RPCf_ri2[j] = detId.ring(); 
          event_.RPCf_st2[j] = detId.station(); 
          event_.RPCf_se2[j] = detId.sector();
          event_.RPCf_la2[j] = detId.layer(); 
          event_.RPCf_su2[j] = detId.subsector(); 
          event_.RPCf_ro2[j] = detId.roll();
          break;
        case 3:
          event_.RPCf_bx3[j] = link.bx(ii); 
          event_.RPCf_strip3[j] = link.strip(ii); 
          event_.RPCf_phi3[j] = phi;
          event_.RPCf_re3[j] = detId.region(); 
          event_.RPCf_ri3[j] = detId.ring(); 
          event_.RPCf_st3[j] = detId.station(); 
          event_.RPCf_se3[j] = detId.sector();
          event_.RPCf_la3[j] = detId.layer(); 
          event_.RPCf_su3[j] = detId.subsector(); 
          event_.RPCf_ro3[j] = detId.roll();
          break;
        case 4:
          event_.RPCf_bx4[j] = link.bx(ii); 
          event_.RPCf_strip4[j] = link.strip(ii); 
          event_.RPCf_phi4[j] = phi;
          event_.RPCf_re4[j] = detId.region(); 
          event_.RPCf_ri4[j] = detId.ring(); 
          event_.RPCf_st4[j] = detId.station(); 
          event_.RPCf_se4[j] = detId.sector();
          event_.RPCf_la4[j] = detId.layer(); 
          event_.RPCf_su4[j] = detId.subsector(); 
          event_.RPCf_ro4[j] = detId.roll();
          break;
        case 5:
          event_.RPCf_bx5[j] = link.bx(ii); 
          event_.RPCf_strip5[j] = link.strip(ii); 
          event_.RPCf_phi5[j] = phi;
          event_.RPCf_re5[j] = detId.region(); 
          event_.RPCf_ri5[j] = detId.ring(); 
          event_.RPCf_st5[j] = detId.station(); 
          event_.RPCf_se5[j] = detId.sector();
          event_.RPCf_la5[j] = detId.layer(); 
          event_.RPCf_su5[j] = detId.subsector(); 
          event_.RPCf_ro5[j] = detId.roll();
          break;
        case 6:
          event_.RPCf_bx6[j] = link.bx(ii); 
          event_.RPCf_strip6[j] = link.strip(ii); 
          event_.RPCf_phi6[j] = phi;
          event_.RPCf_re6[j] = detId.region(); 
          event_.RPCf_ri6[j] = detId.ring(); 
          event_.RPCf_st6[j] = detId.station(); 
          event_.RPCf_se6[j] = detId.sector();
          event_.RPCf_la6[j] = detId.layer(); 
          event_.RPCf_su6[j] = detId.subsector(); 
          event_.RPCf_ro6[j] = detId.roll();
          break;
        };
      } 
    }
  }


  // Store the L1Mu variables
  // also check which DTTF, CSCTF, RPCf, RPCb this L1Mu corresponds to!
  if(verbose) std::cout << "Number of L1Mu candidates before selections " << event_.nL1Mu << std::endl; 

  for (unsigned int i=0; i<l1GmtCands.size(); ++i) {
    auto l1Mu = l1GmtCands[i];
    event_.L1Mu_pt[i] = l1Mu.ptValue();
    event_.L1Mu_eta[i] = l1Mu.etaValue();
    event_.L1Mu_phi[i] = normalizedPhi(l1Mu.phiValue());
    event_.L1Mu_charge[i] = l1Mu.charge();
    event_.L1Mu_quality[i] = l1Mu.quality();
    event_.L1Mu_bx[i] = l1Mu.bx();
  
    if(verbose) {
      cout << "l1Mu " << i << endl; 
      cout << "l1Mu_pt " << event_.L1Mu_pt[i] << endl;
      cout << "l1Mu_eta " << event_.L1Mu_eta[i] << endl;
      cout << "l1Mu_phi " << event_.L1Mu_phi[i] << endl;
      cout << "l1Mu_quality " << event_.L1Mu_quality[i] << endl;
      cout << "l1Mu_charge " << event_.L1Mu_charge[i] << endl;
      cout << "l1Mu_bx " << event_.L1Mu_bx[i] << endl;
    }

    // Matching to DTTF
    double bestDrL1MuL1DTTrack = 99;
    for (unsigned int j=0; j<L1DTTrackPhis.size(); ++j) {   
      if ( ( event_.L1Mu_quality[i] > 0 ) &&
           ( reco::deltaPhi( event_.L1Mu_phi[i], event_.DTTF_phi[j] ) < 0.001 ) &&             
           ( event_.L1Mu_bx[i] == event_.DTTF_bx[j] ) ) {
        double drL1MuL1DTTrack = reco::deltaR(l1Mu.etaValue(), 
                                              normalizedPhi(l1Mu.phiValue()), 
                                              event_.DTTF_eta[j], 
                                              event_.DTTF_phi[j]);
        if (drL1MuL1DTTrack < bestDrL1MuL1DTTrack and drL1MuL1DTTrack < 0.3) {
          bestDrL1MuL1DTTrack = drL1MuL1DTTrack;
          event_.L1Mu_DTTF_index[i] = j;
        }
      }
    }
        
    if(verbose) {  
      int tempIndex = event_.L1Mu_DTTF_index[i]; 
      if (tempIndex != -1) { // and bestDrL1MuL1DTTrack < 0.2
        // Print matching DTTF track
        std::cout << "\tMatching DTTF track" << std::endl;
        std::cout << "\tpt = "  << event_.DTTF_pt[tempIndex]
                  << ", eta = " << event_.DTTF_eta[tempIndex]
                  << ", phi = " << event_.DTTF_phi[tempIndex] 
                  << ", bx = "  << event_.DTTF_bx[tempIndex]
                  << ", quality = " << event_.DTTF_quality[tempIndex] << std::endl;
        
        // Print stubs
        std::cout << "\tNumber of stubs: " << event_.DTTF_nStubs[tempIndex] << std::endl; 
        for (auto stub: L1DTTrackPhis[tempIndex].second) {
          std::cout << "\t\t " << stub << std::endl;
          std::cout << "\t\t phiValue = " << stub.phiValue() << ", phibValue = " << stub.phibValue() << std::endl;
        }  
      }
      else {
        std::cout << "\tNo matching DTTF track" << std::endl;
      }
    }
    
    // Matching to CSCTF     
    double bestDrL1MuL1CSCTrack = 99;
    for (unsigned int j=0; j<l1Tracks.size(); ++j) {
      if ( ( event_.L1Mu_quality[i] > 0 ) &&
           ( reco::deltaPhi( event_.L1Mu_phi[i], event_.CSCTF_phi[j] ) < 0.001 ) &&             
           ( event_.L1Mu_bx[i] == event_.CSCTF_bx[j] ) ) {
        double drL1MuL1CSCTrack = reco::deltaR(l1Mu.etaValue(), 
                                               normalizedPhi(l1Mu.phiValue()), 
                                               event_.CSCTF_eta[j], 
                                               event_.CSCTF_phi[j]);
        if (drL1MuL1CSCTrack < bestDrL1MuL1CSCTrack and drL1MuL1CSCTrack < 0.3) {
          bestDrL1MuL1CSCTrack = drL1MuL1CSCTrack;
          event_.L1Mu_CSCTF_index[i] = j;
        }
      }                
    }
    

    if(verbose) {  
      int tempIndex = event_.L1Mu_CSCTF_index[i]; 
      if (tempIndex != -1) { // and bestDrL1MuL1CSCTrack < 0.2
        // Print matching CSCTF track
        std::cout << "\tMatching CSCTF track" << std::endl;
        std::cout << "\tpt = "  << event_.CSCTF_pt[tempIndex]
                  << ", eta = " << event_.CSCTF_eta[tempIndex]
                  << ", phi = " << event_.CSCTF_phi[tempIndex]
                  << ", bx = "  << event_.CSCTF_bx[tempIndex]
                  << ", quality = " << event_.CSCTF_quality[tempIndex]
                  << std::endl;
        //printCSCStubProperties(tempIndex);

        // // Print stubs
        // std::cout << "\tNumber of stubs: " << event_.CSCTF_nStubs[tempIndex] << std::endl;
        // auto stubCollection(l1Tracks[tempIndex].second);
        // for (auto detUnitIt = stubCollection.begin(); detUnitIt != stubCollection.end(); detUnitIt++) {
        //   const CSCDetId& id = (*detUnitIt).first;
        //   std::cout << "\t\tDetId " << id << std::endl;
        //   const auto range = (*detUnitIt).second;
        //   for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        //     //if (!(*digiIt).isValid()) continue;
        //     GlobalPoint csc_gp = getCSCSpecificPoint2(id.rawId(), *digiIt);
        //     std::cout << "\t\tPosition " << csc_gp.eta() << " " << csc_gp.phi() << std::endl;
        //     std::cout << "\t\t" << *digiIt << std::endl;
        //   }
        // }
      }
      else {
        std::cout << "\tNo matching CSCTF track" << std::endl;
      }
    }

    /* 
       CSCTF stub recovery
       The CSC track-finder may drop certain stubs if they don't match the pattern
       First get the station numbers where stubs are not filled
    */
    // std::cout << "CSC stub recovery" << std::endl;
    if (event_.L1Mu_CSCTF_index[i] != -1 and false) {
      bool stubMissingSt1 = event_.CSCTF_st1[ event_.L1Mu_CSCTF_index[i] ] == 99;
      bool stubMissingSt2 = event_.CSCTF_st2[ event_.L1Mu_CSCTF_index[i] ] == 99;
      bool stubMissingSt3 = event_.CSCTF_st3[ event_.L1Mu_CSCTF_index[i] ] == 99;
      bool stubMissingSt4 = event_.CSCTF_st4[ event_.L1Mu_CSCTF_index[i] ] == 99;
      bool doStubRecovery = stubMissingSt1 or stubMissingSt2 or stubMissingSt3 or stubMissingSt4;
      
      std::vector<float> xs;
      std::vector<float> ys;
      std::vector<float> zs;
      std::cout << "Get stub positions" << std::endl;
      getStubPositions(event_.L1Mu_CSCTF_index[i], xs, ys, zs);
      
      std::cout << "fit stub positions with straight line" << std::endl;
      float alpha_x, beta_x, alpha_y, beta_y;
      fitStraightLine(zs, xs, alpha_x, beta_x); 
      fitStraightLine(zs, ys, alpha_y, beta_y); 

      std::cout << "Get positions stations" << std::endl;
      std::vector<float> allxs;
      std::vector<float> allys;
      int sign_z = int(event_.L1Mu_eta[i]/std::abs(event_.L1Mu_eta[i]));
      getPositionsStations(alpha_x, beta_x, alpha_y, beta_y,
                           allxs, allys, sign_z);

      event_.CSCTF_fitline_x1[ event_.L1Mu_CSCTF_index[i] ] = allxs[0];
      event_.CSCTF_fitline_x2[ event_.L1Mu_CSCTF_index[i] ] = allxs[1];
      event_.CSCTF_fitline_x3[ event_.L1Mu_CSCTF_index[i] ] = allxs[2];
      event_.CSCTF_fitline_x4[ event_.L1Mu_CSCTF_index[i] ] = allxs[3];

      event_.CSCTF_fitline_y1[ event_.L1Mu_CSCTF_index[i] ] = allys[0];
      event_.CSCTF_fitline_y2[ event_.L1Mu_CSCTF_index[i] ] = allys[1];
      event_.CSCTF_fitline_y3[ event_.L1Mu_CSCTF_index[i] ] = allys[2];
      event_.CSCTF_fitline_y4[ event_.L1Mu_CSCTF_index[i] ] = allys[3];

      if (doStubRecovery and event_.L1Mu_CSCTF_index[i] != -1) {
        int triggerSector = (l1Tracks[ event_.L1Mu_CSCTF_index[i] ].first).sector();
        std::cout << "trigger sector " << triggerSector << std::endl;
        
        for (int endcap=1; endcap<=2; endcap++){
          //do not consider stubs in the wrong endcap
          int zendcap(endcap!=1 ? -1 : +1 );
          if (zendcap * event_.L1Mu_eta[i] < 0) continue;
          for (int station=1; station<=4; station++){
            
            // ignore station where a L1Mu stub is present!
            if (not stubMissingSt1 and station==1) continue;
            if (not stubMissingSt2 and station==2) continue;
            if (not stubMissingSt3 and station==3) continue;
            if (not stubMissingSt4 and station==4) continue;
            std::cout << "Attempt recovery  in station " << station << std::endl;
            // temp storage of candidate stubs per station and ring
            CSCCorrelatedLCTDigiId bestMatchingStub;
            int iStub = 0;
            for (int ring=1; ring<=3; ring++){
              if (station!=1 and ring==3) continue;
              std::cout << "Analyzing ring " << ring << std::endl;
              
              for (int chamber=1; chamber<=36; chamber++){
                // do not consider invalid detids
                if ( (station==2 or station==3 or station==4) and 
                     (ring==1) and chamber>18) continue;
                //std::cout << "Analyzing chamber " << chamber << std::endl; 
                // create the detid
                CSCDetId ch_id(endcap, station, ring, chamber);
                //std::cout << "ch_id " <<  ch_id << std::endl;
                // get the stubs in this detid
                auto range = CSCCorrelatedLCTs.get(ch_id);
                for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
                  iStub++; 
                  // trigger sector must be the same
                  if (triggerSector != ch_id.triggerSector()) continue;
                  auto stub(*digiItr);
                  int deltaBX = std::abs(stub.getBX() - (6 + event_.L1Mu_bx[i]));
                  
                  // BXs have to match
                  if (deltaBX > 1) continue;
                  std::cout << ch_id << std::endl;
                  std::cout<<"Candidate " << stub << std::endl;
                  bestMatchingStub = pickBestMatchingStub(allxs[ch_id.station()-1], allys[ch_id.station()-1], 
                                                          bestMatchingStub, std::make_pair(ch_id, stub), 6 + event_.L1Mu_bx[i]);
                }
                // consider the case ME1a
                if (station==1 and ring==1){
                  CSCDetId me1a_id(endcap, station, 4, chamber);
                  auto range = CSCCorrelatedLCTs.get(me1a_id);
                  for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
                    iStub++;
                    // trigger sector must be the same
                    if (triggerSector != me1a_id.triggerSector()) continue;
                    auto stub(*digiItr);
                    int deltaBX = std::abs(stub.getBX() - (6 + event_.L1Mu_bx[i]));
                    
                    // BXs have to match
                    if (deltaBX > 1) continue; 
                    std::cout << me1a_id << std::endl;
                    std::cout<<"Candidate " << stub << std::endl;
                    bestMatchingStub = pickBestMatchingStub(allxs[me1a_id.station()-1], allys[me1a_id.station()-1], 
                                                            bestMatchingStub, std::make_pair(me1a_id, stub), 6 + event_.L1Mu_bx[i]);
                    
                  }
                }
                //std::cout << "Current best: " << bestMatchingStub.second << std::endl;
              }
              //std::cout << "Current best2: " << bestMatchingStub.second << std::endl;
            }
            if (bestMatchingStub.second != CSCCorrelatedLCTDigi()) {
              std::cout << "Best matching stub " << bestMatchingStub.first << " " << bestMatchingStub.second <<std::endl;
            
              // stub position
              auto gp = getCSCSpecificPoint2(bestMatchingStub.first.rawId(), bestMatchingStub.second);

              float z_pos_L3, bestFitPhi, bestFitDPhi;
              fitComparatorsLCT(*hCSCComparators.product(), bestMatchingStub.second, bestMatchingStub.first, 0, z_pos_L3, bestFitPhi, bestFitDPhi);
              
              fillCSCStubProperties(bestMatchingStub.first, bestMatchingStub.second, event_.L1Mu_CSCTF_index[i], 
                                    gp, z_pos_L3, bestFitPhi, bestFitDPhi);
            }
            else{
              if (iStub!=0) std::cout << "No best matching stub " << std::endl;
            }
          }
        }
      }
    }

    if (processRPCb_) {
      // Matching to RPCb cands
      double bestDrL1MuL1RPCb = 99;
      for (unsigned int j=0; j<l1MuRPCbs.size(); ++j) { 
        if ( ( event_.L1Mu_quality[i] > 0 ) &&
             ( reco::deltaPhi( event_.L1Mu_phi[i], event_.RPCb_phi[j] ) < 0.001 ) &&             
             ( event_.L1Mu_bx[i] == event_.RPCb_bx[j] ) ) {
          double drL1MuL1RPCb = reco::deltaR(l1Mu.etaValue(), 
                                             normalizedPhi(l1Mu.phiValue()), 
                                             event_.RPCb_eta[j], 
                                             event_.RPCb_phi[j]);
          if (drL1MuL1RPCb < bestDrL1MuL1RPCb and drL1MuL1RPCb < 0.3) {
            bestDrL1MuL1RPCb = drL1MuL1RPCb;
            event_.L1Mu_RPCb_index[i] = j;
          }
        }                
      }
      
      if(verbose) {  
        int tempIndex = event_.L1Mu_RPCb_index[i]; 
        if (tempIndex != -1) { // and bestDrL1MuL1CSCTrack < 0.2
          // Print matching RPCb track
          std::cout << "\tMatching RPCb track" << std::endl;
          std::cout << "\tpt = "  << event_.RPCb_pt[tempIndex]
                    << ", eta = " << event_.RPCb_eta[tempIndex]
                    << ", phi = " << event_.RPCb_phi[tempIndex]
                    << ", bx = "  << event_.RPCb_bx[tempIndex]
                    << ", quality = " << event_.RPCb_quality[tempIndex]
                    << ", nStubs = " << event_.RPCb_nStubs[tempIndex]
                    << std::endl;
        }
        else {
          std::cout << "\tNo matching RPCb track" << std::endl;
        }
      }
    }
    
    if (processRPCf_) {
      // Matching to RPCf cands
      double bestDrL1MuL1RPCf = 99;
      for (unsigned int j=0; j<l1MuRPCfs.size(); ++j) { 
        if ( ( event_.L1Mu_quality[i] > 0 ) &&
             ( reco::deltaPhi( event_.L1Mu_phi[i], event_.RPCf_phi[j] ) < 0.001 ) &&             
             ( event_.L1Mu_bx[i] == event_.RPCf_bx[j] ) ) {
          double drL1MuL1RPCf = reco::deltaR(l1Mu.etaValue(), 
                                             normalizedPhi(l1Mu.phiValue()), 
                                             event_.RPCf_eta[j], 
                                             event_.RPCf_phi[j]);
          if (drL1MuL1RPCf < bestDrL1MuL1RPCf and drL1MuL1RPCf < 0.3) {
            bestDrL1MuL1RPCf = drL1MuL1RPCf;
            event_.L1Mu_RPCf_index[i] = j;
          }
        }                
      }
      
      if(verbose) {  
        int tempIndex = event_.L1Mu_RPCf_index[i]; 
        if (tempIndex != -1) { // and bestDrL1MuL1CSCTrack < 0.2
          // Print matching RPCf track
          std::cout << "\tMatching RPCf track" << std::endl;
          std::cout << "\tpt = "  << event_.RPCf_pt[tempIndex]
                    << ", eta = " << event_.RPCf_eta[tempIndex]
                    << ", phi = " << event_.RPCf_phi[tempIndex]
                    << ", bx = "  << event_.RPCf_bx[tempIndex]
                    << ", quality = " << event_.RPCf_quality[tempIndex]
                    << ", nStubs = " << event_.RPCf_nStubs[tempIndex]
                    << std::endl;
        }
        else {
          std::cout << "\tNo matching RPCf track" << std::endl;
        }
      }
    }
    
    // // calculate the number of L1Tk within 0.12
    // for (unsigned int j=0; j<TTTracks.size(); ++j) {
    //   auto l1Tk = TTTracks[j];
    //   const double l1Tk_pt = l1Tk.getMomentum().perp();
    //   const double l1Tk_eta = l1Tk.getMomentum().eta();
    //   const double l1Tk_phi = normalizedPhi(l1Tk.getMomentum().phi());
    //   const double l1Tk_charge = l1Tk.getRInv()>0? 1: -1;
    //   const double l1Tk_eta_corr = l1Tk_eta;
    //   const double l1Tk_phi_corr = phiHeavyCorr(l1Tk_pt, l1Tk_eta, l1Tk_phi, l1Tk_charge);

    //   if(verbose and false) {
    //     cout << "l1Tk " << j << endl; 
    //     cout << "l1Tk_pt " << l1Tk_pt << endl;
    //     cout << "l1Tk_eta " << l1Tk_eta << endl;
    //     cout << "l1Tk_phi " << l1Tk_phi << endl;
    //     cout << "l1Tk_phi_corr " << l1Tk_phi_corr << endl;
    //     cout << "l1Tk_charge " << l1Tk_charge << endl;
    //   }
      
    //   double l1Tk_eta_prop = -99;
    //   double l1Tk_phi_prop = -99;
    //   GlobalPoint ex_point(extrapolateGP(l1Tk));
    //   if (!(ex_point == GlobalPoint())) {
    //     l1Tk_eta_prop = ex_point.eta();
    //     l1Tk_phi_prop = ex_point.phi();
    //     if(verbose and false) {
    //       cout << "l1Tk_eta_prop " << l1Tk_eta_prop << endl;
    //       cout << "l1Tk_phi_prop " << l1Tk_phi_prop << endl;
    //     }
    //     const double dR_l1Mu_l1Tk_prop = reco::deltaR(l1Tk_eta_prop, l1Tk_phi_prop, event_.L1Mu_eta[i], event_.L1Mu_phi[i]);
    //     if (dR_l1Mu_l1Tk_prop < event_.L1Mu_L1Tk_dR_prop[i]) {
    //       event_.L1Mu_L1Tk_dR_prop[i] = dR_l1Mu_l1Tk_prop;
    //       event_.L1Mu_L1Tk_pt_prop[i] = l1Tk_pt;
    //     }
    //   }

    //   // TrajectoryStateOnSurface stateAtMB2 = extrapolate(l1Tk);
    //   // if (stateAtMB2.isValid()) {
    //   //   // std::cout << ">>>Final State is valid" << std::endl;
    //   //   l1Tk_eta_prop = stateAtMB2.globalPosition().eta();
    //   //   l1Tk_phi_prop = stateAtMB2.globalPosition().phi();
    //   //   if(verbose) {
    //   //     cout << "l1Tk_eta_prop " << l1Tk_eta_prop << endl;
    //   //     cout << "l1Tk_phi_prop " << l1Tk_phi_prop << endl;
    //   //   }
    //   //   // cout << endl;
    //   //   const double dR_l1Mu_l1Tk_prop = reco::deltaR(l1Tk_eta_prop, l1Tk_phi_prop, event_.L1Mu_eta[i], event_.L1Mu_phi[i]);
    //   //   if (dR_l1Mu_l1Tk_prop < event_.L1Mu_L1Tk_dR_prop[i]) {
    //   //     event_.L1Mu_L1Tk_dR_prop[i] = dR_l1Mu_l1Tk_prop;
    //   //     event_.L1Mu_L1Tk_pt_prop[i] = l1Tk_pt;
    //   //   }
    //   // }
      
    //   const double dR_l1Mu_l1Tk_corr = reco::deltaR(l1Tk_eta_corr, l1Tk_phi_corr, event_.L1Mu_eta[i], event_.L1Mu_phi[i]);
    //   if (dR_l1Mu_l1Tk_corr < event_.L1Mu_L1Tk_dR_corr[i]) {
    //     event_.L1Mu_L1Tk_dR_corr[i] = dR_l1Mu_l1Tk_corr;
    //     event_.L1Mu_L1Tk_pt_corr[i] = l1Tk_pt;
    //   }
    // } // end of loop on TTTracks      
  }
  
  event_tree_->Fill();  
  
  return true;
}


// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedL1MuFilter::beginJob()
{
}
 
// ------------ method called once each job just after ending the event loop  ------------
 void 
   DisplacedL1MuFilter::endJob()
{
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedL1MuFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
DisplacedL1MuFilter::getStubPositions(int index, 
                                      std::vector<float>& x, 
                                      std::vector<float>& y, 
                                      std::vector<float>& z) const
{
  x.clear();
  y.clear();
  z.clear();
  if (event_.CSCTF_st1[index] == 1){
    x.push_back(event_.CSCTF_x1[index]);
    y.push_back(event_.CSCTF_y1[index]);
    z.push_back(event_.CSCTF_z1[index]);
  }
  if (event_.CSCTF_st2[index] == 2){
    x.push_back(event_.CSCTF_x2[index]);
    y.push_back(event_.CSCTF_y2[index]);
    z.push_back(event_.CSCTF_z2[index]);
  }
  if (event_.CSCTF_st3[index] == 3){
    x.push_back(event_.CSCTF_x3[index]);
    y.push_back(event_.CSCTF_y3[index]);
    z.push_back(event_.CSCTF_z3[index]);
  }
  if (event_.CSCTF_st4[index] == 4){
    x.push_back(event_.CSCTF_x4[index]);
    y.push_back(event_.CSCTF_y4[index]);
    z.push_back(event_.CSCTF_z4[index]);
  }
} 


void 
DisplacedL1MuFilter::fitComparatorsLCT(const CSCComparatorDigiCollection& hCSCComparators,
                                       const CSCCorrelatedLCTDigi& stub, CSCDetId ch_id, 
                                       int iMuon, float& fit_z, float& fit_phi, float& fit_dphi) const
{
  bool verbose = false;
  
  auto cscChamber = cscGeometry_->chamber(ch_id);
  
  // fetch the CSC comparator digis in this chamber
  CSCComparatorDigiContainerIds compDigisIds;
  for (int iLayer=1; iLayer<=6; ++iLayer){
    CSCDetId layerId(ch_id.endcap(), ch_id.station(), ch_id.ring(), ch_id.chamber(), iLayer);
    // get the digis per layer
    auto compRange = hCSCComparators.get(layerId);
    CSCComparatorDigiContainer compDigis;
    for (auto compDigiItr = compRange.first; compDigiItr != compRange.second; compDigiItr++) {
      auto compDigi = *compDigiItr;
      //if (stub.getTimeBin() < 4 or stub.getTimeBin() > 8) continue;
      int stubHalfStrip(getHalfStrip(compDigi));
      // these comparator digis never fit the pattern anyway!
      if (std::abs(stubHalfStrip-stub.getStrip())>5) continue;
      // check if this comparator digi fits the pattern
      //if(verbose) std::cout << "Comparator digi L1Mu " << layerId << " " << compDigi << " HS " << stubHalfStrip << " stubKeyHS " << stub.getStrip() << std::endl; 
      if (comparatorInLCTPattern(stub.getStrip(), stub.getPattern(), iLayer, stubHalfStrip)) {
        //if(verbose) std::cout<<"\tACCEPT"<<std::endl;
        compDigis.push_back(compDigi);
      }
      // else{
      //   if(verbose) std::cout<<"\tDECLINE!"<<std::endl;
      // }
    }
    // if(verbose) if (compDigis.size() > 2) std::cout << ">>> INFO: " << compDigis.size() << " matched comp digis in this layer!" << std::endl;
    compDigisIds.push_back(std::make_pair(layerId, compDigis));
  }
  
  // get the z and phi positions
  std::vector<float> phis;
  std::vector<float> zs;
  std::vector<float> ephis;
  std::vector<float> ezs;
  for (auto p: compDigisIds){
    auto detId = p.first;
    for (auto hit: p.second){
      float fractional_strip = getFractionalStrip(hit);
      auto layer_geo = cscChamber->layer(detId.layer())->geometry();
      LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, 20);
      GlobalPoint csc_gp = cscGeometry_->idToDet(detId)->surface().toGlobal(csc_intersect);
      zs.push_back(csc_gp.z());
      ezs.push_back(0);
      phis.push_back(csc_gp.phi());
      ephis.push_back(0.);
      //ephis.push_back(gemvalidation::cscHalfStripWidth(detId)/sqrt(12));
    }
  }
  
  // do a fit to the comparator digis
  float alpha = 0., beta = 0.;
  fitStraightLineErrors(zs, phis, ezs, ephis,
                        alpha, beta, 
                        event_.lumi, event_.run, event_.event, iMuon, ch_id.station(), false);
  
  fit_z = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->centerOfStrip(20).z();
  fit_phi = alpha + beta * fit_z;
  
  if(verbose) {
    std::cout << "Number of comparator digis used in the fit " << ezs.size() << std::endl;
    std::cout << "best CSC stub fit phi position (L1Only) " << fit_z << " " << fit_phi << std::endl;
  }

  // calculate the z position in L1 and L6
  float l1_z = cscChamber->layer(1)->centerOfStrip(20).z();
  float l6_z = cscChamber->layer(6)->centerOfStrip(20).z();
  // calculate the bending angle
  fit_dphi = beta*(l6_z-l1_z);
}


float 
DisplacedL1MuFilter::getGlobalPhi(unsigned int rawid, int stripN)
{  
  const RPCDetId id(rawid);
  std::unique_ptr<const RPCRoll>  roll(rpcGeometry_->roll(id));
  const uint16_t strip = stripN;
  const LocalPoint lp = roll->centreOfStrip(strip);
  const GlobalPoint gp = roll->toGlobal(lp);
  roll.release();
  return gp.phi();
}

double 
DisplacedL1MuFilter::calcCSCSpecificPhi(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const
{
  return getCSCSpecificPoint2(rawId, lct).phi();
}

GlobalPoint 
DisplacedL1MuFilter::getGlobalPointPad(unsigned int rawId, const GEMPadDigi& tp) const
{
  GEMDetId gem_id(rawId);
  // cout << "old " << gem_id << endl;
  // if (gem_id.station()==2){
  //   gem_id = GEMDetId(gem_id.region(), gem_id.ring(), 3, gem_id.layer(), gem_id.chamber(), gem_id.roll());
  // }
  // cout << "new " << gem_id << endl;
  LocalPoint gem_lp = gemGeometry_->etaPartition(gem_id)->centreOfPad(tp.pad());
  GlobalPoint gem_gp = gemGeometry_->idToDet(gem_id)->surface().toGlobal(gem_lp);
  return gem_gp;
}


GlobalPoint
DisplacedL1MuFilter::getCSCSpecificPoint2(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const 
{
  // taken from https://github.com/cms-sw/cmssw/blob/dc9f78b6af4ad56c9342cf14041b6485a60b0691/L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME11GEM.cc
  CSCDetId cscId = CSCDetId(rawId);
  CSCDetId key_id(cscId.endcap(), cscId.station(), cscId.ring(), cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);
  
  auto cscChamber = cscGeometry_->chamber(cscId);
  // "strip" here is actually a half-strip in geometry's terms
  // note that LCT::getStrip() starts from 0
  float fractional_strip = 0.5 * (lct.getStrip() + 1) - 0.25;
  auto layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();
  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  GlobalPoint csc_gp = cscGeometry_->idToDet(key_id)->surface().toGlobal(csc_intersect);
  //std::cout << "\t\t>>> other CSC LCT phi " << csc_gp.phi() << std::endl;
  //return getCSCSpecificPoint(rawId, lct).phi();
  return csc_gp;
}

CSCCorrelatedLCTDigiId 
DisplacedL1MuFilter::pickBestMatchingStub(float xref, float yref,
                                          const CSCCorrelatedLCTDigiId& oldStub,
                                          const CSCCorrelatedLCTDigiId& newStub, 
                                          int refBx) const
{
  bool debug=false;

  if (debug){
  std::cout << "In function pickBestMatchingStub" << std::endl;
  std::cout << "candidate 1 " << oldStub.first << " "  << oldStub.second << std::endl;
  std::cout << "candidate 2 " << newStub.first << " "  << newStub.second << std::endl;
  }

  // check for invalid/valid
  if (oldStub.second == CSCCorrelatedLCTDigi()) {
    if (debug) cout<<"Old stub invalid"<<endl;
    return newStub;
  }
  if (newStub.second == CSCCorrelatedLCTDigi()) {
    if (debug) cout<< "New stub invalid"<<endl;
    return oldStub;
  }
  
  int deltaBXOld = std::abs(oldStub.second.getBX() - refBx);
  int deltaBXNew = std::abs(newStub.second.getBX() - refBx);
  
  if (deltaBXOld==0 and deltaBXNew!=0) {
    if (debug) cout<<"Old stub in time, new stub out of time"<<endl;
    return oldStub;
  }
  if (deltaBXNew==0 and deltaBXOld!=0) {
    if (debug) cout<<"New stub in time, old stub out of time"<<endl;
    return newStub;
  }
  if ( (deltaBXOld!=0 and deltaBXNew!=0) or 
       (deltaBXOld==0 and deltaBXNew==0) ){
    if (debug) cout<<"Both stubs out of time"<<endl;
    // pick the one with the better matching wiregroup and halfstrip
    auto gpOld = getCSCSpecificPoint2(oldStub.first.rawId(), oldStub.second);
    auto gpNew = getCSCSpecificPoint2(newStub.first.rawId(), newStub.second);
    float deltaXYOld = TMath::Sqrt( (xref-gpOld.x())*(xref-gpOld.x()) + (yref-gpOld.y())*(yref-gpOld.y()) );
    float deltaXYNew = TMath::Sqrt( (xref-gpNew.x())*(xref-gpNew.x()) + (yref-gpNew.y())*(yref-gpNew.y()) );
    if (debug) {
      cout<<"xref "<< xref << " yref " << yref << endl;
      cout <<"gpOld " << gpOld << " gpNew " << gpNew << endl; 
      cout << "deltaXYOld " << deltaXYOld << " deltaXYNew " << deltaXYNew <<endl;
    }  
    if (deltaXYOld < deltaXYNew) {
      if (debug) cout<<"Old  stub better matching XY"<<endl;
      return oldStub;
    }
    else {
      if (debug) cout<<"New stub better matching XY"<<endl;
      return newStub;
    }
  }
  std::cout << "All else fails" << std::endl;
  // in case all else fails...
  return CSCCorrelatedLCTDigiId();
}

bool 
DisplacedL1MuFilter::stubInCSCTFTracks(const CSCCorrelatedLCTDigi& candidateStub, const L1CSCTrackCollection& l1Tracks) const
{
  bool isMatched = false;
  for (auto tftrack: l1Tracks){
    auto stubCollection = tftrack.second;
    for (auto detUnitIt = stubCollection.begin(); detUnitIt != stubCollection.end(); detUnitIt++) {
      const auto range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        //if (!(*digiIt).isValid()) continue;
        auto stub = *digiIt;
        if (candidateStub == stub) { 
          isMatched = true;
          break;
        }
      }
    }
  }
  return isMatched;
}


GEMPadDigiId
DisplacedL1MuFilter::pickBestMatchingCoPad(float xref, float yref,
                                           const GEMPadDigiId& oldCoPad,
                                           const GEMPadDigiId& newCoPad, 
                                           int bxref) const
{
  bool debug = true;
  // check for invalid/valid
  if (oldCoPad.second == GEMPadDigi()) {
    if (debug) cout<<"Old copad invalid"<<endl;
    return newCoPad;
  }
  if (newCoPad.second == GEMPadDigi()) {
    if (debug) cout<< "New copad invalid"<<endl;
    return oldCoPad;
  }
  
  // check the timing
  bool oldCoPadInTime = oldCoPad.second.bx() - bxref==0;// and oldCoPad.second.second().bx() - bxref==0;
  bool newCoPadInTime = newCoPad.second.bx() - bxref==0;// and newCoPad.second.second().bx() - bxref==0;
  if (oldCoPadInTime and not newCoPadInTime) {
    if (debug) cout<<"Old copad in time, new copad out of time"<<endl;
    return oldCoPad;
  }
  if (newCoPadInTime and not oldCoPadInTime) {
    if (debug) cout<<"New copad in time, old copad out of time"<<endl;
    return newCoPad;
  }

  // both copads in time, check the closest matching one!
  if ((oldCoPadInTime and newCoPadInTime) or (not oldCoPadInTime and not newCoPadInTime)){
    cout << "check better matching one in space" << endl;
    auto gpOld = getGlobalPointPad(oldCoPad.first.rawId(), oldCoPad.second);
    auto gpNew = getGlobalPointPad(newCoPad.first.rawId(), newCoPad.second);
    float deltaXYOld = TMath::Sqrt( (xref-gpOld.x())*(xref-gpOld.x()) + (yref-gpOld.y())*(yref-gpOld.y()) );
    float deltaXYNew = TMath::Sqrt( (xref-gpNew.x())*(xref-gpNew.x()) + (yref-gpNew.y())*(yref-gpNew.y()) );
    if (debug) {
      cout<<"xref "<< xref << " yref " << yref << endl;
      cout <<"gpOld " << gpOld << " gpNew " << gpNew << endl; 
      cout << "deltaXYOld " << deltaXYOld << " deltaXYNew " << deltaXYNew <<endl;
    }  
    if (deltaXYOld < deltaXYNew) {
      if (debug) cout<<"Old  copad better matching XY"<<endl;
      return oldCoPad;
    }
    else {
      if (debug) cout<<"New copad better matching XY"<<endl;
      return newCoPad;
    }
  }
  // in case all else fails...
  return GEMPadDigiId();
}


GEMPadDigiId
DisplacedL1MuFilter::pickBestMatchingPad(float xref, float yref,
                                           const GEMPadDigiId& oldPad,
                                           const GEMPadDigiId& newPad, 
                                           int bxref) const
{
  bool debug = false;
  // check for invalid/valid
  if (oldPad.second == GEMPadDigi()) {
    if (debug) cout<<"Old pad invalid"<<endl;
    return newPad;
  }
  if (newPad.second == GEMPadDigi()) {
    if (debug) cout<< "New pad invalid"<<endl;
    return oldPad;
  }
  
  // check the timing
  bool oldPadInTime = oldPad.second.bx() - bxref==0 and oldPad.second.bx() - bxref==0;
  bool newPadInTime = newPad.second.bx() - bxref==0 and newPad.second.bx() - bxref==0;
  if (oldPadInTime and not newPadInTime) {
    if (debug) cout<<"Old copad in time, new copad out of time"<<endl;
    return oldPad;
  }
  if (newPadInTime and not oldPadInTime) {
    if (debug) cout<<"New copad in time, old copad out of time"<<endl;
    return newPad;
  }

  // both copads in time, check the closest matching one!
  if ((oldPadInTime and newPadInTime) or (not oldPadInTime and  not newPadInTime)){
    cout << "check better matching one in space" << endl;
    auto gpOld = getGlobalPointPad(oldPad.first.rawId(), oldPad.second);
    auto gpNew = getGlobalPointPad(newPad.first.rawId(), newPad.second);
    float deltaXYOld = TMath::Sqrt( (xref-gpOld.x())*(xref-gpOld.x()) + (yref-gpOld.y())*(yref-gpOld.y()) );
    float deltaXYNew = TMath::Sqrt( (xref-gpNew.x())*(xref-gpNew.x()) + (yref-gpNew.y())*(yref-gpNew.y()) );
    if (debug) {
      cout<<"xref "<< xref << " yref " << yref << endl;
      cout <<"gpOld " << gpOld << " gpNew " << gpNew << endl; 
      cout << "deltaXYOld " << deltaXYOld << " deltaXYNew " << deltaXYNew <<endl;
    }  
    if (deltaXYOld < deltaXYNew) {
      if (debug) cout<<"Old  copad better matching XY"<<endl;
      return oldPad;
    }
    else {
      if (debug) cout<<"New copad better matching XY"<<endl;
      return newPad;
    }
  }
  // in case all else fails...
  return GEMPadDigiId();
}


void 
DisplacedL1MuFilter::fillCSCStubProperties(const CSCDetId& ch_id,
                                           const CSCCorrelatedLCTDigi& stub,
                                           int index,
                                           const GlobalPoint& gp,
                                           float z_pos_L3, float bestFitPhi, float bestFitDPhi)
{
  double csc_x = gp.x();
  double csc_y = gp.y();
  double csc_z = gp.z();
  double csc_R = TMath::Sqrt(gp.y()*gp.y() + gp.x()*gp.x());
  float radius = csc_R;

  // std::cout << "Printing stub properties" << std::endl 
  //           << "Id " << ch_id << " stub " << stub 
  //           << "GP " << gp << " eta " << gp.eta() << " phi " << gp.phi() << std::endl << std::endl;  

  switch(ch_id.station()) {
  case 1:
    event_.CSCTF_id1[index] = ch_id.rawId();
    event_.CSCTF_st1[index] = ch_id.station();
    event_.CSCTF_ri1[index] = ch_id.ring(); 
    event_.CSCTF_ch1[index] = ch_id.chamber();
    event_.CSCTF_en1[index] = ch_id.zendcap();
    event_.CSCTF_trk1[index] = stub.getTrknmb(); 
    event_.CSCTF_quality1[index] = stub.getQuality();
    event_.CSCTF_wg1[index] = stub.getKeyWG();
    event_.CSCTF_hs1[index] = stub.getStrip();
    event_.CSCTF_pat1[index] = stub.getPattern();
    event_.CSCTF_bend1[index] = stub.getBend();
    event_.CSCTF_bx1[index] = stub.getBX();
    event_.CSCTF_clctpat1[index] = stub.getCLCTPattern();
    event_.CSCTF_val1[index] = stub.isValid();
    event_.CSCTF_phi1[index] = gp.phi();
    event_.CSCTF_eta1[index] = gp.eta();
    //event_.CSCTF_gemdphi1[index] = stub.getGEMDPhi();
    event_.CSCTF_R1[index] = csc_R;
    event_.CSCTF_x1[index] = csc_x;
    event_.CSCTF_y1[index] = csc_y;
    event_.CSCTF_z1[index] = csc_z;
    // fitted positions
    event_.CSCTF_fit_phi1[index] = bestFitPhi;
    event_.CSCTF_fit_dphi1[index] = bestFitDPhi;
    event_.CSCTF_fit_x1[index] = radius*cos(bestFitPhi);
    event_.CSCTF_fit_y1[index] = radius*sin(bestFitPhi);
    event_.CSCTF_fit_z1[index] = z_pos_L3;
    event_.CSCTF_fit_R1[index] = radius;
    break;
  case 2:
    event_.CSCTF_id2[index] = ch_id.rawId();
    event_.CSCTF_st2[index] = ch_id.station();
    event_.CSCTF_ri2[index] = ch_id.ring(); 
    event_.CSCTF_ch2[index] = ch_id.chamber();
    event_.CSCTF_en2[index] = ch_id.zendcap();
    event_.CSCTF_trk2[index] = stub.getTrknmb(); 
    event_.CSCTF_quality2[index] = stub.getQuality();
    event_.CSCTF_wg2[index] = stub.getKeyWG();
    event_.CSCTF_hs2[index] = stub.getStrip();
    event_.CSCTF_pat2[index] = stub.getPattern();
    event_.CSCTF_bend2[index] = stub.getBend();
    event_.CSCTF_bx2[index] = stub.getBX();
    event_.CSCTF_clctpat2[index] = stub.getCLCTPattern();
    event_.CSCTF_val2[index] = stub.isValid();
    event_.CSCTF_phi2[index] = gp.phi();
    event_.CSCTF_eta2[index] = gp.eta();
    //event_.CSCTF_gemdphi2[index] = stub.getGEMDPhi();
    event_.CSCTF_R2[index] = csc_R;
    event_.CSCTF_x2[index] = csc_x;
    event_.CSCTF_y2[index] = csc_y;
    event_.CSCTF_z2[index] = csc_z;
    // fitted positions
    event_.CSCTF_fit_phi2[index] = bestFitPhi;
    event_.CSCTF_fit_dphi2[index] = bestFitDPhi;
    event_.CSCTF_fit_x2[index] = radius*cos(bestFitPhi);
    event_.CSCTF_fit_y2[index] = radius*sin(bestFitPhi);
    event_.CSCTF_fit_z2[index] = z_pos_L3;
    event_.CSCTF_fit_R2[index] = radius;
    break;
  case 3:
    event_.CSCTF_id3[index] = ch_id.rawId();
    event_.CSCTF_st3[index] = ch_id.station();
    event_.CSCTF_ri3[index] = ch_id.ring(); 
    event_.CSCTF_ch3[index] = ch_id.chamber();
    event_.CSCTF_en3[index] = ch_id.zendcap();
    event_.CSCTF_trk3[index] = stub.getTrknmb(); 
    event_.CSCTF_quality3[index] = stub.getQuality();
    event_.CSCTF_wg3[index] = stub.getKeyWG();
    event_.CSCTF_hs3[index] = stub.getStrip();
    event_.CSCTF_pat3[index] = stub.getPattern();
    event_.CSCTF_bend3[index] = stub.getBend();
    event_.CSCTF_bx3[index] = stub.getBX();
    event_.CSCTF_clctpat3[index] = stub.getCLCTPattern();
    event_.CSCTF_val3[index] = stub.isValid();
    event_.CSCTF_phi3[index] = gp.phi();
    event_.CSCTF_eta3[index] = gp.eta();
    event_.CSCTF_R3[index] = csc_R;
    event_.CSCTF_x3[index] = csc_x;
    event_.CSCTF_y3[index] = csc_y;
    event_.CSCTF_z3[index] = csc_z;
    // fitted positions
    event_.CSCTF_fit_phi3[index] = bestFitPhi;
    event_.CSCTF_fit_dphi3[index] = bestFitDPhi;
    event_.CSCTF_fit_x3[index] = radius*cos(bestFitPhi);
    event_.CSCTF_fit_y3[index] = radius*sin(bestFitPhi);
    event_.CSCTF_fit_z3[index] = z_pos_L3;
    event_.CSCTF_fit_R3[index] = radius;
    break;
  case 4:
    event_.CSCTF_id4[index] = ch_id.rawId();
    event_.CSCTF_st4[index] = ch_id.station();
    event_.CSCTF_ri4[index] = ch_id.ring(); 
    event_.CSCTF_ch4[index] = ch_id.chamber();
    event_.CSCTF_en4[index] = ch_id.zendcap();
    event_.CSCTF_trk4[index] = stub.getTrknmb(); 
    event_.CSCTF_quality4[index] = stub.getQuality();
    event_.CSCTF_wg4[index] = stub.getKeyWG();
    event_.CSCTF_hs4[index] = stub.getStrip();
    event_.CSCTF_pat4[index] = stub.getPattern();
    event_.CSCTF_bend4[index] = stub.getBend();
    event_.CSCTF_bx4[index] = stub.getBX();
    event_.CSCTF_clctpat4[index] = stub.getCLCTPattern();
    event_.CSCTF_val4[index] = stub.isValid();
    event_.CSCTF_phi4[index] = gp.phi();
    event_.CSCTF_eta4[index] = gp.eta();
    event_.CSCTF_R4[index] = csc_R;
    event_.CSCTF_x4[index] = csc_x;
    event_.CSCTF_y4[index] = csc_y;
    event_.CSCTF_z4[index] = csc_z;
    // fitted positions
    event_.CSCTF_fit_phi4[index] = bestFitPhi;
    event_.CSCTF_fit_dphi4[index] = bestFitDPhi;
    event_.CSCTF_fit_x4[index] = radius*cos(bestFitPhi);
    event_.CSCTF_fit_y4[index] = radius*sin(bestFitPhi);
    event_.CSCTF_fit_z4[index] = z_pos_L3;
    event_.CSCTF_fit_R4[index] = radius;
    break;
  };
}


void 
DisplacedL1MuFilter::printCSCStubProperties(int index) const
{
  
  cout << "id1 "<< event_.CSCTF_id1[index] << endl;
  cout << "st1 "<< event_.CSCTF_st1[index] << endl;
  cout << "ri1 "<< event_.CSCTF_ri1[index] << endl; 
  cout << "ch1 "<< event_.CSCTF_ch1[index] << endl;
  cout << "en1 "<< event_.CSCTF_en1[index] << endl;
  cout << "trk1 "<< event_.CSCTF_trk1[index]<< endl; 
  cout << "quality1 "<< event_.CSCTF_quality1[index] << endl;
  cout << "wg1 "<< event_.CSCTF_wg1[index] << endl;
  cout << "hs1 "<< event_.CSCTF_hs1[index] << endl;
  cout << "pat1 "<< event_.CSCTF_pat1[index] << endl;
  cout << "bend1 "<< event_.CSCTF_bend1[index]<< endl;
  cout << "bx1 "<< event_.CSCTF_bx1[index]<< endl;
  cout << "clctpat1 "<< event_.CSCTF_clctpat1[index]<< endl;
  cout << "val1 "<< event_.CSCTF_val1[index]<< endl;
  cout << "phi1 "<< event_.CSCTF_phi1[index]<< endl;
  cout << "eta1 "<< event_.CSCTF_eta1[index]<< endl;
  cout << "gemdphi1 "<< event_.CSCTF_gemdphi1[index]<< endl;
  cout << "R1 "<< event_.CSCTF_R1[index]<< endl;
  cout << "x1 "<< event_.CSCTF_x1[index]<< endl;
  cout << "y1 "<< event_.CSCTF_y1[index]<< endl;
  cout << "z1 "<< event_.CSCTF_z1[index]<< endl;
  cout << "fit phi1 "<< event_.CSCTF_fit_phi1[index]<< endl;
  cout << "fit dphi1 "<< event_.CSCTF_fit_dphi1[index]<< endl;
  cout << "fit x1 "<< event_.CSCTF_fit_x1[index]<< endl;
  cout << "fit y1 "<< event_.CSCTF_fit_y1[index]<< endl;
  cout << "fit z1 "<< event_.CSCTF_fit_z1[index]<< endl;
  cout << "fit R1 "<< event_.CSCTF_fit_R1[index]<< endl;

  cout << "id2 "<< event_.CSCTF_id2[index] << endl;
  cout << "st2 "<< event_.CSCTF_st2[index] << endl;
  cout << "ri2 "<< event_.CSCTF_ri2[index] << endl; 
  cout << "ch2 "<< event_.CSCTF_ch2[index] << endl;
  cout << "en2 "<< event_.CSCTF_en2[index] << endl;
  cout << "trk2 "<< event_.CSCTF_trk2[index]<< endl; 
  cout << "quality2 "<< event_.CSCTF_quality2[index] << endl;
  cout << "wg2 "<< event_.CSCTF_wg2[index] << endl;
  cout << "hs2 "<< event_.CSCTF_hs2[index] << endl;
  cout << "pat2 "<< event_.CSCTF_pat2[index] << endl;
  cout << "bend2 "<< event_.CSCTF_bend2[index]<< endl;
  cout << "bx2 "<< event_.CSCTF_bx2[index]<< endl;
  cout << "clctpat2 "<< event_.CSCTF_clctpat2[index]<< endl;
  cout << "val2 "<< event_.CSCTF_val2[index]<< endl;
  cout << "phi2 "<< event_.CSCTF_phi2[index]<< endl;
  cout << "eta2 "<< event_.CSCTF_eta2[index]<< endl;
  cout << "gemdphi2 "<< event_.CSCTF_gemdphi2[index]<< endl;
  cout << "R2 "<< event_.CSCTF_R2[index]<< endl;
  cout << "x2 "<< event_.CSCTF_x2[index]<< endl;
  cout << "y2 "<< event_.CSCTF_y2[index]<< endl;
  cout << "z2 "<< event_.CSCTF_z2[index]<< endl;
  cout << "fit phi2 "<< event_.CSCTF_fit_phi2[index]<< endl;
  cout << "fit dphi2 "<< event_.CSCTF_fit_dphi2[index]<< endl;
  cout << "fit x2 "<< event_.CSCTF_fit_x2[index]<< endl;
  cout << "fit y2 "<< event_.CSCTF_fit_y2[index]<< endl;
  cout << "fit z2 "<< event_.CSCTF_fit_z2[index]<< endl;
  cout << "fit R2 "<< event_.CSCTF_fit_R2[index]<< endl;

  cout << "id3 "<< event_.CSCTF_id3[index] << endl;
  cout << "st3 "<< event_.CSCTF_st3[index] << endl;
  cout << "ri3 "<< event_.CSCTF_ri3[index] << endl; 
  cout << "ch3 "<< event_.CSCTF_ch3[index] << endl;
  cout << "en3 "<< event_.CSCTF_en3[index] << endl;
  cout << "trk3 "<< event_.CSCTF_trk3[index]<< endl; 
  cout << "quality3 "<< event_.CSCTF_quality3[index] << endl;
  cout << "wg3 "<< event_.CSCTF_wg3[index] << endl;
  cout << "hs3 "<< event_.CSCTF_hs3[index] << endl;
  cout << "pat3 "<< event_.CSCTF_pat3[index] << endl;
  cout << "bend3 "<< event_.CSCTF_bend3[index]<< endl;
  cout << "bx3 "<< event_.CSCTF_bx3[index]<< endl;
  cout << "clctpat3 "<< event_.CSCTF_clctpat3[index]<< endl;
  cout << "val3 "<< event_.CSCTF_val3[index]<< endl;
  cout << "phi3 "<< event_.CSCTF_phi3[index]<< endl;
  cout << "eta3 "<< event_.CSCTF_eta3[index]<< endl;
  cout << "R3 "<< event_.CSCTF_R3[index]<< endl;
  cout << "x3 "<< event_.CSCTF_x3[index]<< endl;
  cout << "y3 "<< event_.CSCTF_y3[index]<< endl;
  cout << "z3 "<< event_.CSCTF_z3[index]<< endl;
  cout << "fit phi3 "<< event_.CSCTF_fit_phi3[index]<< endl;
  cout << "fit dphi3 "<< event_.CSCTF_fit_dphi3[index]<< endl;
  cout << "fit x3 "<< event_.CSCTF_fit_x3[index]<< endl;
  cout << "fit y3 "<< event_.CSCTF_fit_y3[index]<< endl;
  cout << "fit z3 "<< event_.CSCTF_fit_z3[index]<< endl;
  cout << "fit R3 "<< event_.CSCTF_fit_R3[index]<< endl;

  cout << "id4 "<< event_.CSCTF_id4[index] << endl;
  cout << "st4 "<< event_.CSCTF_st4[index] << endl;
  cout << "ri4 "<< event_.CSCTF_ri4[index] << endl; 
  cout << "ch4 "<< event_.CSCTF_ch4[index] << endl;
  cout << "en4 "<< event_.CSCTF_en4[index] << endl;
  cout << "trk4 "<< event_.CSCTF_trk4[index]<< endl; 
  cout << "quality4 "<< event_.CSCTF_quality4[index] << endl;
  cout << "wg4 "<< event_.CSCTF_wg4[index] << endl;
  cout << "hs4 "<< event_.CSCTF_hs4[index] << endl;
  cout << "pat4 "<< event_.CSCTF_pat4[index] << endl;
  cout << "bend4 "<< event_.CSCTF_bend4[index]<< endl;
  cout << "bx4 "<< event_.CSCTF_bx4[index]<< endl;
  cout << "clctpat4 "<< event_.CSCTF_clctpat4[index]<< endl;
  cout << "val4 "<< event_.CSCTF_val4[index]<< endl;
  cout << "phi4 "<< event_.CSCTF_phi4[index]<< endl;
  cout << "eta4 "<< event_.CSCTF_eta4[index]<< endl;
  cout << "R4 "<< event_.CSCTF_R4[index]<< endl;
  cout << "x4 "<< event_.CSCTF_x4[index]<< endl;
  cout << "y4 "<< event_.CSCTF_y4[index]<< endl;
  cout << "z4 "<< event_.CSCTF_z4[index]<< endl;
  cout << "fit phi4 "<< event_.CSCTF_fit_phi4[index]<< endl;
  cout << "fit dphi4 "<< event_.CSCTF_fit_dphi4[index]<< endl;
  cout << "fit x4 "<< event_.CSCTF_fit_x4[index]<< endl;
  cout << "fit y4 "<< event_.CSCTF_fit_y4[index]<< endl;
  cout << "fit z4 "<< event_.CSCTF_fit_z4[index]<< endl;
  cout << "fit R4 "<< event_.CSCTF_fit_R4[index]<< endl;

}


GlobalPoint 
DisplacedL1MuFilter::getCSCSpecificPoint(unsigned int rawId, const CSCCorrelatedLCTDigi& tp) const 
{
  const CSCDetId id(rawId); 
  // we should change this to weak_ptrs at some point
  // requires introducing std::shared_ptrs to geometry
  std::unique_ptr<const CSCChamber> chamb(cscGeometry_->chamber(id));
  std::unique_ptr<const CSCLayerGeometry> layer_geom(
                                                     chamb->layer(CSCConstants::KEY_ALCT_LAYER)->geometry()
                                                     );
  std::unique_ptr<const CSCLayer> layer(
                                        chamb->layer(CSCConstants::KEY_ALCT_LAYER)
                                        );
  
  const uint16_t halfstrip = tp.getStrip();
  const uint16_t pattern = tp.getPattern();
  const uint16_t keyWG = tp.getKeyWG(); 
  //const unsigned maxStrips = layer_geom->numberOfStrips();  

  // so we can extend this later 
  // assume TMB2007 half-strips only as baseline
  double offset = 0.0;
  switch(1) {
  case 1:
    offset = CSCPatternLUT::get2007Position(pattern);
  }
  const unsigned halfstrip_offs = unsigned(0.5 + halfstrip + offset);
  const unsigned strip = halfstrip_offs/2 + 1; // geom starts from 1

  // the rough location of the hit at the ALCT key layer
  // we will refine this using the half strip information
  const LocalPoint coarse_lp = 
    layer_geom->stripWireGroupIntersection(strip,keyWG);  
  const GlobalPoint coarse_gp = layer->surface().toGlobal(coarse_lp);  
  
  // the strip width/4.0 gives the offset of the half-strip
  // center with respect to the strip center
  const double hs_offset = layer_geom->stripPhiPitch()/4.0;
  
  // determine handedness of the chamber
  const bool ccw = isCSCCounterClockwise(layer);
  // we need to subtract the offset of even half strips and add the odd ones
  const double phi_offset = ( ( halfstrip_offs%2 ? 1 : -1)*
                              ( ccw ? -hs_offset : hs_offset ) );
  
  // the global eta calculation uses the middle of the strip
  // so no need to increment it
  const GlobalPoint final_gp( GlobalPoint::Polar( coarse_gp.theta(),
                                                  (coarse_gp.phi().value() + 
                                                   phi_offset),
                                                  coarse_gp.mag() ) );
    
  // We need to add in some notion of the 'error' on trigger primitives
  // like the width of the wire group by the width of the strip
  // or something similar      

  // release ownership of the pointers
  chamb.release();
  layer_geom.release();
  layer.release();
  
  return final_gp;
}

bool 
DisplacedL1MuFilter::isCSCCounterClockwise(const std::unique_ptr<const CSCLayer>& layer) const {
  const int nStrips = layer->geometry()->numberOfStrips();
  const double phi1 = layer->centerOfStrip(1).phi();
  const double phiN = layer->centerOfStrip(nStrips).phi();
  return ( (std::abs(phi1 - phiN) < M_PI  && phi1 >= phiN) || 
           (std::abs(phi1 - phiN) >= M_PI && phi1 < phiN)     );  
}


void DisplacedL1MuFilter::bookL1MuTree()
{
  edm::Service< TFileService > fs;
  event_tree_ = fs->make<TTree>("L1MuTree", "L1MuTree");
  event_tree_->Branch("lumi", &event_.lumi);
  event_tree_->Branch("run", &event_.run);
  event_tree_->Branch("event", &event_.event);

  // Beam spot
  event_tree_->Branch("beamSpot_x",    &event_.beamSpot_x,    "beamSpot_x/F");
  event_tree_->Branch("beamSpot_y",    &event_.beamSpot_y,    "beamSpot_y/F");
  event_tree_->Branch("beamSpot_z",    &event_.beamSpot_z,    "beamSpot_z/F");

  event_tree_->Branch("nL1Mu", &event_.nL1Mu);

  event_tree_->Branch("L1Mu_pt",event_.L1Mu_pt,"L1Mu_pt[nL1Mu]/F");
  event_tree_->Branch("L1Mu_eta", event_.L1Mu_eta,"L1Mu_eta[nL1Mu]/F");
  event_tree_->Branch("L1Mu_phi", event_.L1Mu_phi,"L1Mu_phi[nL1Mu]/F");
  event_tree_->Branch("L1Mu_bx", event_.L1Mu_bx,"L1Mu_bx[nL1Mu]/I");
  event_tree_->Branch("L1Mu_charge", event_.L1Mu_charge,"L1Mu_charge[nL1Mu]/I");
  event_tree_->Branch("L1Mu_quality", event_.L1Mu_quality,"L1Mu_quality[nL1Mu]/I");

  event_tree_->Branch("nDTTF", &event_.nDTTF);
  event_tree_->Branch("L1Mu_DTTF_index", event_.L1Mu_DTTF_index,"L1Mu_DTTF_index[nL1Mu]/I");

  event_tree_->Branch("DTTF_pt", event_.DTTF_pt,"DTTF_pt[nDTTF]/F");
  event_tree_->Branch("DTTF_eta", event_.DTTF_eta,"DTTF_eta[nDTTF]/F");
  event_tree_->Branch("DTTF_phi", event_.DTTF_phi,"DTTF_phi[nDTTF]/F");
  event_tree_->Branch("DTTF_bx", event_.DTTF_bx,"DTTF_bx[nDTTF]/F");
  event_tree_->Branch("DTTF_nStubs", event_.DTTF_nStubs,"DTTF_nStubs[nDTTF]/F");

  event_tree_->Branch("DTTF_phi1", event_.DTTF_phi1,"DTTF_phi1[nDTTF]/F");
  event_tree_->Branch("DTTF_phib1", event_.DTTF_phib1,"DTTF_phib1[nDTTF]/F");
  event_tree_->Branch("DTTF_quality1", event_.DTTF_quality1,"DTTF_quality1[nDTTF]/F");
  event_tree_->Branch("DTTF_bx1", event_.DTTF_bx1,"DTTF_bx1[nDTTF]/F");
  event_tree_->Branch("DTTF_wh1", event_.DTTF_wh1,"DTTF_wh1[nDTTF]/F");
  event_tree_->Branch("DTTF_se1", event_.DTTF_se1,"DTTF_se1[nDTTF]/F");
  event_tree_->Branch("DTTF_st1", event_.DTTF_st1,"DTTF_st1[nDTTF]/F");

  event_tree_->Branch("DTTF_phi2", event_.DTTF_phi2,"DTTF_phi2[nDTTF]/F");
  event_tree_->Branch("DTTF_phib2", event_.DTTF_phib2,"DTTF_phib2[nDTTF]/F");
  event_tree_->Branch("DTTF_quality2", event_.DTTF_quality2,"DTTF_quality2[nDTTF]/F");
  event_tree_->Branch("DTTF_bx2", event_.DTTF_bx2,"DTTF_bx2[nDTTF]/F");
  event_tree_->Branch("DTTF_wh2", event_.DTTF_wh2,"DTTF_wh2[nDTTF]/F");
  event_tree_->Branch("DTTF_se2", event_.DTTF_se2,"DTTF_se2[nDTTF]/F");
  event_tree_->Branch("DTTF_st2", event_.DTTF_st2,"DTTF_st2[nDTTF]/F");

  event_tree_->Branch("DTTF_phi3", event_.DTTF_phi3,"DTTF_phi3[nDTTF]/F");
  event_tree_->Branch("DTTF_phib3", event_.DTTF_phib3,"DTTF_phib3[nDTTF]/F");
  event_tree_->Branch("DTTF_quality3", event_.DTTF_quality3,"DTTF_quality3[nDTTF]/F");
  event_tree_->Branch("DTTF_bx3", event_.DTTF_bx3,"DTTF_bx3[nDTTF]/F");
  event_tree_->Branch("DTTF_wh3", event_.DTTF_wh3,"DTTF_wh3[nDTTF]/F");
  event_tree_->Branch("DTTF_se3", event_.DTTF_se3,"DTTF_se3[nDTTF]/F");
  event_tree_->Branch("DTTF_st3", event_.DTTF_st3,"DTTF_st3[nDTTF]/F");

  event_tree_->Branch("DTTF_phi4", event_.DTTF_phi4,"DTTF_phi4[nDTTF]/F");
  event_tree_->Branch("DTTF_phib4", event_.DTTF_phib4,"DTTF_phib4[nDTTF]/F");
  event_tree_->Branch("DTTF_quality4", event_.DTTF_quality4,"DTTF_quality4[nDTTF]/F");
  event_tree_->Branch("DTTF_bx4", event_.DTTF_bx4,"DTTF_bx4[nDTTF]/F");
  event_tree_->Branch("DTTF_wh4", event_.DTTF_wh4,"DTTF_wh4[nDTTF]/F");
  event_tree_->Branch("DTTF_se4", event_.DTTF_se4,"DTTF_se4[nDTTF]/F");
  event_tree_->Branch("DTTF_st4", event_.DTTF_st4,"DTTF_st4[nDTTF]/F");


  event_tree_->Branch("nCSCTF", &event_.nCSCTF);
  event_tree_->Branch("L1Mu_CSCTF_index", event_.L1Mu_CSCTF_index,"L1Mu_CSCTF_index[nL1Mu]/I");

  event_tree_->Branch("CSCTF_pt", event_.CSCTF_pt,"CSCTF_pt[nCSCTF]/F");
  event_tree_->Branch("CSCTF_eta", event_.CSCTF_eta,"CSCTF_eta[nCSCTF]/F");
  event_tree_->Branch("CSCTF_phi", event_.CSCTF_phi,"CSCTF_phi[nCSCTF]/F");
  event_tree_->Branch("CSCTF_bx", event_.CSCTF_bx,"CSCTF_bx[nCSCTF]/I");
  event_tree_->Branch("CSCTF_nStubs", event_.CSCTF_nStubs,"CSCTF_nStubs[nCSCTF]/I");
  event_tree_->Branch("CSCTF_quality", event_.CSCTF_quality,"CSCTF_quality[nCSCTF]/I");

  event_tree_->Branch("CSCTF_st1", event_.CSCTF_st1,"CSCTF_st1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ri1", event_.CSCTF_ri1,"CSCTF_ri1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ch1", event_.CSCTF_ch1,"CSCTF_ch1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_en1", event_.CSCTF_en1,"CSCTF_en1[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_trk1", event_.CSCTF_trk1,"CSCTF_trk1[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_quality1", event_.CSCTF_quality1,"CSCTF_quality1[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_wg1", event_.CSCTF_wg1,"CSCTF_wg1[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_hs1", event_.CSCTF_hs1,"CSCTF_hs1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_pat1", event_.CSCTF_pat1,"CSCTF_pat1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bend1", event_.CSCTF_bend1,"CSCTF_bend1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bx1", event_.CSCTF_bx1,"CSCTF_bx1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_clctpat1", event_.CSCTF_clctpat1,"CSCTF_clctpat1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_val1", event_.CSCTF_val1,"CSCTF_val1[nCSCTF]/I");
  event_tree_->Branch("CSCTF_phi1", event_.CSCTF_phi1,"CSCTF_phi1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_phib1", event_.CSCTF_phib1,"CSCTF_phib1[nCSCTF]/F");

  event_tree_->Branch("CSCTF_st2", event_.CSCTF_st2,"CSCTF_st2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ri2", event_.CSCTF_ri2,"CSCTF_ri2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ch2", event_.CSCTF_ch2,"CSCTF_ch2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_en2", event_.CSCTF_en2,"CSCTF_en2[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_trk2", event_.CSCTF_trk2,"CSCTF_trk2[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_quality2", event_.CSCTF_quality2,"CSCTF_quality2[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_wg2", event_.CSCTF_wg2,"CSCTF_wg2[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_hs2", event_.CSCTF_hs2,"CSCTF_hs2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_pat2", event_.CSCTF_pat2,"CSCTF_pat2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bend2", event_.CSCTF_bend2,"CSCTF_bend2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bx2", event_.CSCTF_bx2,"CSCTF_bx2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_clctpat2", event_.CSCTF_clctpat2,"CSCTF_clctpat2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_val2", event_.CSCTF_val2,"CSCTF_val2[nCSCTF]/I");
  event_tree_->Branch("CSCTF_phi2", event_.CSCTF_phi2,"CSCTF_phi2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_phib2", event_.CSCTF_phib2,"CSCTF_phib2[nCSCTF]/F");

  event_tree_->Branch("CSCTF_st3", event_.CSCTF_st3,"CSCTF_st3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ri3", event_.CSCTF_ri3,"CSCTF_ri3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ch3", event_.CSCTF_ch3,"CSCTF_ch3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_en3", event_.CSCTF_en3,"CSCTF_en3[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_trk3", event_.CSCTF_trk3,"CSCTF_trk3[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_quality3", event_.CSCTF_quality3,"CSCTF_quality3[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_wg3", event_.CSCTF_wg3,"CSCTF_wg3[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_hs3", event_.CSCTF_hs3,"CSCTF_hs3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_pat3", event_.CSCTF_pat3,"CSCTF_pat3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bend3", event_.CSCTF_bend3,"CSCTF_bend3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bx3", event_.CSCTF_bx3,"CSCTF_bx3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_clctpat3", event_.CSCTF_clctpat3,"CSCTF_clctpat3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_val3", event_.CSCTF_val3,"CSCTF_val3[nCSCTF]/I");
  event_tree_->Branch("CSCTF_phi3", event_.CSCTF_phi3,"CSCTF_phi3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_phib3", event_.CSCTF_phib3,"CSCTF_phib3[nCSCTF]/F");

  event_tree_->Branch("CSCTF_st4", event_.CSCTF_st4,"CSCTF_st4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ri4", event_.CSCTF_ri4,"CSCTF_ri4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_ch4", event_.CSCTF_ch4,"CSCTF_ch4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_en4", event_.CSCTF_en4,"CSCTF_en4[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_trk4", event_.CSCTF_trk4,"CSCTF_trk4[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_quality4", event_.CSCTF_quality4,"CSCTF_quality4[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_wg4", event_.CSCTF_wg4,"CSCTF_wg4[nCSCTF]/I");
  // event_tree_->Branch("CSCTF_hs4", event_.CSCTF_hs4,"CSCTF_hs4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_pat4", event_.CSCTF_pat4,"CSCTF_pat4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bend4", event_.CSCTF_bend4,"CSCTF_bend4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_bx4", event_.CSCTF_bx4,"CSCTF_bx4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_clctpat4", event_.CSCTF_clctpat4,"CSCTF_clctpat4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_val4", event_.CSCTF_val4,"CSCTF_val4[nCSCTF]/I");
  event_tree_->Branch("CSCTF_phi4", event_.CSCTF_phi4,"CSCTF_phi4[nCSCTF]/F");
  event_tree_->Branch("CSCTF_phib4", event_.CSCTF_phib4,"CSCTF_phib4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_eta1", event_.CSCTF_eta1,"CSCTF_eta1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_eta2", event_.CSCTF_eta2,"CSCTF_eta2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_eta3", event_.CSCTF_eta3,"CSCTF_eta3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_eta4", event_.CSCTF_eta4,"CSCTF_eta4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_gemdphi1", event_.CSCTF_gemdphi1,"CSCTF_gemdphi1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_gemdphi2", event_.CSCTF_gemdphi2,"CSCTF_gemdphi2[nCSCTF]/F");

  // true positions
  event_tree_->Branch("CSCTF_R1", event_.CSCTF_R1,"CSCTF_R1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_x1", event_.CSCTF_x1,"CSCTF_x1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_y1", event_.CSCTF_y1,"CSCTF_y1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_z1", event_.CSCTF_z1,"CSCTF_z1[nCSCTF]/F");

  event_tree_->Branch("CSCTF_R2", event_.CSCTF_R2,"CSCTF_R2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_x2", event_.CSCTF_x2,"CSCTF_x2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_y2", event_.CSCTF_y2,"CSCTF_y2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_z2", event_.CSCTF_z2,"CSCTF_z2[nCSCTF]/F");

  event_tree_->Branch("CSCTF_R3", event_.CSCTF_R3,"CSCTF_R3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_x3", event_.CSCTF_x3,"CSCTF_x3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_y3", event_.CSCTF_y3,"CSCTF_y3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_z3", event_.CSCTF_z3,"CSCTF_z3[nCSCTF]/F");

  event_tree_->Branch("CSCTF_R4", event_.CSCTF_R4,"CSCTF_R4[nCSCTF]/F");
  event_tree_->Branch("CSCTF_x4", event_.CSCTF_x4,"CSCTF_x4[nCSCTF]/F");
  event_tree_->Branch("CSCTF_y4", event_.CSCTF_y4,"CSCTF_y4[nCSCTF]/F");
  event_tree_->Branch("CSCTF_z4", event_.CSCTF_z4,"CSCTF_z4[nCSCTF]/F");


  event_tree_->Branch("CSCTF_fit_phi1", event_.CSCTF_fit_phi1,"CSCTF_fit_phi1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_phi2", event_.CSCTF_fit_phi2,"CSCTF_fit_phi2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_phi3", event_.CSCTF_fit_phi3,"CSCTF_fit_phi3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_phi4", event_.CSCTF_fit_phi4,"CSCTF_fit_phi4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fit_dphi1", event_.CSCTF_fit_dphi1,"CSCTF_fit_dphi1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_dphi2", event_.CSCTF_fit_dphi2,"CSCTF_fit_dphi2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_dphi3", event_.CSCTF_fit_dphi3,"CSCTF_fit_dphi3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_dphi4", event_.CSCTF_fit_dphi4,"CSCTF_fit_dphi4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fit_R1", event_.CSCTF_fit_R1,"CSCTF_fit_R1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_R2", event_.CSCTF_fit_R2,"CSCTF_fit_R2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_R3", event_.CSCTF_fit_R3,"CSCTF_fit_R3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_R4", event_.CSCTF_fit_R4,"CSCTF_fit_R4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fit_x1", event_.CSCTF_fit_x1,"CSCTF_fit_x1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_x2", event_.CSCTF_fit_x2,"CSCTF_fit_x2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_x3", event_.CSCTF_fit_x3,"CSCTF_fit_x3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_x4", event_.CSCTF_fit_x4,"CSCTF_fit_x4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fit_y1", event_.CSCTF_fit_y1,"CSCTF_fit_y1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_y2", event_.CSCTF_fit_y2,"CSCTF_fit_y2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_y3", event_.CSCTF_fit_y3,"CSCTF_fit_y3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_y4", event_.CSCTF_fit_y4,"CSCTF_fit_y4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fit_z1", event_.CSCTF_fit_z1,"CSCTF_fit_z1[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_z2", event_.CSCTF_fit_z2,"CSCTF_fit_z2[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_z3", event_.CSCTF_fit_z3,"CSCTF_fit_z3[nCSCTF]/F");
  event_tree_->Branch("CSCTF_fit_z4", event_.CSCTF_fit_z4,"CSCTF_fit_z4[nCSCTF]/F");

  event_tree_->Branch("CSCTF_fitline_x1", event_.CSCTF_fitline_x1,"CSCTF_fitline_x1[4]/F");
  event_tree_->Branch("CSCTF_fitline_x2", event_.CSCTF_fitline_x2,"CSCTF_fitline_x2[4]/F");
  event_tree_->Branch("CSCTF_fitline_x3", event_.CSCTF_fitline_x3,"CSCTF_fitline_x3[4]/F");
  event_tree_->Branch("CSCTF_fitline_x4", event_.CSCTF_fitline_x4,"CSCTF_fitline_x4[4]/F");

  event_tree_->Branch("CSCTF_fitline_y1", event_.CSCTF_fitline_y1,"CSCTF_fitline_y1[4]/F");
  event_tree_->Branch("CSCTF_fitline_y2", event_.CSCTF_fitline_y2,"CSCTF_fitline_y2[4]/F");
  event_tree_->Branch("CSCTF_fitline_y3", event_.CSCTF_fitline_y3,"CSCTF_fitline_y3[4]/F");
  event_tree_->Branch("CSCTF_fitline_y4", event_.CSCTF_fitline_y4,"CSCTF_fitline_y4[4]/F");


  if (processRPCb_) {
  event_tree_->Branch("nRPCb", &event_.nRPCb);
  event_tree_->Branch("L1Mu_RPCb_index", event_.L1Mu_RPCb_index,"L1Mu_RPCb_index[nL1Mu]/I");

  event_tree_->Branch("RPCb_pt", event_.RPCb_pt,"RPCb_pt[nRPCb]/F");
  event_tree_->Branch("RPCb_eta", event_.RPCb_eta,"RPCb_eta[nRPCb]/F");
  event_tree_->Branch("RPCb_phi", event_.RPCb_phi,"RPCb_phi[nRPCb]/F");
  event_tree_->Branch("RPCb_bx", event_.RPCb_bx,"RPCb_bx[nRPCb]/I");
  event_tree_->Branch("RPCb_nStubs", event_.RPCb_nStubs,"RPCb_nStubs[nRPCb]/I");
  event_tree_->Branch("RPCb_quality", event_.RPCb_quality,"RPCb_quality[nRPCb]/I");

  event_tree_->Branch("RPCb_bx1", event_.RPCb_bx1,"RPCb_bx1[nRPCb]/I");
  event_tree_->Branch("RPCb_strip1", event_.RPCb_strip1,"RPCb_strip1[nRPCb]/I");
  event_tree_->Branch("RPCb_phi1", event_.RPCb_phi1,"RPCb_phi1[nRPCb]/F");
  event_tree_->Branch("RPCb_re1", event_.RPCb_re1,"RPCb_re1[nRPCb]/I");
  event_tree_->Branch("RPCb_ri1", event_.RPCb_ri1,"RPCb_ri1[nRPCb]/I");
  event_tree_->Branch("RPCb_st1", event_.RPCb_st1,"RPCb_st1[nRPCb]/I");
  event_tree_->Branch("RPCb_se1", event_.RPCb_se1,"RPCb_se1[nRPCb]/I");
  event_tree_->Branch("RPCb_la1", event_.RPCb_la1,"RPCb_la1[nRPCb]/I");
  event_tree_->Branch("RPCb_su1", event_.RPCb_su1,"RPCb_su1[nRPCb]/I");
  event_tree_->Branch("RPCb_ro1", event_.RPCb_ro1,"RPCb_ro1[nRPCb]/I");
  
  event_tree_->Branch("RPCb_bx2", event_.RPCb_bx2,"RPCb_bx2[nRPCb]/I");
  event_tree_->Branch("RPCb_strip2", event_.RPCb_strip2,"RPCb_strip2[nRPCb]/I");
  event_tree_->Branch("RPCb_phi2", event_.RPCb_phi2,"RPCb_phi2[nRPCb]/F");
  event_tree_->Branch("RPCb_re2", event_.RPCb_re2,"RPCb_re2[nRPCb]/I");
  event_tree_->Branch("RPCb_ri2", event_.RPCb_ri2,"RPCb_ri2[nRPCb]/I");
  event_tree_->Branch("RPCb_st2", event_.RPCb_st2,"RPCb_st2[nRPCb]/I");
  event_tree_->Branch("RPCb_se2", event_.RPCb_se2,"RPCb_se2[nRPCb]/I");
  event_tree_->Branch("RPCb_la2", event_.RPCb_la2,"RPCb_la2[nRPCb]/I");
  event_tree_->Branch("RPCb_su2", event_.RPCb_su2,"RPCb_su2[nRPCb]/I");
  event_tree_->Branch("RPCb_ro2", event_.RPCb_ro2,"RPCb_ro2[nRPCb]/I");

  event_tree_->Branch("RPCb_bx3", event_.RPCb_bx3,"RPCb_bx3[nRPCb]/I");
  event_tree_->Branch("RPCb_strip3", event_.RPCb_strip3,"RPCb_strip3[nRPCb]/I");
  event_tree_->Branch("RPCb_phi3", event_.RPCb_phi3,"RPCb_phi3[nRPCb]/F");
  event_tree_->Branch("RPCb_re3", event_.RPCb_re3,"RPCb_re3[nRPCb]/I");
  event_tree_->Branch("RPCb_ri3", event_.RPCb_ri3,"RPCb_ri3[nRPCb]/I");
  event_tree_->Branch("RPCb_st3", event_.RPCb_st3,"RPCb_st3[nRPCb]/I");
  event_tree_->Branch("RPCb_se3", event_.RPCb_se3,"RPCb_se3[nRPCb]/I");
  event_tree_->Branch("RPCb_la3", event_.RPCb_la3,"RPCb_la3[nRPCb]/I");
  event_tree_->Branch("RPCb_su3", event_.RPCb_su3,"RPCb_su3[nRPCb]/I");
  event_tree_->Branch("RPCb_ro3", event_.RPCb_ro3,"RPCb_ro3[nRPCb]/I");

  event_tree_->Branch("RPCb_bx4", event_.RPCb_bx4,"RPCb_bx4[nRPCb]/I");
  event_tree_->Branch("RPCb_strip4", event_.RPCb_strip4,"RPCb_strip4[nRPCb]/I");
  event_tree_->Branch("RPCb_phi4", event_.RPCb_phi4,"RPCb_phi4[nRPCb]/F");
  event_tree_->Branch("RPCb_re4", event_.RPCb_re4,"RPCb_re4[nRPCb]/I");
  event_tree_->Branch("RPCb_ri4", event_.RPCb_ri4,"RPCb_ri4[nRPCb]/I");
  event_tree_->Branch("RPCb_st4", event_.RPCb_st4,"RPCb_st4[nRPCb]/I");
  event_tree_->Branch("RPCb_se4", event_.RPCb_se4,"RPCb_se4[nRPCb]/I");
  event_tree_->Branch("RPCb_la4", event_.RPCb_la4,"RPCb_la4[nRPCb]/I");
  event_tree_->Branch("RPCb_su4", event_.RPCb_su4,"RPCb_su4[nRPCb]/I");
  event_tree_->Branch("RPCb_ro4", event_.RPCb_ro4,"RPCb_ro4[nRPCb]/I");

  event_tree_->Branch("RPCb_bx5", event_.RPCb_bx5,"RPCb_bx5[nRPCb]/I");
  event_tree_->Branch("RPCb_strip5", event_.RPCb_strip5,"RPCb_strip5[nRPCb]/I");
  event_tree_->Branch("RPCb_phi5", event_.RPCb_phi5,"RPCb_phi5[nRPCb]/F");
  event_tree_->Branch("RPCb_re5", event_.RPCb_re5,"RPCb_re5[nRPCb]/I");
  event_tree_->Branch("RPCb_ri5", event_.RPCb_ri5,"RPCb_ri5[nRPCb]/I");
  event_tree_->Branch("RPCb_st5", event_.RPCb_st5,"RPCb_st5[nRPCb]/I");
  event_tree_->Branch("RPCb_se5", event_.RPCb_se5,"RPCb_se5[nRPCb]/I");
  event_tree_->Branch("RPCb_la5", event_.RPCb_la5,"RPCb_la5[nRPCb]/I");
  event_tree_->Branch("RPCb_su5", event_.RPCb_su5,"RPCb_su5[nRPCb]/I");
  event_tree_->Branch("RPCb_ro5", event_.RPCb_ro5,"RPCb_ro5[nRPCb]/I");

  event_tree_->Branch("RPCb_bx6", event_.RPCb_bx6,"RPCb_bx6[nRPCb]/I");
  event_tree_->Branch("RPCb_strip6", event_.RPCb_strip6,"RPCb_strip6[nRPCb]/I");
  event_tree_->Branch("RPCb_phi6", event_.RPCb_phi6,"RPCb_phi6[nRPCb]/F");
  event_tree_->Branch("RPCb_re6", event_.RPCb_re6,"RPCb_re6[nRPCb]/I");
  event_tree_->Branch("RPCb_ri6", event_.RPCb_ri6,"RPCb_ri6[nRPCb]/I");
  event_tree_->Branch("RPCb_st6", event_.RPCb_st6,"RPCb_st6[nRPCb]/I");
  event_tree_->Branch("RPCb_se6", event_.RPCb_se6,"RPCb_se6[nRPCb]/I");
  event_tree_->Branch("RPCb_la6", event_.RPCb_la6,"RPCb_la6[nRPCb]/I");
  event_tree_->Branch("RPCb_su6", event_.RPCb_su6,"RPCb_su6[nRPCb]/I");
  event_tree_->Branch("RPCb_ro6", event_.RPCb_ro6,"RPCb_ro6[nRPCb]/I");
  }

  if (processRPCf_) {
  event_tree_->Branch("nRPCf", &event_.nRPCf);
  event_tree_->Branch("L1Mu_RPCf_index", event_.L1Mu_RPCf_index,"L1Mu_RPCf_index[nL1Mu]/I");

  event_tree_->Branch("RPCf_pt", event_.RPCf_pt,"RPCf_pt[nRPCf]/F");
  event_tree_->Branch("RPCf_eta", event_.RPCf_eta,"RPCf_eta[nRPCf]/F");
  event_tree_->Branch("RPCf_phi", event_.RPCf_phi,"RPCf_phi[nRPCf]/F");
  event_tree_->Branch("RPCf_bx", event_.RPCf_bx,"RPCf_bx[nRPCf]/I");
  event_tree_->Branch("RPCf_nStubs", event_.RPCf_nStubs,"RPCf_nStubs[nRPCf]/I");
  event_tree_->Branch("RPCf_quality", event_.RPCf_quality,"RPCf_quality[nRPCf]/I");

  event_tree_->Branch("RPCf_bx1", event_.RPCf_bx1,"RPCf_bx1[nRPCf]/I");
  event_tree_->Branch("RPCf_strip1", event_.RPCf_strip1,"RPCf_strip1[nRPCf]/I");
  event_tree_->Branch("RPCf_phi1", event_.RPCf_phi1,"RPCf_phi1[nRPCf]/F");
  event_tree_->Branch("RPCf_re1", event_.RPCf_re1,"RPCf_re1[nRPCf]/I");
  event_tree_->Branch("RPCf_ri1", event_.RPCf_ri1,"RPCf_ri1[nRPCf]/I");
  event_tree_->Branch("RPCf_st1", event_.RPCf_st1,"RPCf_st1[nRPCf]/I");
  event_tree_->Branch("RPCf_se1", event_.RPCf_se1,"RPCf_se1[nRPCf]/I");
  event_tree_->Branch("RPCf_la1", event_.RPCf_la1,"RPCf_la1[nRPCf]/I");
  event_tree_->Branch("RPCf_su1", event_.RPCf_su1,"RPCf_su1[nRPCf]/I");
  event_tree_->Branch("RPCf_ro1", event_.RPCf_ro1,"RPCf_ro1[nRPCf]/I");
  
  event_tree_->Branch("RPCf_bx2", event_.RPCf_bx2,"RPCf_bx2[nRPCf]/I");
  event_tree_->Branch("RPCf_strip2", event_.RPCf_strip2,"RPCf_strip2[nRPCf]/I");
  event_tree_->Branch("RPCf_phi2", event_.RPCf_phi2,"RPCf_phi2[nRPCf]/F");
  event_tree_->Branch("RPCf_re2", event_.RPCf_re2,"RPCf_re2[nRPCf]/I");
  event_tree_->Branch("RPCf_ri2", event_.RPCf_ri2,"RPCf_ri2[nRPCf]/I");
  event_tree_->Branch("RPCf_st2", event_.RPCf_st2,"RPCf_st2[nRPCf]/I");
  event_tree_->Branch("RPCf_se2", event_.RPCf_se2,"RPCf_se2[nRPCf]/I");
  event_tree_->Branch("RPCf_la2", event_.RPCf_la2,"RPCf_la2[nRPCf]/I");
  event_tree_->Branch("RPCf_su2", event_.RPCf_su2,"RPCf_su2[nRPCf]/I");
  event_tree_->Branch("RPCf_ro2", event_.RPCf_ro2,"RPCf_ro2[nRPCf]/I");

  event_tree_->Branch("RPCf_bx3", event_.RPCf_bx3,"RPCf_bx3[nRPCf]/I");
  event_tree_->Branch("RPCf_strip3", event_.RPCf_strip3,"RPCf_strip3[nRPCf]/I");
  event_tree_->Branch("RPCf_phi3", event_.RPCf_phi3,"RPCf_phi3[nRPCf]/F");
  event_tree_->Branch("RPCf_re3", event_.RPCf_re3,"RPCf_re3[nRPCf]/I");
  event_tree_->Branch("RPCf_ri3", event_.RPCf_ri3,"RPCf_ri3[nRPCf]/I");
  event_tree_->Branch("RPCf_st3", event_.RPCf_st3,"RPCf_st3[nRPCf]/I");
  event_tree_->Branch("RPCf_se3", event_.RPCf_se3,"RPCf_se3[nRPCf]/I");
  event_tree_->Branch("RPCf_la3", event_.RPCf_la3,"RPCf_la3[nRPCf]/I");
  event_tree_->Branch("RPCf_su3", event_.RPCf_su3,"RPCf_su3[nRPCf]/I");
  event_tree_->Branch("RPCf_ro3", event_.RPCf_ro3,"RPCf_ro3[nRPCf]/I");

  event_tree_->Branch("RPCf_bx4", event_.RPCf_bx4,"RPCf_bx4[nRPCf]/I");
  event_tree_->Branch("RPCf_strip4", event_.RPCf_strip4,"RPCf_strip4[nRPCf]/I");
  event_tree_->Branch("RPCf_phi4", event_.RPCf_phi4,"RPCf_phi4[nRPCf]/F");
  event_tree_->Branch("RPCf_re4", event_.RPCf_re4,"RPCf_re4[nRPCf]/I");
  event_tree_->Branch("RPCf_ri4", event_.RPCf_ri4,"RPCf_ri4[nRPCf]/I");
  event_tree_->Branch("RPCf_st4", event_.RPCf_st4,"RPCf_st4[nRPCf]/I");
  event_tree_->Branch("RPCf_se4", event_.RPCf_se4,"RPCf_se4[nRPCf]/I");
  event_tree_->Branch("RPCf_la4", event_.RPCf_la4,"RPCf_la4[nRPCf]/I");
  event_tree_->Branch("RPCf_su4", event_.RPCf_su4,"RPCf_su4[nRPCf]/I");
  event_tree_->Branch("RPCf_ro4", event_.RPCf_ro4,"RPCf_ro4[nRPCf]/I");

  event_tree_->Branch("RPCf_bx5", event_.RPCf_bx5,"RPCf_bx5[nRPCf]/I");
  event_tree_->Branch("RPCf_strip5", event_.RPCf_strip5,"RPCf_strip5[nRPCf]/I");
  event_tree_->Branch("RPCf_phi5", event_.RPCf_phi5,"RPCf_phi5[nRPCf]/F");
  event_tree_->Branch("RPCf_re5", event_.RPCf_re5,"RPCf_re5[nRPCf]/I");
  event_tree_->Branch("RPCf_ri5", event_.RPCf_ri5,"RPCf_ri5[nRPCf]/I");
  event_tree_->Branch("RPCf_st5", event_.RPCf_st5,"RPCf_st5[nRPCf]/I");
  event_tree_->Branch("RPCf_se5", event_.RPCf_se5,"RPCf_se5[nRPCf]/I");
  event_tree_->Branch("RPCf_la5", event_.RPCf_la5,"RPCf_la5[nRPCf]/I");
  event_tree_->Branch("RPCf_su5", event_.RPCf_su5,"RPCf_su5[nRPCf]/I");
  event_tree_->Branch("RPCf_ro5", event_.RPCf_ro5,"RPCf_ro5[nRPCf]/I");

  event_tree_->Branch("RPCf_bx6", event_.RPCf_bx6,"RPCf_bx6[nRPCf]/I");
  event_tree_->Branch("RPCf_strip6", event_.RPCf_strip6,"RPCf_strip6[nRPCf]/I");
  event_tree_->Branch("RPCf_phi6", event_.RPCf_phi6,"RPCf_phi6[nRPCf]/F");
  event_tree_->Branch("RPCf_re6", event_.RPCf_re6,"RPCf_re6[nRPCf]/I");
  event_tree_->Branch("RPCf_ri6", event_.RPCf_ri6,"RPCf_ri6[nRPCf]/I");
  event_tree_->Branch("RPCf_st6", event_.RPCf_st6,"RPCf_st6[nRPCf]/I");
  event_tree_->Branch("RPCf_se6", event_.RPCf_se6,"RPCf_se6[nRPCf]/I");
  event_tree_->Branch("RPCf_la6", event_.RPCf_la6,"RPCf_la6[nRPCf]/I");
  event_tree_->Branch("RPCf_su6", event_.RPCf_su6,"RPCf_su6[nRPCf]/I");
  event_tree_->Branch("RPCf_ro6", event_.RPCf_ro6,"RPCf_ro6[nRPCf]/I");
  }

  event_tree_->Branch("nGEM", &event_.nGEM);
  event_tree_->Branch("GE11_phi_L1", event_.GE11_phi_L1,"GE11_phi_L1[4]/F");
  event_tree_->Branch("GE11_phi_L2", event_.GE11_phi_L2,"GE11_phi_L2[4]/F");
  event_tree_->Branch("GE21_phi_L1", event_.GE21_phi_L1,"GE21_phi_L1[4]/F");
  event_tree_->Branch("GE21_phi_L2", event_.GE21_phi_L2,"GE21_phi_L2[4]/F");
  event_tree_->Branch("GE11_bx_L1", event_.GE11_bx_L1,"GE11_bx_L1[4]/I");
  event_tree_->Branch("GE11_bx_L2", event_.GE11_bx_L2,"GE11_bx_L2[4]/I");
  event_tree_->Branch("GE21_bx_L1", event_.GE21_bx_L1,"GE21_bx_L1[4]/I");
  event_tree_->Branch("GE21_bx_L2", event_.GE21_bx_L2,"GE21_bx_L2[4]/I");
  event_tree_->Branch("GE11_ch_L1", event_.GE11_ch_L1,"GE11_ch_L1[4]/I");
  event_tree_->Branch("GE11_ch_L2", event_.GE11_ch_L2,"GE11_ch_L2[4]/I");
  event_tree_->Branch("GE21_ch_L1", event_.GE21_ch_L1,"GE21_ch_L1[4]/I");
  event_tree_->Branch("GE21_ch_L2", event_.GE21_ch_L2,"GE21_ch_L2[4]/I");
  event_tree_->Branch("GE11_z_L1", event_.GE11_z_L1,"GE11_z_L1[4]/F");
  event_tree_->Branch("GE11_z_L2", event_.GE11_z_L2,"GE11_z_L2[4]/F");
  event_tree_->Branch("GE21_z_L1", event_.GE21_z_L1,"GE21_z_L1[4]/F");
  event_tree_->Branch("GE21_z_L2", event_.GE21_z_L2,"GE21_z_L2[4]/F");

  event_tree_->Branch("GE0_phi", event_.GE0_phi,"GE0_phi[4]/F");
  event_tree_->Branch("GE0_phib", event_.GE0_phib,"GE0_phib[4]/F");

  event_tree_->Branch("GE21_pad1_phi_L1", event_.GE21_pad1_phi_L1,"GE21_pad1_phi_L1[4]/F");
  event_tree_->Branch("GE21_pad1_phi_L2", event_.GE21_pad1_phi_L2,"GE21_pad1_phi_L2[4]/F");
  event_tree_->Branch("GE21_pad2_phi_L1", event_.GE21_pad2_phi_L1,"GE21_pad2_phi_L1[4]/F");
  event_tree_->Branch("GE21_pad2_phi_L2", event_.GE21_pad2_phi_L2,"GE21_pad2_phi_L2[4]/F");
  event_tree_->Branch("GE21_pad4_phi_L1", event_.GE21_pad4_phi_L1,"GE21_pad4_phi_L1[4]/F");
  event_tree_->Branch("GE21_pad4_phi_L2", event_.GE21_pad4_phi_L2,"GE21_pad4_phi_L2[4]/F");
  event_tree_->Branch("GE21_pad8_phi_L1", event_.GE21_pad8_phi_L1,"GE21_pad8_phi_L1[4]/F");
  event_tree_->Branch("GE21_pad8_phi_L2", event_.GE21_pad8_phi_L2,"GE21_pad8_phi_L2[4]/F");
}



void 
DisplacedL1MuFilter::clearBranches()
{
  event_.lumi = -99;
  event_.run = -99;
  event_.event = -99;

  event_.nL1Mu = 0;
  event_.nL1Tk = 0;


  for (int i=0; i<kMaxL1Mu; ++i){
    event_.L1Mu_pt[i] = -99;
    event_.L1Mu_eta[i] = -99;
    event_.L1Mu_phi[i] = -99;
    event_.L1Mu_charge[i] = -99;
    event_.L1Mu_bx[i] = -99;
    event_.L1Mu_quality[i] = -99;
    event_.L1Mu_DTTF_index[i] = -1;
    event_.L1Mu_CSCTF_index[i] = -1;
    event_.L1Mu_RPCb_index[i] = -1;
    event_.L1Mu_RPCf_index[i] = -1;
  }

  event_.nDTTF = 0;

  for (int i=0; i<kMaxDTTF; ++i){
    event_.DTTF_pt[i] = 99;
    event_.DTTF_eta[i] = 99;
    event_.DTTF_phi[i] = 99;
    event_.DTTF_nStubs[i] = 0;

    event_.DTTF_phi1[i] = 99;
    event_.DTTF_phib1[i] = 99;
    event_.DTTF_quality1[i] = 99;
    event_.DTTF_bx1[i] = 99;
    event_.DTTF_wh1[i] = 99;
    event_.DTTF_se1[i] = 99;
    event_.DTTF_st1[i] = 99;

    event_.DTTF_phi2[i] = 99;
    event_.DTTF_phib2[i] = 99;
    event_.DTTF_quality2[i] = 99;
    event_.DTTF_bx2[i] = 99;
    event_.DTTF_wh2[i] = 99;
    event_.DTTF_se2[i] = 99;
    event_.DTTF_st2[i] = 99;

    event_.DTTF_phi3[i] = 99;
    event_.DTTF_phib3[i] = 99;
    event_.DTTF_quality3[i] = 99;
    event_.DTTF_bx3[i] = 99;
    event_.DTTF_wh3[i] = 99;
    event_.DTTF_se3[i] = 99;
    event_.DTTF_st3[i] = 99;

    event_.DTTF_phi4[i] = 99;
    event_.DTTF_phib4[i] = 99;
    event_.DTTF_quality4[i] = 99;
    event_.DTTF_bx4[i] = 99;
    event_.DTTF_wh4[i] = 99;
    event_.DTTF_se4[i] = 99;
    event_.DTTF_st4[i] = 99;
  }


  event_.nCSCTF = 0;

  for (int i=0; i<kMaxCSCTF; ++i){
    event_.CSCTF_pt[i] = 99;
    event_.CSCTF_eta[i] = 99;
    event_.CSCTF_phi[i] = 99;
    event_.CSCTF_nStubs[i] = 0;

    event_.CSCTF_st1[i] = 99; 
    event_.CSCTF_ri1[i] = 99; 
    event_.CSCTF_ch1[i] = 99; 
    event_.CSCTF_en1[i] = 99;
    event_.CSCTF_trk1[i] = 99; 
    event_.CSCTF_quality1[i] = 99; 
    event_.CSCTF_wg1[i] = 99; 
    event_.CSCTF_hs1[i] = 99; 
    event_.CSCTF_pat1[i] = 99; 
    event_.CSCTF_bend1[i] = 99; 
    event_.CSCTF_bx1[i] = 99; 
    event_.CSCTF_clctpat1[i] = 99;
    event_.CSCTF_val1[i] = 99;
    event_.CSCTF_phi1[i] = 99;
    event_.CSCTF_phib1[i] = 99;

    event_.CSCTF_st2[i] = 99; 
    event_.CSCTF_ri2[i] = 99; 
    event_.CSCTF_ch2[i] = 99; 
    event_.CSCTF_en2[i] = 99;
    event_.CSCTF_trk2[i] = 99; 
    event_.CSCTF_quality2[i] = 99; 
    event_.CSCTF_wg2[i] = 99; 
    event_.CSCTF_hs2[i] = 99; 
    event_.CSCTF_pat2[i] = 99; 
    event_.CSCTF_bend2[i] = 99; 
    event_.CSCTF_bx2[i] = 99; 
    event_.CSCTF_clctpat2[i] = 99;
    event_.CSCTF_val2[i] = 99;
    event_.CSCTF_phi2[i] = 99;
    event_.CSCTF_phib2[i] = 99;

    event_.CSCTF_st3[i] = 99; 
    event_.CSCTF_ri3[i] = 99; 
    event_.CSCTF_ch3[i] = 99; 
    event_.CSCTF_en3[i] = 99;
    event_.CSCTF_trk3[i] = 99; 
    event_.CSCTF_quality3[i] = 99; 
    event_.CSCTF_wg3[i] = 99; 
    event_.CSCTF_hs3[i] = 99; 
    event_.CSCTF_pat3[i] = 99; 
    event_.CSCTF_bend3[i] = 99; 
    event_.CSCTF_bx3[i] = 99; 
    event_.CSCTF_clctpat3[i] = 99;
    event_.CSCTF_val3[i] = 99;
    event_.CSCTF_phi3[i] = 99;
    event_.CSCTF_phib3[i] = 99;

    event_.CSCTF_st4[i] = 99; 
    event_.CSCTF_ri4[i] = 99; 
    event_.CSCTF_ch4[i] = 99; 
    event_.CSCTF_en4[i] = 99;
    event_.CSCTF_trk4[i] = 99; 
    event_.CSCTF_quality4[i] = 99; 
    event_.CSCTF_wg4[i] = 99; 
    event_.CSCTF_hs4[i] = 99; 
    event_.CSCTF_pat4[i] = 99; 
    event_.CSCTF_bend4[i] = 99; 
    event_.CSCTF_bx4[i] = 99; 
    event_.CSCTF_clctpat4[i] = 99;
    event_.CSCTF_val4[i] = 99;
    event_.CSCTF_phi4[i] = 99;
    event_.CSCTF_phib4[i] = 99;

    event_.CSCTF_eta1[i] = 99;
    event_.CSCTF_eta2[i] = 99;
    event_.CSCTF_eta3[i] = 99;
    event_.CSCTF_eta4[i] = 99;

    event_.CSCTF_gemdphi1[i] = 99;
    event_.CSCTF_gemdphi2[i] = 99;

    event_.CSCTF_R1[i] = 99;
    event_.CSCTF_x1[i] = 99;
    event_.CSCTF_y1[i] = 99;
    event_.CSCTF_z1[i] = 99;
    event_.CSCTF_R2[i] = 99;
    event_.CSCTF_x2[i] = 99;
    event_.CSCTF_y2[i] = 99;
    event_.CSCTF_z2[i] = 99;
    event_.CSCTF_R3[i] = 99;
    event_.CSCTF_x3[i] = 99;
    event_.CSCTF_y3[i] = 99;
    event_.CSCTF_z3[i] = 99;
    event_.CSCTF_R4[i] = 99;
    event_.CSCTF_x4[i] = 99;
    event_.CSCTF_y4[i] = 99;
    event_.CSCTF_z4[i] = 99;

    event_.CSCTF_fit_R1[i] = 99;
    event_.CSCTF_fit_x1[i] = 99;
    event_.CSCTF_fit_y1[i] = 99;
    event_.CSCTF_fit_z1[i] = 99;
    event_.CSCTF_fit_R2[i] = 99;
    event_.CSCTF_fit_x2[i] = 99;
    event_.CSCTF_fit_y2[i] = 99;
    event_.CSCTF_fit_z2[i] = 99;
    event_.CSCTF_fit_R3[i] = 99;
    event_.CSCTF_fit_x3[i] = 99;
    event_.CSCTF_fit_y3[i] = 99;
    event_.CSCTF_fit_z3[i] = 99;
    event_.CSCTF_fit_R4[i] = 99;
    event_.CSCTF_fit_x4[i] = 99;
    event_.CSCTF_fit_y4[i] = 99;
    event_.CSCTF_fit_z4[i] = 99;


    event_.CSCTF_fit_phi1[i] = 99;
    event_.CSCTF_fit_phi2[i] = 99;
    event_.CSCTF_fit_phi3[i] = 99;
    event_.CSCTF_fit_phi4[i] = 99;

    event_.CSCTF_fit_dphi1[i] = 99;
    event_.CSCTF_fit_dphi2[i] = 99;
    event_.CSCTF_fit_dphi3[i] = 99;
    event_.CSCTF_fit_dphi4[i] = 99;

  }

  event_.nRPCb = 0;

  for (int i=0; i<kMaxRPCb; ++i){
    event_.RPCb_pt[i] = 99;
    event_.RPCb_eta[i] = 99;
    event_.RPCb_phi[i] = 99;
    event_.RPCb_bx[i] = 99;
    event_.RPCb_nStubs[i] = 0;
    event_.RPCb_quality[i] = 99;
    
    event_.RPCb_bx1[i] = 99;
    event_.RPCb_strip1[i] = 99;
    event_.RPCb_phi1[i] = 99;
    event_.RPCb_re1[i] = 99;
    event_.RPCb_ri1[i] = 99;
    event_.RPCb_st1[i] = 99;
    event_.RPCb_se1[i] = 99;
    event_.RPCb_la1[i] = 99;
    event_.RPCb_su1[i] = 99;
    event_.RPCb_ro1[i] = 99;

    event_.RPCb_bx2[i] = 99;
    event_.RPCb_strip2[i] = 99;
    event_.RPCb_phi2[i] = 99;
    event_.RPCb_re2[i] = 99;
    event_.RPCb_ri2[i] = 99;
    event_.RPCb_st2[i] = 99;
    event_.RPCb_se2[i] = 99;
    event_.RPCb_la2[i] = 99;
    event_.RPCb_su2[i] = 99;
    event_.RPCb_ro2[i] = 99;

    event_.RPCb_bx3[i] = 99;
    event_.RPCb_strip3[i] = 99;
    event_.RPCb_phi3[i] = 99;
    event_.RPCb_re3[i] = 99;
    event_.RPCb_ri3[i] = 99;
    event_.RPCb_st3[i] = 99;
    event_.RPCb_se3[i] = 99;
    event_.RPCb_la3[i] = 99;
    event_.RPCb_su3[i] = 99;
    event_.RPCb_ro3[i] = 99;

    event_.RPCb_bx4[i] = 99;
    event_.RPCb_strip4[i] = 99;
    event_.RPCb_phi4[i] = 99;
    event_.RPCb_re4[i] = 99;
    event_.RPCb_ri4[i] = 99;
    event_.RPCb_st4[i] = 99;
    event_.RPCb_se4[i] = 99;
    event_.RPCb_la4[i] = 99;
    event_.RPCb_su4[i] = 99;
    event_.RPCb_ro4[i] = 99;

    event_.RPCb_bx5[i] = 99;
    event_.RPCb_strip5[i] = 99;
    event_.RPCb_phi5[i] = 99;
    event_.RPCb_re5[i] = 99;
    event_.RPCb_ri5[i] = 99;
    event_.RPCb_st5[i] = 99;
    event_.RPCb_se5[i] = 99;
    event_.RPCb_la5[i] = 99;
    event_.RPCb_su5[i] = 99;
    event_.RPCb_ro5[i] = 99;

    event_.RPCb_bx6[i] = 99;
    event_.RPCb_strip6[i] = 99;
    event_.RPCb_phi6[i] = 99;
    event_.RPCb_re6[i] = 99;
    event_.RPCb_ri6[i] = 99;
    event_.RPCb_st6[i] = 99;
    event_.RPCb_se6[i] = 99;
    event_.RPCb_la6[i] = 99;
    event_.RPCb_su6[i] = 99;
    event_.RPCb_ro6[i] = 99;
  }

  event_.nRPCf = 0;

  for (int i=0; i<kMaxRPCf; ++i){
    event_.RPCf_pt[i] = 99;
    event_.RPCf_eta[i] = 99;
    event_.RPCf_phi[i] = 99;
    event_.RPCf_bx[i] = 99;
    event_.RPCf_nStubs[i] = 0;
    event_.RPCf_quality[i] = 99;
    
    event_.RPCf_bx1[i] = 99;
    event_.RPCf_strip1[i] = 99;
    event_.RPCf_phi1[i] = 99;
    event_.RPCf_re1[i] = 99;
    event_.RPCf_ri1[i] = 99;
    event_.RPCf_st1[i] = 99;
    event_.RPCf_se1[i] = 99;
    event_.RPCf_la1[i] = 99;
    event_.RPCf_su1[i] = 99;
    event_.RPCf_ro1[i] = 99;

    event_.RPCf_bx2[i] = 99;
    event_.RPCf_strip2[i] = 99;
    event_.RPCf_phi2[i] = 99;
    event_.RPCf_re2[i] = 99;
    event_.RPCf_ri2[i] = 99;
    event_.RPCf_st2[i] = 99;
    event_.RPCf_se2[i] = 99;
    event_.RPCf_la2[i] = 99;
    event_.RPCf_su2[i] = 99;
    event_.RPCf_ro2[i] = 99;

    event_.RPCf_bx3[i] = 99;
    event_.RPCf_strip3[i] = 99;
    event_.RPCf_phi3[i] = 99;
    event_.RPCf_re3[i] = 99;
    event_.RPCf_ri3[i] = 99;
    event_.RPCf_st3[i] = 99;
    event_.RPCf_se3[i] = 99;
    event_.RPCf_la3[i] = 99;
    event_.RPCf_su3[i] = 99;
    event_.RPCf_ro3[i] = 99;

    event_.RPCf_bx4[i] = 99;
    event_.RPCf_strip4[i] = 99;
    event_.RPCf_phi4[i] = 99;
    event_.RPCf_re4[i] = 99;
    event_.RPCf_ri4[i] = 99;
    event_.RPCf_st4[i] = 99;
    event_.RPCf_se4[i] = 99;
    event_.RPCf_la4[i] = 99;
    event_.RPCf_su4[i] = 99;
    event_.RPCf_ro4[i] = 99;

    event_.RPCf_bx5[i] = 99;
    event_.RPCf_strip5[i] = 99;
    event_.RPCf_phi5[i] = 99;
    event_.RPCf_re5[i] = 99;
    event_.RPCf_ri5[i] = 99;
    event_.RPCf_st5[i] = 99;
    event_.RPCf_se5[i] = 99;
    event_.RPCf_la5[i] = 99;
    event_.RPCf_su5[i] = 99;
    event_.RPCf_ro5[i] = 99;

    event_.RPCf_bx6[i] = 99;
    event_.RPCf_strip6[i] = 99;
    event_.RPCf_phi6[i] = 99;
    event_.RPCf_re6[i] = 99;
    event_.RPCf_ri6[i] = 99;
    event_.RPCf_st6[i] = 99;
    event_.RPCf_se6[i] = 99;
    event_.RPCf_la6[i] = 99;
    event_.RPCf_su6[i] = 99;
    event_.RPCf_ro6[i] = 99;
  }


  event_.nGEM = 4;
  for (int i=0; i<4; ++i){
    event_.GE11_phi_L1[i] = 99.;
    event_.GE11_phi_L2[i] = 99.;
    event_.GE21_phi_L1[i] = 99.;
    event_.GE21_phi_L2[i] = 99.;
    event_.GE11_bx_L1[i] = 99;
    event_.GE11_bx_L2[i] = 99;
    event_.GE21_bx_L1[i] = 99;
    event_.GE21_bx_L2[i] = 99;
    event_.GE11_ch_L1[i] = 99;
    event_.GE11_ch_L2[i] = 99;
    event_.GE21_ch_L1[i] = 99;
    event_.GE21_ch_L2[i] = 99;
    event_.GE11_z_L1[i] = 99;
    event_.GE11_z_L2[i] = 99;
    event_.GE21_z_L1[i] = 99;
    event_.GE21_z_L2[i] = 99;
    event_.GE0_phi[i] = 99;
    event_.GE0_phib[i] = 99;

    event_.GE21_pad1_phi_L1[i] = 99.; 
    event_.GE21_pad1_phi_L2[i] = 99.;
    event_.GE21_pad2_phi_L1[i] = 99.; 
    event_.GE21_pad2_phi_L2[i] = 99.;
    event_.GE21_pad4_phi_L1[i] = 99.; 
    event_.GE21_pad4_phi_L2[i] = 99.;
    event_.GE21_pad8_phi_L1[i] = 99.; 
    event_.GE21_pad8_phi_L2[i] = 99.;

    event_.CSCTF_fitline_x1[i] = 99;
    event_.CSCTF_fitline_x2[i] = 99;
    event_.CSCTF_fitline_x3[i] = 99;
    event_.CSCTF_fitline_x4[i] = 99;
    
    event_.CSCTF_fitline_y1[i] = 99;
    event_.CSCTF_fitline_y2[i] = 99;
    event_.CSCTF_fitline_y3[i] = 99;
    event_.CSCTF_fitline_y4[i] = 99;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedL1MuFilter);

//  LocalWords:  kMaxCSCTF
