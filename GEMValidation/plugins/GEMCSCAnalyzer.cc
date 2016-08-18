/**\class GEMCSCAnalyzer

 Description:

 Analyzer of correlations of signals in CSC & GEM using SimTracks
 Needed for the GEM-CSC triggering algorithm development.

 Original Author:  "Vadim Khotilovich"
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"
#include "GEMCode/GEMValidation/interface/Helpers.h"
#include "GEMCode/GEMValidation/interface/Ptassignment.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/Randomize.h"

#include "TTree.h"

#include <iomanip>
#include <sstream>
#include <memory>
#include <math.h>
#include <bitset>
#include <vector>

using namespace std;
using namespace matching;


// "signed" LCT bend pattern
const int LCT_BEND_PATTERN[11] = { -99,  -5,  4, -4,  3, -3,  2, -2,  1, -1,  0};


struct MyTrackChamberDelta
{
  Bool_t odd;
  Char_t charge;
  Char_t chamber;
  Char_t endcap;
  Char_t roll;
  Char_t bend;
  Float_t pt, eta, phi;
  Float_t csc_sh_phi;
  Float_t csc_dg_phi;
  Float_t gem_sh_phi;
  Float_t gem_dg_phi;
  Float_t gem_pad_phi;
  Float_t dphi_sh;
  Float_t dphi_dg;
  Float_t dphi_pad;
  Float_t csc_sh_eta;
  Float_t csc_dg_eta;
  Float_t gem_sh_eta;
  Float_t gem_dg_eta;
  Float_t gem_pad_eta;
  Float_t deta_sh;
  Float_t deta_dg;
  Float_t deta_pad;
  Float_t csc_lct_phi;
  Float_t dphi_lct_pad;
  Float_t csc_lct_eta;
  Float_t deta_lct_pad;
  Float_t dphi_gem_sh_csc_sh;
  Float_t dphi_gem_dg_csc_dg;
  Float_t dphi_gem_pad_csc_lct;
  Float_t dphi_gem_rh_csc_seg;
};


struct MyTrackEff
{
  void init(); // initialize to default values
  TTree* book(TTree *t, const std::string & name = "trk_eff");

  Int_t lumi;
  Int_t run;
  Int_t event;

  Float_t pt, eta, phi, pz;
  Char_t charge;
  Char_t endcap;

  // Dark photon  
  Float_t genGd_m;
  Float_t genGd_E;
  Float_t genGd_p;
  Float_t genGd_pt;
  Float_t genGd_px;
  Float_t genGd_py;
  Float_t genGd_pz;
  Float_t genGd_eta;
  Float_t genGd_phi;
  Float_t genGd_vx;
  Float_t genGd_vy;
  Float_t genGd_vz;
  Float_t genGd_vLx;
  Float_t genGd_vLy;
  Float_t genGd_vLz;
  Float_t genGd_lxy;
  Float_t genGd_l;
  Float_t genGd_dxy;
  Float_t genGd0Gd1_dR;
  Float_t genGd0Gd1_m;

  Float_t genGdMu_dxy_max;
  Float_t genGdMu_eta_max;

  // Gen level muon
  Int_t genGd_index;
  Int_t genGdMu_index;
  Float_t genGdMu_p;
  Float_t genGdMu_pt;
  Float_t genGdMu_px;
  Float_t genGdMu_py;
  Float_t genGdMu_pz;
  Float_t genGdMu_eta;
  Float_t genGdMu_phi;
  Float_t genGdMu_phi_corr;
  Float_t genGdMu_vx;
  Float_t genGdMu_vy;
  Float_t genGdMu_vz;
  Float_t genGdMu_dxy;
  Float_t genGdMu_dR;

  Float_t beamSpot_x;
  Float_t beamSpot_y;
  Float_t beamSpot_z;


  Char_t chamber_ME1_csc_sh;//bit1:odd, bit2:even
  Char_t chamber_ME2_csc_sh;
  Char_t chamber_odd; // bit1: has GEM pad   bit2: has CSC LCT
  Char_t chamber_even; // bit1: has GEM pad   bit2: has CSC LCT
  Float_t bending_sh;
  Float_t phi_cscsh_even, phi_cscsh_odd, eta_cscsh_even, eta_cscsh_odd;
  Float_t phi_layer1_sh_even, eta_layer1_sh_even, phi_layer1_sh_odd,eta_layer1_sh_odd, perp_layer1_sh_odd, perp_layer1_sh_even;
  Float_t phi_layer6_sh_even, eta_layer6_sh_even, phi_layer6_sh_odd,eta_layer6_sh_odd, perp_layer6_sh_odd, perp_layer6_sh_even;
  Float_t perp_gemsh_even, perp_gemsh_odd, perp_cscsh_even, perp_cscsh_odd,centralperp_gemsh_even, centralperp_gemsh_odd;
  Float_t dphi_sh_even,dphi_sh_odd;
  Float_t dphipositionpt_cscsh_even, dphipositionpt_cscsh_odd, dphipositionpt_gemsh_even,dphipositionpt_gemsh_odd;

  Char_t has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
  Char_t has_csc_strips; // #layers with comparator digis > minHitsChamber    bit1: in odd, bit2: even
  Char_t has_csc_wires; // #layers with wire digis > minHitsChamber    bit1: in odd, bit2: even
  Float_t phi_layer1_fit_even, phi_layer3_fit_even, phi_layer6_fit_even, phi_layer1_fit_odd, phi_layer3_fit_odd, phi_layer6_fit_odd;

  Char_t has_clct; // bit1: in odd, bit2: even
  Char_t has_alct; // bit1: in odd, bit2: even
  Char_t has_lct; // bit1: in odd, bit2: even

  Char_t bend_lct_odd;
  Char_t bend_lct_even;
  Char_t bx_lct_odd;
  Char_t bx_lct_even;


  UChar_t hs_lct_odd;
  UChar_t wg_lct_odd;
  UChar_t hs_lct_even;
  UChar_t wg_lct_even;

  Float_t phi_lct_odd;
  Float_t phi_lct_even;
  Float_t eta_lct_odd;
  Float_t eta_lct_even;
  Float_t dphi_lct_odd; // dphi stored as data member in LCT
  Float_t dphi_lct_even;
  Float_t perp_lct_odd;
  Float_t perp_lct_even;
  Bool_t passdphi_odd;
  Bool_t passdphi_even;

  Int_t wiregroup_odd;
  Int_t wiregroup_even;
  Int_t halfstrip_odd;
  Int_t halfstrip_even;

  Int_t quality_clct_odd;
  Int_t quality_clct_even;
  Int_t quality_alct_odd;
  Int_t quality_alct_even;

  Int_t nlayers_csc_sh_odd;
  Int_t nlayers_wg_dg_odd;
  Int_t nlayers_st_dg_odd;
  Int_t nlayers_csc_sh_even;
  Int_t nlayers_wg_dg_even;
  Int_t nlayers_st_dg_even;
  Int_t pad_odd;
  Int_t pad_even;
  Int_t Copad_odd;
  Int_t Copad_even;
  Int_t hsfromgem_odd;
  Int_t hsfromgem_even;

  Char_t has_gem_sh; // bit1: in odd, bit2: even
  Char_t has_gem_sh2; // has SimHits in 2 layers  bit1: in odd, bit2: even
  Char_t has_gem_dg; // bit1: in odd, bit2: even
  Char_t has_gem_dg2; // has pads in 2 layers  bit1: in odd, bit2: even
  Char_t has_gem_pad; // bit1: in odd, bit2: even
  Char_t has_gem_pad2; // has pads in 2 layers  bit1: in odd, bit2: even
  Char_t has_gem_copad; // bit1: in odd, bit2: even

  Float_t strip_gemsh_odd; // average hits' strip
  Float_t strip_gemsh_even;
  Float_t eta_gemsh_odd;
  Float_t eta_gemsh_even;
  Float_t phi_gemsh_odd;
  Float_t phi_gemsh_even;
  Int_t strip_gemdg_odd; // median digis' strip
  Int_t strip_gemdg_even;
   

  Char_t has_rpc_sh; // bit1: in odd, bit2: even
  Char_t has_rpc_dg; // bit1: in odd, bit2: even
  Int_t strip_rpcdg_odd; // median digis' strip
  Int_t strip_rpcdg_even;

  Char_t bx_pad_odd;
  Char_t bx_pad_even;
  Float_t phi_pad_odd;
  Float_t phi_pad_even;
  Float_t eta_pad_odd;
  Float_t eta_pad_even;

  Float_t dphi_pad_odd;
  Float_t dphi_pad_even;
  Float_t deta_pad_odd;
  Float_t deta_pad_even;

  Int_t quality_odd;
  Int_t quality_even;

  Int_t hsfromrpc_odd; // extraplotate hs from rpc
  Int_t hsfromrpc_even;

  Char_t bx_rpcstrip_odd;
  Char_t bx_rpcstrip_even;
  Float_t phi_rpcstrip_odd;
  Float_t phi_rpcstrip_even;
  Float_t eta_rpcstrip_odd;
  Float_t eta_rpcstrip_even;

  Float_t dphi_rpcstrip_odd;
  Float_t dphi_rpcstrip_even;
  Float_t deta_rpcstrip_odd;
  Float_t deta_rpcstrip_even;

  // Track properties
  Int_t has_tfTrack;
  Int_t has_tfCand;
  Int_t has_gmtRegCand;
  Int_t has_gmtCand;
 
  //csctf
  Float_t trackpt, tracketa, trackphi;
  UInt_t quality_packed, pt_packed, eta_packed, phi_packed;
  UInt_t chargesign;
  UInt_t rank;
  UInt_t nstubs;
  UInt_t deltaphi12, deltaphi23; 
  Bool_t hasME1,hasME2;
  Char_t chamberME1,chamberME2;//bit1: odd, bit2: even
  Int_t ME1_ring, ME2_ring;
  Int_t ME1_hs, ME2_hs, ME1_wg,ME2_wg;
  Float_t dphiGE11,dphiGE21;
  Bool_t passGE11,passGE21;
  Bool_t passGE11_pt5, passGE11_pt7, passGE11_pt10, passGE11_pt15, passGE11_pt20, passGE11_pt30, passGE11_pt40;
  Bool_t passGE21_pt5, passGE21_pt7, passGE21_pt10, passGE21_pt15, passGE21_pt20, passGE21_pt30, passGE21_pt40;
  Bool_t passGE11_simpt, passGE21_simpt;
  Float_t deltaR;
  Float_t lctdphi12;
  Float_t eta_propagated_ME1;
  Float_t eta_propagated_ME2;
  Float_t eta_propagated_ME3;
  Float_t eta_propagated_ME4;
  Float_t phi_propagated_ME1;
  Float_t phi_propagated_ME2;
  Float_t phi_propagated_ME3;
  Float_t phi_propagated_ME4;
  Float_t eta_ME1_TF;
  Float_t eta_ME2_TF;
  Float_t eta_ME3_TF;
  Float_t eta_ME4_TF;
  Float_t phi_ME1_TF;
  Float_t phi_ME2_TF;
  Float_t phi_ME3_TF;
  Float_t phi_ME4_TF;

  Float_t eta_interStat12;
  Float_t phi_interStat12;
  Float_t eta_interStat23;
  Float_t phi_interStat23;
  Float_t eta_interStat13;
  Float_t phi_interStat13;

  Bool_t allstubs_matched_TF;

  Int_t has_l1Extra;
  Float_t l1Extra_pt;
  Float_t l1Extra_eta;
  Float_t l1Extra_phi;
  Float_t l1Extra_dR;
  Int_t has_recoTrackExtra;
  Float_t recoTrackExtra_pt_inner;
  Float_t recoTrackExtra_eta_inner;
  Float_t recoTrackExtra_phi_inner;  
  Float_t recoTrackExtra_pt_outer;
  Float_t recoTrackExtra_eta_outer;
  Float_t recoTrackExtra_phi_outer;
  Int_t has_recoTrack;
  Float_t recoTrack_pt_outer;
  Float_t recoTrack_eta_outer;
  Float_t recoTrack_phi_outer;
  Int_t has_recoChargedCandidate;
  Float_t recoChargedCandidate_pt;
  Float_t recoChargedCandidate_eta;
  Float_t recoChargedCandidate_phi;

  Int_t recoChargedCandidate_nValidCSCHits;
  Int_t recoChargedCandidate_nValidRPCHits;
  Int_t recoChargedCandidate_nValidDTHits;


  // pt assginment
  Float_t phiM_gemcsc_even,phiM_gemcsc_odd,phiM_gemcsc_central_even, phiM_gemcsc_central_odd;
  Float_t phiM_fitxyz_even,phiM_fitxyz_odd;
  Float_t pt_sh_even,pt_sh_odd,ptphi_sh_even,ptphi_sh_odd,pteta_sh_even,pteta_sh_odd;
  Float_t ptphi_gemsh_odd, ptphi_gemsh_even;
  Float_t csc_bending_angle12_gemcsc,csc_bending_angle12_gemcsc_central;
  Float_t csc_bending_angle12_xfactor,csc_bending_angle12_xfactor_smear1,csc_bending_angle12_xfactor_smear2,csc_bending_angle12_xfactor_smear0;
  Float_t csc_bending_angle12_xfactor_smear3; 
  Float_t csc_bending_angle12_xfactor_L1_1, csc_bending_angle12_xfactor_L1_2;
  //Float_t ptphi_diff_sh_11,ptphi_diff_sh_12,ptphi_diff_sh_21,ptphi_diff_sh_22; 
  Float_t ptphi_diff_sh;
  Float_t eta_st1_sh, eta_st2_sh;
  Float_t pt_position_sh, pt_direction_sh;
  Float_t pt_position, pt_direction_gemcsc, pt_direction_gemcsc_central, pt_direction_xfactor;
  Float_t pt_direction_xfactor_smear0, pt_direction_xfactor_smear1,pt_direction_xfactor_smear2,pt_direction_xfactor_smear3;
  Float_t pt_direction_xfactor_L1_1,pt_direction_xfactor_L1_2;
  Float_t pt_position_smeared, pt_direction_gemcsc_smeared;
  Int_t npar;
  Int_t npar_jose;
  Bool_t hasSt1St2St3;
  Bool_t hasSt3orSt4_sh;
};

void MyTrackEff::init()
{
  lumi = -99;
  run = -99;
  event = -99;

  pt = 0.;
  phi = 0.;
  eta = -9.;
  charge = -9;
  endcap = -9;

  // Dark photon  
  genGd_m =0.0;
  genGd_E =0.0;
  genGd_p =0.0;
  genGd_pt =0.0;
  genGd_px =0.0;
  genGd_py =0.0;
  genGd_pz =0.0;
  genGd_eta =0.0;
  genGd_phi =0.0;
  genGd_vx =0.0;
  genGd_vy =0.0;
  genGd_vz =0.0;
  genGd_vLx = -0.;
  genGd_vLy = -0;
  genGd_vLz= -0;
  genGd_lxy = -0.;
  genGd_l =-0.;
  genGd_dxy =-1;
  genGd0Gd1_dR =-1;
  genGd0Gd1_m = 0;

  genGdMu_dxy_max =0.0;
  genGdMu_eta_max =0.0;

  // Gen level muon
  genGd_index = -99;
  genGdMu_index = -99;
  genGdMu_p =0.0;
  genGdMu_pt =0.0;
  genGdMu_px =0.0;
  genGdMu_py =0.0;
  genGdMu_pz =0.0;
  genGdMu_eta =0.0;
  genGdMu_phi =0.0;
  genGdMu_phi_corr =0.0;
  genGdMu_vx =0.0;
  genGdMu_vy =0.0;
  genGdMu_vz =0.0;
  genGdMu_dxy =0.0;
  genGdMu_dR = 10.0;

  beamSpot_x =0.0;
  beamSpot_y =0.0;
  beamSpot_z =0.0;


  chamber_ME1_csc_sh=0;
  chamber_ME2_csc_sh=0;
  chamber_odd = 0;
  chamber_even = 0;
  quality_odd = 0;
  quality_even = 0;
  bending_sh = -10;
  perp_cscsh_odd = -0.0;
  perp_cscsh_even = -0.0;
  perp_gemsh_odd = -0.0;
  perp_gemsh_even = -0.0;
  centralperp_gemsh_odd = -0.0;
  centralperp_gemsh_even = -0.0;
  phi_cscsh_even = -9.0;
  phi_cscsh_odd = -9.0;
  eta_cscsh_even = -9.0;
  eta_cscsh_odd = -9.0;
  phi_layer1_sh_even =-9.0;
  phi_layer1_sh_odd =-9.0;
  eta_layer1_sh_even =-9.0;
  eta_layer1_sh_odd =-9.0;
  perp_layer1_sh_even = -1;
  perp_layer1_sh_odd = -1;
  phi_layer6_sh_even =-9.0;
  phi_layer6_sh_odd =-9.0;
  eta_layer6_sh_even =-9.0;
  eta_layer6_sh_odd =-9.0;
  perp_layer6_sh_even = -1;
  perp_layer6_sh_odd = -1;
  pt_sh_even = -9.0;
  pteta_sh_even = 0;
  ptphi_sh_even = -9.0;
  pt_sh_odd = -9.0;
  pteta_sh_odd = 0;
  ptphi_sh_odd = -9.0;
  ptphi_gemsh_even = -9.0;
  ptphi_gemsh_odd = -9.0;

  eta_st1_sh = -9;
  eta_st2_sh = -9;
  //ptphi_diff_sh_11 = -9;
  //ptphi_diff_sh_12 = -9;
  //ptphi_diff_sh_21 = -9;
  //ptphi_diff_sh_22 = -9;
  ptphi_diff_sh = -9;
  csc_bending_angle12_gemcsc = -9;
  csc_bending_angle12_xfactor = -9;
  csc_bending_angle12_xfactor_smear0 = -9;
  csc_bending_angle12_xfactor_smear1 = -9;
  csc_bending_angle12_xfactor_smear2 = -9;
  csc_bending_angle12_xfactor_smear3 = -9;
  csc_bending_angle12_xfactor_L1_1 = -9;
  csc_bending_angle12_xfactor_L1_2 = -9;

  phi_layer1_fit_even =-9.0;
  phi_layer1_fit_odd =-9.0;
  phi_layer3_fit_even =-9.0;
  phi_layer3_fit_odd =-9.0;
  phi_layer6_fit_even =-9.0;
  phi_layer6_fit_odd =-9.0;
  has_csc_sh = 0;
  has_csc_strips = 0;
  has_csc_wires = 0;
  has_alct = 0;
  has_clct = 0;
  has_lct = 0;
  bend_lct_odd = -9;
  bend_lct_even = -9;
  bx_lct_odd = -9;
  bx_lct_even = -9;
  hs_lct_odd = 0;
  hs_lct_even = 0;
  wg_lct_odd = 0;
  wg_lct_even = 0;
  phi_lct_odd = -9.;
  phi_lct_even = -9.;
  eta_lct_odd = -9.;
  eta_lct_even = -9.;
  dphi_lct_odd = -9.;
  dphi_lct_even = -9.;
  passdphi_odd = false;
  passdphi_even = false;

  wiregroup_odd = -1;
  wiregroup_even =-1; 
  halfstrip_odd =-1;
  halfstrip_even = -1;
  quality_clct_odd = -1;
  quality_clct_even = -1;
  quality_alct_odd = -1;
  quality_alct_even = -1;
  nlayers_csc_sh_odd = -1;
  nlayers_wg_dg_odd = -1;
  nlayers_st_dg_odd = -1;
  nlayers_csc_sh_even = -1;
  nlayers_wg_dg_even = -1;
  nlayers_st_dg_even = -1;
  pad_odd = -1;
  pad_even = -1;
  Copad_odd = -1;
  Copad_even = -1;

  hsfromgem_odd = -1;
  hsfromgem_even = -1;

  has_gem_sh = 0;
  has_gem_sh2 = 0;
  has_gem_dg = 0;
  has_gem_dg2 = 0;
  has_gem_pad = 0;
  has_gem_pad2 = 0;
  has_gem_copad = 0;
  strip_gemsh_odd = -9.;
  strip_gemsh_even = -9.;
  eta_gemsh_odd = -9.;
  eta_gemsh_even = -9.;
  phi_gemsh_odd = -9.;
  phi_gemsh_even = -9.;
  dphi_sh_odd = -9;
  dphi_sh_even = -9;
  dphipositionpt_gemsh_even = -9;
  dphipositionpt_gemsh_odd = -9;
  dphipositionpt_cscsh_even = -9;
  dphipositionpt_cscsh_odd = -9;
  phiM_gemcsc_odd = -9;
  phiM_gemcsc_even =-9;
  phiM_gemcsc_central_odd = -9;
  phiM_gemcsc_central_even =-9;
  phiM_fitxyz_odd = -9;
  phiM_fitxyz_even =-9;
  strip_gemdg_odd = -9;
  strip_gemdg_even = -9;
 
  has_rpc_sh = 0;
  has_rpc_dg = 0; // bit1: in odd, bit2: even
  strip_rpcdg_odd = -1;
  strip_rpcdg_even = -1;

  hsfromrpc_odd = 0;
  hsfromrpc_even = 0;

  bx_pad_odd = -9;
  bx_pad_even = -9;
  phi_pad_odd = -9.;
  phi_pad_even = -9.;
  eta_pad_odd = -9.;
  eta_pad_even = -9.;
  dphi_pad_odd = -9.;
  dphi_pad_even = -9.;
  deta_pad_odd = -9.;
  deta_pad_even = -9.;

  bx_rpcstrip_odd = -9;
  bx_rpcstrip_even = -9;
  phi_rpcstrip_odd = -9.;
  phi_rpcstrip_even = -9.;
  eta_rpcstrip_odd = -9.;
  eta_rpcstrip_even = -9.;
  dphi_rpcstrip_odd = -9.;
  dphi_rpcstrip_even = -9.;
  deta_rpcstrip_odd = -9.;
  deta_rpcstrip_even = -9.;

  // Track properties
  has_tfTrack = -99;
  has_tfCand = -99;
  has_gmtRegCand = -99;
  has_gmtCand = -99;

  //csctf
  trackpt = 0 ;
  tracketa = 0;
  trackphi = -9;
  quality_packed = 0;
  pt_packed = 0;
  eta_packed = 0;
  phi_packed = 0;
  ME1_hs = -1;
  ME1_wg = -1;
  ME2_hs = -1;
  ME2_wg = -1;
  chargesign =99;
  rank = 0;
  deltaphi12 = 0;
  deltaphi23 = 0;; 
  hasME1 = false;
  hasME2 = false;
  ME1_ring = -1;
  ME2_ring = -1;
  chamberME1 = 0;
  chamberME2 = 0;
  dphiGE11 = -99.0;
  dphiGE21 = -99.0;
  passGE11 = false;
  passGE11_pt5 = false;
  passGE11_pt7 = false;
  passGE11_pt10 = false;
  passGE11_pt15 = false;
  passGE11_pt20 = false;
  passGE11_pt30 = false;
  passGE11_pt40 = false;
  passGE21 = false;
  passGE21_pt5 = false;
  passGE21_pt7 = false;
  passGE21_pt10 = false;
  passGE21_pt15 = false;
  passGE21_pt20 = false;
  passGE21_pt30 = false;
  passGE21_pt40 = false;
  passGE11_simpt = false;
  passGE21_simpt = false;//to debug dphi cut eff
  nstubs = 0;
  deltaR = 10;
  lctdphi12 = -99;

  eta_propagated_ME1 = -9;
  eta_propagated_ME2 = -9;
  eta_propagated_ME3 = -9;
  eta_propagated_ME4 = -9;
  phi_propagated_ME1 = -9;
  phi_propagated_ME2 = -9;
  phi_propagated_ME3 = -9;
  phi_propagated_ME4 = -9;
  eta_ME1_TF = -9;
  eta_ME2_TF = -9;
  eta_ME3_TF = -9;
  eta_ME4_TF = -9;
  phi_ME1_TF = -9;
  phi_ME2_TF = -9;
  phi_ME3_TF = -9;
  phi_ME4_TF = -9;
 
  eta_interStat12 = -9;
  phi_interStat12 = -9;
  eta_interStat23 = -9;
  phi_interStat23 = -9;
  eta_interStat13 = -9;
  phi_interStat13 = -9;
  
  allstubs_matched_TF = false; 

  has_l1Extra = 0;
  l1Extra_pt = -99;
  l1Extra_eta = -99;
  l1Extra_phi = -99;
  l1Extra_dR = -99;
  has_recoTrackExtra = 0;
  recoTrackExtra_pt_inner = - 99.;
  recoTrackExtra_eta_inner = - 99.;
  recoTrackExtra_phi_inner = - 99.;  
  recoTrackExtra_pt_outer = - 99.;
  recoTrackExtra_eta_outer = - 99.;
  recoTrackExtra_phi_outer = - 99.;
  has_recoTrack = 0;
  recoTrack_pt_outer = - 99.;
  recoTrack_eta_outer = - 99.;
  recoTrack_phi_outer = - 99.;
  has_recoChargedCandidate = 0;
  recoChargedCandidate_pt = - 99.;
  recoChargedCandidate_eta = - 99.;
  recoChargedCandidate_phi = - 99.;

  recoChargedCandidate_nValidDTHits = 0;
  recoChargedCandidate_nValidCSCHits = 0;
  recoChargedCandidate_nValidRPCHits = 0;

  npar = -1;
  npar_jose = -1;
  pt_position_sh=-1;
  pt_position=-1;
  pt_position_smeared=-1;
  pt_direction_gemcsc = -1;
  pt_direction_gemcsc_central = -1;
  pt_direction_sh=-1;
  pt_direction_gemcsc_smeared=-1;
  pt_direction_xfactor = -1;
  pt_direction_xfactor_smear0 = -1;
  pt_direction_xfactor_smear1 = -1;
  pt_direction_xfactor_smear2 = -1;
  pt_direction_xfactor_smear3 = -1;
  pt_direction_xfactor_L1_1 = -1;
  pt_direction_xfactor_L1_2 = -1;
  hasSt1St2St3=false;
  hasSt3orSt4_sh=false;
}


TTree* MyTrackEff::book(TTree *t, const std::string & name)
{
  edm::Service< TFileService > fs;
  t = fs->make<TTree>(name.c_str(), name.c_str());

  t->Branch("lumi", &lumi);
  t->Branch("run", &run);
  t->Branch("event", &event);

  t->Branch("pt", &pt);
  t->Branch("pz", &pz);
  t->Branch("eta", &eta);
  t->Branch("phi", &phi);
  t->Branch("charge", &charge);
  t->Branch("endcap", &endcap);

  t->Branch("genGd_m",   &genGd_m,   "genGd_m/F");
  t->Branch("genGd_E",   &genGd_E,   "genGd_E/F");
  t->Branch("genGd_p",   &genGd_p,   "genGd_p/F");
  t->Branch("genGd_pt",  &genGd_pt,  "genGd_pt/F");
  t->Branch("genGd_px",  &genGd_px,  "genGd_px/F");
  t->Branch("genGd_py",  &genGd_py,  "genGd_py/F");
  t->Branch("genGd_pz",  &genGd_pz,  "genGd_pz/F");
  t->Branch("genGd_eta", &genGd_eta, "genGd_eta/F");
  t->Branch("genGd_phi", &genGd_phi, "genGd_phi/F");
  t->Branch("genGd_vx",  &genGd_vx,  "genGd_vx/F");
  t->Branch("genGd_vy",  &genGd_vy,  "genGd_vy/F");
  t->Branch("genGd_vz",  &genGd_vz,  "genGd_vz/F");
  t->Branch("genGd_vLx",  &genGd_vLx,  "genGd_vLx/F");
  t->Branch("genGd_vLy",  &genGd_vLy,  "genGd_vLy/F");
  t->Branch("genGd_vLz",  &genGd_vLz,  "genGd_vLz/F");
  t->Branch("genGd_lxy",  &genGd_lxy,  "genGd_lxy/F");
  t->Branch("genGd_l",  &genGd_l,  "genGd_l/F");
  t->Branch("genGd_dxy",  &genGd_dxy,  "genGd_dxy/F");
  t->Branch("genGd0Gd1_dR",  &genGd0Gd1_dR,  "genGd0Gd1_dR/F");
  t->Branch("genGd0Gd1_m",  &genGd0Gd1_m,  "genGd0Gd1_m/F");
  t->Branch("genGdMu_dxy_max",  &genGdMu_dxy_max,  "genGdMu_dxy_max/F");
  t->Branch("genGdMu_eta_max",  &genGdMu_eta_max,  "genGdMu_eta_max/F");

  t->Branch("genGd_index", &genGd_index, "genGd_index/F");
  t->Branch("genGdMu_index", &genGdMu_index, "genGdMu_index/F");
  t->Branch("genGdMu_p", &genGdMu_p, "genGdMu_p/F");
  t->Branch("genGdMu_pt", &genGdMu_pt, "genGdMu_pt/F");
  t->Branch("genGdMu_px", &genGdMu_px, "genGdMu_px/F");
  t->Branch("genGdMu_py", &genGdMu_py, "genGdMu_py/F");
  t->Branch("genGdMu_pz", &genGdMu_pz, "genGdMu_pz/F");
  t->Branch("genGdMu_eta",& genGdMu_eta, "genGdMu_eta/F");
  t->Branch("genGdMu_phi", &genGdMu_phi, "genGdMu_phi/F");
  t->Branch("genGdMu_phi_corr", &genGdMu_phi_corr, "genGdMu_phi_corr/F");
  t->Branch("genGdMu_vx", &genGdMu_vx, "genGdMu_vx/F");
  t->Branch("genGdMu_vy", &genGdMu_vy, "genGdMu_vy/F");
  t->Branch("genGdMu_vz", &genGdMu_vz, "genGdMu_vz/F");
  t->Branch("genGdMu_dxy", &genGdMu_dxy, "genGdMu_dxy/F");
  t->Branch("genGdMu_dR", &genGdMu_dR, "genGdMu_dR/F");


  t->Branch("chamber_ME1_csc_sh", &chamber_ME1_csc_sh);
  t->Branch("chamber_ME2_csc_sh", &chamber_ME2_csc_sh);
  t->Branch("chamber_odd", &chamber_odd);
  t->Branch("chamber_even", &chamber_even);
  t->Branch("quality_odd", &quality_odd);
  t->Branch("quality_even", &quality_even);
  t->Branch("bending_sh", &bending_sh);
  t->Branch("perp_cscsh_odd", &perp_cscsh_odd);
  t->Branch("perp_cscsh_even", &perp_cscsh_even);
  t->Branch("perp_gemsh_odd", &perp_gemsh_odd);
  t->Branch("perp_gemsh_even", &perp_gemsh_even);
  t->Branch("centralperp_gemsh_odd", &centralperp_gemsh_odd);
  t->Branch("centralperp_gemsh_even", &centralperp_gemsh_even);
  t->Branch("phi_cscsh_even", &phi_cscsh_even);
  t->Branch("phi_cscsh_odd", &phi_cscsh_odd);
  t->Branch("phi_layer1_sh_even", &phi_layer1_sh_even);
  t->Branch("phi_layer1_sh_odd", &phi_layer1_sh_odd);
  t->Branch("eta_layer1_sh_even", &eta_layer1_sh_even);
  t->Branch("eta_layer1_sh_odd", &eta_layer1_sh_odd);
  t->Branch("perp_layer1_sh_even", &perp_layer1_sh_even);
  t->Branch("perp_layer1_sh_odd", &perp_layer1_sh_odd);
  t->Branch("phi_layer6_sh_even", &phi_layer6_sh_even);
  t->Branch("phi_layer6_sh_odd", &phi_layer6_sh_odd);
  t->Branch("eta_layer6_sh_even", &eta_layer6_sh_even);
  t->Branch("eta_layer6_sh_odd", &eta_layer6_sh_odd);
  t->Branch("perp_layer6_sh_even", &perp_layer6_sh_even);
  t->Branch("perp_layer6_sh_odd", &perp_layer6_sh_odd);
  t->Branch("eta_cscsh_even", &eta_cscsh_even);
  t->Branch("eta_cscsh_odd", &eta_cscsh_odd);
  t->Branch("pt_sh_even", &pt_sh_even);
  t->Branch("pteta_sh_even", &pteta_sh_even);
  t->Branch("ptphi_sh_even", &ptphi_sh_even);
  t->Branch("pt_sh_odd", &pt_sh_odd);
  t->Branch("pteta_sh_odd", &pteta_sh_odd);
  t->Branch("ptphi_sh_odd", &ptphi_sh_odd);
  t->Branch("eta_st1_sh", &eta_st1_sh);
  t->Branch("eta_st2_sh", &eta_st2_sh);
  t->Branch("ptphi_gemsh_even", &ptphi_gemsh_even);
  t->Branch("ptphi_gemsh_odd", &ptphi_gemsh_odd);
  //t->Branch("ptphi_diff_sh_11", &ptphi_diff_sh_11);
  //t->Branch("ptphi_diff_sh_12", &ptphi_diff_sh_12);
  //t->Branch("ptphi_diff_sh_21", &ptphi_diff_sh_21);
  //t->Branch("ptphi_diff_sh_22", &ptphi_diff_sh_22);
  t->Branch("ptphi_diff_sh", &ptphi_diff_sh);
  t->Branch("csc_bending_angle12_gemcsc", &csc_bending_angle12_gemcsc);
  t->Branch("csc_bending_angle12_xfactor", &csc_bending_angle12_xfactor);
  t->Branch("csc_bending_angle12_xfactor_smear0", &csc_bending_angle12_xfactor_smear0);
  t->Branch("csc_bending_angle12_xfactor_smear1", &csc_bending_angle12_xfactor_smear1);
  t->Branch("csc_bending_angle12_xfactor_smear2", &csc_bending_angle12_xfactor_smear2);
  t->Branch("csc_bending_angle12_xfactor_smear3", &csc_bending_angle12_xfactor_smear3);
  t->Branch("csc_bending_angle12_xfactor_L1_1", &csc_bending_angle12_xfactor_L1_1);
  t->Branch("csc_bending_angle12_xfactor_L1_2", &csc_bending_angle12_xfactor_L1_2);

  t->Branch("phi_layer1_fit_even", &phi_layer1_fit_even);
  t->Branch("phi_layer1_fit_odd", &phi_layer1_fit_odd);
  t->Branch("phi_layer3_fit_even", &phi_layer3_fit_even);
  t->Branch("phi_layer3_fit_odd", &phi_layer3_fit_odd);
  t->Branch("phi_layer6_fit_even", &phi_layer6_fit_even);
  t->Branch("phi_layer6_fit_odd", &phi_layer6_fit_odd);
  t->Branch("has_csc_sh", &has_csc_sh);
  t->Branch("has_csc_strips", &has_csc_strips);
  t->Branch("has_csc_wires", &has_csc_wires);
  t->Branch("has_clct", &has_clct);
  t->Branch("has_alct", &has_alct);
  t->Branch("has_lct", &has_lct);
  t->Branch("bend_lct_odd", &bend_lct_odd);
  t->Branch("bend_lct_even", &bend_lct_even);
  t->Branch("bx_lct_odd", &bx_lct_odd);
  t->Branch("bx_lct_even", &bx_lct_even);
  t->Branch("hs_lct_odd", &hs_lct_odd);
  t->Branch("hs_lct_even", &hs_lct_even);
  t->Branch("wg_lct_even", &wg_lct_even);
  t->Branch("wg_lct_odd", &wg_lct_odd);
  t->Branch("phi_lct_odd", &phi_lct_odd);
  t->Branch("phi_lct_even", &phi_lct_even);
  t->Branch("eta_lct_odd", &eta_lct_odd);
  t->Branch("eta_lct_even", &eta_lct_even);
  t->Branch("perp_lct_odd", &perp_lct_odd);
  t->Branch("perp_lct_even", &perp_lct_even);
  t->Branch("dphi_lct_odd", &dphi_lct_odd);
  t->Branch("dphi_lct_even", &dphi_lct_even);
  t->Branch("passdphi_odd", &passdphi_odd);
  t->Branch("passdphi_even", &passdphi_even);
  
  t->Branch("wiregroup_odd", &wiregroup_odd);
  t->Branch("wiregroup_even", &wiregroup_even);
  t->Branch("halfstrip_odd", &halfstrip_odd);
  t->Branch("halfstrip_even", &halfstrip_even);
  t->Branch("quality_clct_odd", &quality_clct_odd);
  t->Branch("quality_clct_even", &quality_clct_even);
  t->Branch("quality_alct_odd", &quality_alct_odd);
  t->Branch("quality_alct_even", &quality_alct_even);
  t->Branch("nlayers_csc_sh_odd", &nlayers_csc_sh_odd);
  t->Branch("nlayers_csc_sh_even", &nlayers_csc_sh_even);
  t->Branch("nlayers_wg_dg_odd", &nlayers_wg_dg_odd);
  t->Branch("nlayers_wg_dg_even", &nlayers_wg_dg_even);
  t->Branch("nlayers_st_dg_odd", &nlayers_st_dg_odd);
  t->Branch("nlayers_st_dg_even", &nlayers_st_dg_even);

  t->Branch("pad_odd", &pad_odd);
  t->Branch("pad_even", &pad_even);
  t->Branch("Copad_odd", &Copad_odd);
  t->Branch("copad_even", &Copad_even);

  t->Branch("hsfromgem_odd", &hsfromgem_odd);
  t->Branch("hsfromgem_even", &hsfromgem_even);

  t->Branch("has_gem_sh", &has_gem_sh);
  t->Branch("has_gem_sh2", &has_gem_sh2);
  t->Branch("has_gem_dg", &has_gem_dg);
  t->Branch("has_gem_dg2", &has_gem_dg2);
  t->Branch("has_gem_pad", &has_gem_pad);
  t->Branch("has_gem_pad2", &has_gem_pad2);
  t->Branch("has_gem_copad", &has_gem_copad);
  t->Branch("strip_gemsh_odd", &strip_gemsh_odd);
  t->Branch("strip_gemsh_even", &strip_gemsh_even);
  t->Branch("eta_gemsh_odd", &eta_gemsh_odd);
  t->Branch("eta_gemsh_even", &eta_gemsh_even);
  t->Branch("phi_gemsh_odd", &phi_gemsh_odd);
  t->Branch("phi_gemsh_even", &phi_gemsh_even);
  t->Branch("dphi_sh_odd", &dphi_sh_odd);
  t->Branch("dphi_sh_even", &dphi_sh_even);
  t->Branch("dphipositionpt_gemsh_even", &dphipositionpt_gemsh_even);
  t->Branch("dphipositionpt_gemsh_odd", &dphipositionpt_gemsh_odd);
  t->Branch("dphipositionpt_cscsh_even", &dphipositionpt_cscsh_even);
  t->Branch("dphipositionpt_cscsh_odd", &dphipositionpt_cscsh_odd);
  t->Branch("phiM_gemcsc_odd", &phiM_gemcsc_odd);
  t->Branch("phiM_gemcsc_even", &phiM_gemcsc_even);
  t->Branch("phiM_gemcsc_central_odd", &phiM_gemcsc_central_odd);
  t->Branch("phiM_gemcsc_central_even", &phiM_gemcsc_central_even);
  t->Branch("phiM_fitxyz_odd", &phiM_fitxyz_odd);
  t->Branch("phiM_fitxyz_even", &phiM_fitxyz_even);
  t->Branch("strip_gemdg_odd", &strip_gemdg_odd);
  t->Branch("strip_gemdg_even", &strip_gemdg_even);

  t->Branch("has_rpc_sh", &has_rpc_sh);
  t->Branch("has_rpc_dg", &has_rpc_dg);
  t->Branch("strip_rpcdg_odd", &strip_rpcdg_odd);
  t->Branch("strip_rpcdg_even", &strip_rpcdg_even);
  t->Branch("hsfromrpc_odd", &hsfromrpc_odd);
  t->Branch("hsfromrpc_even", &hsfromrpc_even);

  t->Branch("bx_pad_odd", &bx_pad_odd);
  t->Branch("bx_pad_even", &bx_pad_even);
  t->Branch("phi_pad_odd", &phi_pad_odd);
  t->Branch("phi_pad_even", &phi_pad_even);
  t->Branch("eta_pad_odd", &eta_pad_odd);
  t->Branch("eta_pad_even", &eta_pad_even);
  t->Branch("dphi_pad_odd", &dphi_pad_odd);
  t->Branch("dphi_pad_even", &dphi_pad_even);
  t->Branch("deta_pad_odd", &deta_pad_odd);
  t->Branch("deta_pad_even", &deta_pad_even);

  t->Branch("bx_rpcstrip_odd", &bx_rpcstrip_odd);
  t->Branch("bx_rpcstrip_even", &bx_rpcstrip_even);
  t->Branch("phi_rpcstrip_odd", &phi_rpcstrip_odd);
  t->Branch("phi_rpcstrip_even", &phi_rpcstrip_even);
  t->Branch("eta_rpcstrip_odd", &eta_rpcstrip_odd);
  t->Branch("eta_rpcstrip_even", &eta_rpcstrip_even);
  t->Branch("dphi_rpcstrip_odd", &dphi_rpcstrip_odd);
  t->Branch("dphi_rpcstrip_even", &dphi_rpcstrip_even);
  t->Branch("deta_rpcstrip_odd", &deta_rpcstrip_odd);
  t->Branch("deta_rpcstrip_even", &deta_rpcstrip_even);

  //t->Branch("", &);
  t->Branch("has_tfTrack", &has_tfTrack);
  t->Branch("has_tfCand", &has_tfCand);
  t->Branch("has_gmtRegCand", &has_gmtRegCand);
  t->Branch("has_gmtCand", &has_gmtCand);
  
  //csctftrack
  t->Branch("trackpt", &trackpt);
  t->Branch("tracketa", &tracketa);
  t->Branch("trackphi", &trackphi);
  t->Branch("quality_packed",&quality_packed);
  t->Branch("rank",&rank);
  t->Branch("pt_packed",&pt_packed);
  t->Branch("eta_packed",&eta_packed);
  t->Branch("phi_packed",&phi_packed);
  t->Branch("chargesign",&chargesign);
  t->Branch("deltaphi12",&deltaphi12);
  t->Branch("deltaphi23",&deltaphi23);
  t->Branch("hasME1",&hasME1);
  t->Branch("hasME2",&hasME2);
  t->Branch("ME1_ring",&ME1_ring);
  t->Branch("ME2_ring",&ME2_ring);
  t->Branch("chamberME1",&chamberME1);
  t->Branch("chamberME2",&chamberME2);
  t->Branch("ME1_hs",&ME1_hs);
  t->Branch("ME1_wg",&ME1_wg);
  t->Branch("ME2_hs",&ME2_hs);
  t->Branch("ME2_wg",&ME2_wg);
  t->Branch("dphiGE11",&dphiGE11);
  t->Branch("dphiGE21",&dphiGE21);
  t->Branch("passGE11",&passGE11);
  t->Branch("passGE11_pt5",&passGE11_pt5);
  t->Branch("passGE11_pt7",&passGE11_pt7);
  t->Branch("passGE11_pt10",&passGE11_pt10);
  t->Branch("passGE11_pt15",&passGE11_pt15);
  t->Branch("passGE11_pt20",&passGE11_pt20);
  t->Branch("passGE11_pt30",&passGE11_pt30);
  t->Branch("passGE11_pt40",&passGE11_pt40);
  t->Branch("passGE21",&passGE21);
  t->Branch("passGE21_pt5",&passGE21_pt5);
  t->Branch("passGE21_pt7",&passGE21_pt7);
  t->Branch("passGE21_pt10",&passGE21_pt10);
  t->Branch("passGE21_pt15",&passGE21_pt15);
  t->Branch("passGE21_pt20",&passGE21_pt20);
  t->Branch("passGE21_pt30",&passGE21_pt30);
  t->Branch("passGE21_pt40",&passGE21_pt40);
  t->Branch("passGE11_simpt",&passGE11_simpt);
  t->Branch("passGE21_simpt",&passGE21_simpt);
  t->Branch("nstubs",&nstubs);
  t->Branch("deltaR",&deltaR);
  t->Branch("lctdphi12",&lctdphi12);
   
  t->Branch("eta_propagated_ME1",&eta_propagated_ME1);
  t->Branch("eta_propagated_ME2",&eta_propagated_ME2);
  t->Branch("eta_propagated_ME3",&eta_propagated_ME3);
  t->Branch("eta_propagated_ME4",&eta_propagated_ME4);
  t->Branch("phi_propagated_ME1",&phi_propagated_ME1);
  t->Branch("phi_propagated_ME2",&phi_propagated_ME2);
  t->Branch("phi_propagated_ME3",&phi_propagated_ME3);
  t->Branch("phi_propagated_ME4",&phi_propagated_ME4);
  t->Branch("eta_ME1_TF",&eta_ME1_TF);
  t->Branch("eta_ME2_TF",&eta_ME2_TF);
  t->Branch("eta_ME3_TF",&eta_ME3_TF);
  t->Branch("eta_ME4_TF",&eta_ME4_TF);
  t->Branch("phi_ME1_TF",&phi_ME1_TF);
  t->Branch("phi_ME2_TF",&phi_ME2_TF);
  t->Branch("phi_ME3_TF",&phi_ME3_TF);
  t->Branch("phi_ME4_TF",&phi_ME4_TF);

  t->Branch("eta_interStat12",&eta_interStat12);
  t->Branch("phi_interStat12",&phi_interStat12);
  t->Branch("eta_interStat23",&eta_interStat23);
  t->Branch("phi_interStat23",&phi_interStat23);
  t->Branch("eta_interStat13",&eta_interStat13);
  t->Branch("phi_interStat13",&phi_interStat13);
  
  t->Branch("allstubs_matched_TF",&allstubs_matched_TF);

  t->Branch("has_l1Extra", &has_l1Extra);
  t->Branch("l1Extra_pt", &l1Extra_pt);
  t->Branch("l1Extra_eta", &l1Extra_eta);
  t->Branch("l1Extra_phi", &l1Extra_phi);
  t->Branch("l1Extra_dR", &l1Extra_dR);
  t->Branch("has_recoTrackExtra", &has_recoTrackExtra);
  t->Branch("recoTrackExtra_pt_inner", &recoTrackExtra_pt_inner);
  t->Branch("recoTrackExtra_eta_inner", &recoTrackExtra_eta_inner);
  t->Branch("recoTrackExtra_phi_inner", &recoTrackExtra_phi_inner);
  t->Branch("recoTrackExtra_pt_outer", &recoTrackExtra_pt_outer);
  t->Branch("recoTrackExtra_eta_outer", &recoTrackExtra_eta_outer);
  t->Branch("recoTrackExtra_phi_outer", &recoTrackExtra_phi_outer);
  t->Branch("has_recoTrack", &has_recoTrack);
  t->Branch("recoTrack_pt_outer", &recoTrack_pt_outer);
  t->Branch("recoTrack_eta_outer", &recoTrack_eta_outer);
  t->Branch("recoTrack_phi_outer", &recoTrack_phi_outer);
  t->Branch("has_recoChargedCandidate", &has_recoChargedCandidate);
  t->Branch("recoChargedCandidate_pt", &recoChargedCandidate_pt);
  t->Branch("recoChargedCandidate_eta", &recoChargedCandidate_eta);
  t->Branch("recoChargedCandidate_phi", &recoChargedCandidate_phi); 

  t->Branch("recoChargedCandidate_nValidDTHits", &recoChargedCandidate_nValidDTHits); 
  t->Branch("recoChargedCandidate_nValidCSCHits", &recoChargedCandidate_nValidCSCHits); 
  t->Branch("recoChargedCandidate_nValidRPCHits", &recoChargedCandidate_nValidRPCHits); 

  
  t->Branch("npar", &npar); 
  t->Branch("npar_jose", &npar_jose); 
  t->Branch("pt_position_sh", &pt_position_sh); 
  t->Branch("pt_position", &pt_position); 
  t->Branch("pt_position_smeared", &pt_position_smeared); 
  t->Branch("pt_direction_sh", &pt_direction_sh); 
  t->Branch("pt_direction_gemcsc", &pt_direction_gemcsc); 
  t->Branch("pt_direction_xfactor", &pt_direction_xfactor); 
  t->Branch("pt_direction_xfactor_smear0", &pt_direction_xfactor_smear0); 
  t->Branch("pt_direction_xfactor_smear1", &pt_direction_xfactor_smear1); 
  t->Branch("pt_direction_xfactor_smear2", &pt_direction_xfactor_smear2); 
  t->Branch("pt_direction_xfactor_smear3", &pt_direction_xfactor_smear3); 
  t->Branch("pt_direction_xfactor_L1_1", &pt_direction_xfactor_L1_1); 
  t->Branch("pt_direction_xfactor_L1_2", &pt_direction_xfactor_L1_2); 
  t->Branch("pt_direction_gemcsc_central", &pt_direction_gemcsc_central); 
  t->Branch("pt_direction_gemcsc_smeared", &pt_direction_gemcsc_smeared); 
  t->Branch("hasSt1St2St3", &hasSt1St2St3); 
  t->Branch("hasSt3orSt4_sh", &hasSt3orSt4_sh); 

  return t;
}

// --------------------------- GEMCSCAnalyzer ---------------------------

class GEMCSCAnalyzer : public edm::EDAnalyzer
{
public:

  explicit GEMCSCAnalyzer(const edm::ParameterSet&);

  ~GEMCSCAnalyzer() {}
  
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  void bookSimTracksDeltaTree();

  void analyzeTrackChamberDeltas(SimTrackMatchManager& match, int trk_no);
  void analyzeTrackEff(SimTrackMatchManager& match, int trk_no);
  void printout(SimTrackMatchManager& match, int trk_no, const char msg[300]);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);
  
  edm::ParameterSet cfg_;
  edm::InputTag simInputLabel_;
  int verboseSimTrack_;
  double simTrackMinPt_;
  double simTrackMinEta_;
  double simTrackMaxEta_;
  double simTrackOnlyMuon_;
  int verbose_;
  bool ntupleTrackChamberDelta_;
  bool ntupleTrackEff_;
  bool matchprint_;
  double bendingcutPt_;
  std::vector<string> cscStations_;
  std::vector<std::pair<int,int> > cscStationsCo_;
  std::set<int> stations_to_use_;

  TTree *tree_eff_[12]; // for up to 9 stations
  TTree *tree_delta_;
  
  MyTrackEff  etrk_[12];
  MyTrackChamberDelta dtrk_;

  int minNHitsChamberCSCSimHit_;
  int minNHitsChamberCSCWireDigi_;
  int minNHitsChamberCSCStripDigi_;
  int minNHitsChamberCLCT_;
  int minNHitsChamberALCT_;
  int minNHitsChamberLCT_;
  int minNHitsChamberMPLCT_;
};


GEMCSCAnalyzer::GEMCSCAnalyzer(const edm::ParameterSet& ps)
: cfg_(ps.getParameterSet("simTrackMatching"))
, verbose_(ps.getUntrackedParameter<int>("verbose", 0))
{
  cscStations_ = cfg_.getParameter<std::vector<string> >("cscStations");
  ntupleTrackChamberDelta_ = cfg_.getParameter<bool>("ntupleTrackChamberDelta");
  ntupleTrackEff_ = cfg_.getParameter<bool>("ntupleTrackEff");
  //matchprint_ = false; 
  matchprint_ =  cfg_.getParameter<bool>("matchprint");
  bendingcutPt_ = cfg_.getUntrackedParameter<double>("bendingcutPt",10);

  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrack.getParameter<int>("verbose");
  simInputLabel_ = edm::InputTag("g4SimHits"); //simTrack.getParameter<edm::InputTag>("input");
  simTrackMinPt_ = simTrack.getParameter<double>("minPt");
  simTrackMinEta_ = simTrack.getParameter<double>("minEta");
  simTrackMaxEta_ = simTrack.getParameter<double>("maxEta");
  simTrackOnlyMuon_ = simTrack.getParameter<bool>("onlyMuon");
    
  auto cscSimHit = cfg_.getParameter<edm::ParameterSet>("cscSimHit");
  minNHitsChamberCSCSimHit_ = cscSimHit.getParameter<int>("minNHitsChamber");

  auto cscWireDigi = cfg_.getParameter<edm::ParameterSet>("cscWireDigi");
  minNHitsChamberCSCWireDigi_ = cscWireDigi.getParameter<int>("minNHitsChamber");

  auto cscComparatorDigi = cfg_.getParameter<edm::ParameterSet>("cscStripDigi");
  minNHitsChamberCSCStripDigi_ = cscComparatorDigi.getParameter<int>("minNHitsChamber");

  auto cscCLCT = cfg_.getParameter<edm::ParameterSet>("cscCLCT");
  minNHitsChamberCLCT_ = cscCLCT.getParameter<int>("minNHitsChamber");

  auto cscALCT = cfg_.getParameter<edm::ParameterSet>("cscALCT");
  minNHitsChamberALCT_ = cscALCT.getParameter<int>("minNHitsChamber");

  auto cscLCT = cfg_.getParameter<edm::ParameterSet>("cscLCT");
  minNHitsChamberLCT_ = cscLCT.getParameter<int>("minNHitsChamber");

  auto cscMPLCT = cfg_.getParameter<edm::ParameterSet>("cscMPLCT");
  minNHitsChamberMPLCT_ = cscMPLCT.getParameter<int>("minNHitsChamber");

  if (ntupleTrackChamberDelta_) bookSimTracksDeltaTree();
  if (ntupleTrackEff_)
  {
    vector<int> stations = cfg_.getParameter<vector<int> >("cscStationsToUse");
    copy(stations.begin(), stations.end(), inserter(stations_to_use_, stations_to_use_.end()) );

    for(auto s: stations_to_use_)
    {
      stringstream ss;
      ss << "trk_eff_"<< cscStations_[s];
      std::cout <<"station to use "<< cscStations_[s]  << std::endl;
      tree_eff_[s] = etrk_[s].book(tree_eff_[s], ss.str());
    }
  }

  cscStationsCo_.push_back(std::make_pair(-99,-99));
  cscStationsCo_.push_back(std::make_pair(1,-99));
  cscStationsCo_.push_back(std::make_pair(1,4));
  cscStationsCo_.push_back(std::make_pair(1,1));
  cscStationsCo_.push_back(std::make_pair(1,2));
  cscStationsCo_.push_back(std::make_pair(1,3));
  cscStationsCo_.push_back(std::make_pair(2,1));
  cscStationsCo_.push_back(std::make_pair(2,2));
  cscStationsCo_.push_back(std::make_pair(3,1));
  cscStationsCo_.push_back(std::make_pair(3,2));
  cscStationsCo_.push_back(std::make_pair(4,1));
  cscStationsCo_.push_back(std::make_pair(4,2));
}


int GEMCSCAnalyzer::detIdToMEStation(int st, int ri)
{
  auto p(std::make_pair(st, ri));
  return std::find(cscStationsCo_.begin(), cscStationsCo_.end(), p) - cscStationsCo_.begin();
}


void GEMCSCAnalyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
}


bool GEMCSCAnalyzer::isSimTrackGood(const SimTrack &t)
{
  // SimTrack selection
  if (t.noVertex()) return false;
  if (t.noGenpart()) return false;
  // only muons 
  if (std::abs(t.type()) != 13 and simTrackOnlyMuon_) return false;
  // pt selection
  if (t.momentum().pt() < simTrackMinPt_) return false;
  // eta selection
  const float eta(std::abs(t.momentum().eta()));
  if (eta > simTrackMaxEta_ || eta < simTrackMinEta_) return false; 
  return true;
}


void GEMCSCAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  edm::Handle<edm::SimTrackContainer> sim_tracks;
  ev.getByLabel(simInputLabel_, sim_tracks);
  const edm::SimTrackContainer & sim_track = *sim_tracks.product();

  edm::Handle<edm::SimVertexContainer> sim_vertices;
  ev.getByLabel(simInputLabel_, sim_vertices);
  const edm::SimVertexContainer & sim_vert = *sim_vertices.product();

  if (verboseSimTrack_){
    std::cout << "Total number of SimTrack in this event: " << sim_track.size() << std::endl;      
  }
    
  int trk_no=0;
  for (auto& t: sim_track)
  {
    if (!isSimTrackGood(t)) continue;
    if (verboseSimTrack_){
      std::cout << "Processing SimTrack " << trk_no + 1 << std::endl;      
      std::cout << "pt(GeV/c) = " << t.momentum().pt() << ", eta = " << t.momentum().eta()  
                << ", phi = " << t.momentum().phi() << ", Q = " << t.charge() << std::endl;
    }
    
//    std::cout<< " initialize SimTrackMatcherManager "<< std::endl;  
    // match hits and digis to this SimTrack
    SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);

    if (ntupleTrackChamberDelta_) analyzeTrackChamberDeltas(match, trk_no);
    if (ntupleTrackEff_) analyzeTrackEff(match, trk_no);
    ++trk_no;

    // if (matchprint_) printout(match, trk_no);    
        
   // bool has_tftracks(etrk_[0].has_tfTrack>0 && !etrk_[0].allstubs_matched_TF && etrk_[0].nstubs>2);
    //bool phi_diff_ME1(abs(etrk_[0].phi_propagated_ME1-etrk_[0].phi_ME1_TF) < 0.02); 
    //bool phi_diff_ME2(abs(etrk_[0].phi_interStat12-etrk_[0].phi_ME2_TF) < 0.02); 
    //bool phi_diff_ME3(abs(etrk_[0].phi_interStat23-etrk_[0].phi_ME3_TF) < 0.02); 
    // if (has_csc_sh_odd || has_csc_sh_even)  std::cout <<"st1 has_csc_sh " << std::endl;
    // if (has_alct_odd || has_alct_even)   std::cout <<"  st1 has_alct " << std::endl;
    //bool Debug (etrk_[10].has_alct>0 && (etrk_[10].has_clct>0 || etrk_[10].has_rpc_dg>0) && etrk_[10].has_lct == 0);
    //bool Debug (etrk_[0].allstubs_matched_TF==0);
    //bool Debug (etrk_[1].pt<10 and (fabs(etrk_[1].phi_lct_even-etrk_[1].phi_cscsh_even)>0.008 and etrk_[1].phi_lct_even>-4 and etrk_[1].phi_cscsh_even>-4) or (fabs(etrk_[1].phi_lct_odd-etrk_[1].phi_cscsh_odd)>0.008 and etrk_[1].phi_lct_odd>-4 and etrk_[1].phi_cscsh_odd>-4));
    bool Debug ((fabs(etrk_[1].dphi_sh_odd)>0.5 and fabs(etrk_[1].dphi_sh_odd)<9) or (fabs(etrk_[1].dphi_sh_even)>0.5 and fabs(etrk_[1].dphi_sh_even)<9));
    if (matchprint_ and Debug){
    	std::cout <<"ME11 phi_cscsh even "<<etrk_[1].phi_cscsh_even <<" odd "<<etrk_[1].phi_cscsh_odd<<" phi_gemsh even "<< etrk_[1].phi_gemsh_even <<" odd "<< etrk_[1].phi_gemsh_odd<<" dphi_sh even "<< etrk_[1].dphi_sh_even<<" odd "<< etrk_[1].dphi_sh_odd <<std::endl;
	printout(match, trk_no, "to debug dephi at sim level");
       } 
    
  }
}



void GEMCSCAnalyzer::analyzeTrackEff(SimTrackMatchManager& match, int trk_no)
{

  const DisplacedGENMuonMatcher& match_gen = match.genMuons();
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const RPCDigiMatcher& match_rd = match.rpcDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  const TrackMatcher& match_track = match.tracks();
  const L1GlobalMuonTriggerMatcher& match_l1_gmt = match.l1GMTCands();
  const HLTTrackMatcher& match_hlt_track = match.hltTracks();
  const SimTrack &t = match_sh.trk();
   

  auto matchedDarkBoson(match_gen.getMatchedDarkBoson());
  auto matchedGENMuon(match_gen.getMatchedGENMuon());
  /*
  if (match_gen.checkRunOK()){
  	std::cout <<" matched dark boson mass "<< matchedDarkBoson->mass() << std::endl;
 	std::cout <<"matched muon pt "<< matchedGENMuon->pt()<<" phi "<< matchedGENMuon->phi()<<" eta "<< matchedGENMuon->eta() <<" dR "<< match_gen.matchedGenMudR()<<" dxy "<< match_gen.matchedGenMudxy() <<std::endl;
  }
  */
  std::cout <<"Sim trk_no " << trk_no <<" eta "<< t.momentum().eta() << " phi "<< t.momentum().phi() << " pt "<< t.momentum().pt()<<" pz "<<  t.momentum().pz()<<std::endl;

  for (auto s: stations_to_use_)
  {

    etrk_[s].init();
    etrk_[s].run = match.simhits().event().id().run();
    etrk_[s].lumi = match.simhits().event().id().luminosityBlock();
    etrk_[s].event = match.simhits().event().id().event();
    etrk_[s].pt = t.momentum().pt();
    etrk_[s].pz = t.momentum().pz();
    etrk_[s].phi = t.momentum().phi();
    etrk_[s].eta = t.momentum().eta();
    etrk_[s].charge = t.charge();
    etrk_[s].endcap = (etrk_[s].eta > 0.) ? 1 : -1;

  if (match_gen.checkRunOK()){
    // Dark photon  
    etrk_[s].genGd_m = matchedDarkBoson->mass();
    etrk_[s].genGd_E = matchedDarkBoson->energy();
    etrk_[s].genGd_p = matchedDarkBoson->p();
    etrk_[s].genGd_pt = matchedDarkBoson->pt();
    etrk_[s].genGd_px = matchedDarkBoson->px();
    etrk_[s].genGd_py = matchedDarkBoson->py();
    etrk_[s].genGd_pz = matchedDarkBoson->pz();
    etrk_[s].genGd_eta = matchedDarkBoson->eta();
    etrk_[s].genGd_phi = matchedDarkBoson->phi();
    etrk_[s].genGd_vx = matchedDarkBoson->vx();
    etrk_[s].genGd_vy = matchedDarkBoson->vy();
    etrk_[s].genGd_vz = matchedDarkBoson->vz();

     etrk_[s].genGd_index = match_gen.darkBosonIndex();
     etrk_[s].genGdMu_index = match_gen.genMuonIndex();
     etrk_[s].genGdMu_p = matchedGENMuon->p();
     etrk_[s].genGdMu_pt = matchedGENMuon->pt();
     etrk_[s].genGdMu_px = matchedGENMuon->px();
     etrk_[s].genGdMu_py = matchedGENMuon->py();
     etrk_[s].genGdMu_pz = matchedGENMuon->pz();
     etrk_[s].genGdMu_eta = matchedGENMuon->eta();
     etrk_[s].genGdMu_phi = matchedGENMuon->phi();
     etrk_[s].genGdMu_vx = matchedGENMuon->vx();
     etrk_[s].genGdMu_vy = matchedGENMuon->vy();
     etrk_[s].genGdMu_vz = matchedGENMuon->vz();
     etrk_[s].genGdMu_dxy = match_gen.matchedGenMudxy();
     etrk_[s].genGdMu_dR = match_gen.matchedGenMudR();
  }

  }

  int chargesign = (t.charge()>0? 1:0);
  float pt = t.momentum().pt();
  // SimHits
  auto csc_simhits(match_sh.chamberIdsCSC(0));
  GlobalPoint gp_sh_odd[12];
  GlobalPoint gp_sh_even[12];
  GlobalVector gv_sh_odd[12];
  GlobalVector gv_sh_even[12];
  for(auto d: csc_simhits)
  {

    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers(match_sh.nLayersWithHitsInSuperChamber(d));
    if (id.station() == 1 and id.chamber()%2 == 1) etrk_[0].chamber_ME1_csc_sh |= 1;
    if (id.station() == 1 and id.chamber()%2 == 0) etrk_[0].chamber_ME1_csc_sh |= 2;
    if (id.station() == 2 and id.chamber()%2 == 1) etrk_[0].chamber_ME2_csc_sh |= 1;
    if (id.station() == 2 and id.chamber()%2 == 0) etrk_[0].chamber_ME2_csc_sh |= 2;
    // case ME11
    if (id.station()==1 and (id.ring()==4 or id.ring()==1)){
      // get the detId of the pairing subchamber
      int other_ring(id.ring()==4 ? 1 : 4);
      CSCDetId co_id(id.endcap(), id.station(), other_ring, id.chamber());
      // check if co_id occurs in the list
      // add the hit layers
     
      auto rawId(co_id.rawId());
      if (csc_simhits.find(rawId) != csc_simhits.end()) {
	nlayers = nlayers+match_sh.nLayersWithHitsInSuperChamber(rawId);

      } 
      
    }
    
    if (nlayers < minNHitsChamberCSCSimHit_) continue;
    GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    GlobalPoint gp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    etrk_[st].bending_sh = match_sh.LocalBendingInChamber(d);
    const CSCDetId id1(id.endcap(), id.station(), id.ring(), id.chamber(), 1);
    const CSCDetId id6(id.endcap(), id.station(), id.ring(), id.chamber(), 6);
    GlobalPoint gp1 = match_sh.simHitsMeanPosition(match_sh.hitsInDetId(id1.rawId()));
    GlobalPoint gp6 = match_sh.simHitsMeanPosition(match_sh.hitsInDetId(id6.rawId()));
    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_csc_sh |= 1;
    else etrk_[st].has_csc_sh |= 2;
    if (odd){
    	etrk_[st].pteta_sh_odd = ym.eta();
    	etrk_[st].ptphi_sh_odd = ym.phi();
    	etrk_[st].pt_sh_odd = ym.perp();
    }else{
    	etrk_[st].pteta_sh_even = ym.eta();
    	etrk_[st].ptphi_sh_even = ym.phi();
    	etrk_[st].pt_sh_even = ym.perp();
    }
    if (odd){
    	if (match_sh.hitsInDetId(id1.rawId()).size()>0){
		etrk_[st].eta_layer1_sh_odd = gp1.eta();
		etrk_[st].phi_layer1_sh_odd = gp1.phi();
		etrk_[st].perp_layer1_sh_odd = gp1.perp();
		//std::cout <<"layer1 id "<< id1 <<" phi "<< gp1.phi() <<" eta "<< gp1.eta() <<" x "<< gp1.x()<<" y "<< gp1.y()<<" z "<< gp1.z() <<" perp "<< gp1.perp() << std::endl;
		}
    	if (match_sh.hitsInDetId(id6.rawId()).size()>0){
		etrk_[st].eta_layer6_sh_odd = gp6.eta();
		etrk_[st].phi_layer6_sh_odd = gp6.phi();
		etrk_[st].perp_layer6_sh_odd = gp6.perp();
		//std::cout <<"layer6 id "<< id6 <<" phi "<< gp6.phi() <<" eta "<< gp6.eta() <<" x "<< gp6.x()<<" y "<< gp6.y()<<" z "<< gp6.z() <<" perp "<< gp6.perp() << std::endl;
		}
    }else{
    	if (match_sh.hitsInDetId(id1.rawId()).size()>0){
		etrk_[st].eta_layer1_sh_even = gp1.eta();
		etrk_[st].phi_layer1_sh_even = gp1.phi();
		etrk_[st].perp_layer1_sh_even = gp1.perp();
		//std::cout <<"layer1 id "<< id1 <<" phi "<< gp1.phi() <<" eta "<< gp1.eta() <<" x "<< gp1.x()<<" y "<< gp1.y()<<" z "<< gp1.z() <<" perp "<< gp1.perp() << std::endl;
		}
    	if (match_sh.hitsInDetId(id6.rawId()).size()>0){
		etrk_[st].eta_layer6_sh_even = gp6.eta();
		etrk_[st].phi_layer6_sh_even = gp6.phi();
		etrk_[st].perp_layer6_sh_even = gp6.perp();
		//std::cout <<"layer6 id "<< id6 <<" phi "<< gp6.phi() <<" eta "<< gp6.eta() <<" x "<< gp6.x()<<" y "<< gp6.y()<<" z "<< gp6.z() <<" perp "<< gp6.perp() << std::endl;
		}
    }

    if (odd) etrk_[st].nlayers_csc_sh_odd = nlayers;
    else etrk_[st].nlayers_csc_sh_even = nlayers;
    
    if (odd) gp_sh_odd[st] = gp;
    else gp_sh_even[st] = gp;
    if (odd) gv_sh_odd[st] = ym;
    else gv_sh_even[st] = ym;
    //std::cout <<" CSC id "<< id <<" gp.eta "<< gp.eta() <<" gp.phi "<< gp.phi() << std::endl;
    //std::cout <<"csc id "<< id << " simhit direction eta "<< ym.eta()<< " phi "<< ym.phi() <<" perp "<< ym.perp()<<" x "<< ym.x() <<" y "<< ym.y() <<" z "<< ym.z() << std::endl; 
    for (int layer=1; layer <= 6; layer++){
	//should be just layer3 or layer 4
    	const CSCDetId csckeyid(id.endcap(), id.station(), id.ring(), id.chamber(), layer); 
        GlobalPoint keygp = match_sh.simHitsMeanPosition(match_sh.hitsInDetId(csckeyid.rawId()));	
	if (match_sh.hitsInDetId(csckeyid.rawId()).size()>0){
	    	//std::cout <<"simhit csc id "<< csckeyid << " eta "<< keygp.eta()<< " phi "<< keygp.phi() <<" perp "<< keygp.perp() <<" x "<< keygp.x() <<" y "<< keygp.y() <<" z "<< keygp.z()<< std::endl; 
    		if (odd) etrk_[st].eta_cscsh_odd = keygp.eta();
    		else     etrk_[st].eta_cscsh_even = keygp.eta();
    		if (odd) etrk_[st].phi_cscsh_odd = keygp.phi();
    		else     etrk_[st].phi_cscsh_even = keygp.phi();
    		if (odd) etrk_[st].perp_cscsh_odd = keygp.perp();
    		else     etrk_[st].perp_cscsh_even = keygp.perp();
    		if (odd) etrk_[st].dphipositionpt_cscsh_odd = deltaPhi(keygp.phi(),etrk_[st].ptphi_sh_odd);
    		else     etrk_[st].dphipositionpt_cscsh_even = deltaPhi(keygp.phi(),etrk_[st].ptphi_sh_even);

		if (st==2 or st==3){
      			if (odd) etrk_[1].eta_cscsh_odd = keygp.eta();
      			else     etrk_[1].eta_cscsh_even = keygp.eta();
      			if (odd) etrk_[1].phi_cscsh_odd = keygp.phi();
      			else     etrk_[1].phi_cscsh_even = keygp.phi();
      			if (odd) etrk_[1].perp_cscsh_odd = keygp.perp();
      			else     etrk_[1].perp_cscsh_even = keygp.perp();
    			if (odd) etrk_[1].dphipositionpt_cscsh_odd = deltaPhi(keygp.phi(),etrk_[st].ptphi_sh_odd);
    			else     etrk_[1].dphipositionpt_cscsh_even = deltaPhi(keygp.phi(),etrk_[st].ptphi_sh_even);
		}
		//std::cout <<" csckeyid "<< csckeyid<<" simhits size "<< match_sh.hitsInDetId(csckeyid.rawId()).size() <<" keygp.eta "<< keygp.eta() << " keygp.phi "<< keygp.phi() << std::endl;
		if (layer >= 3){
		    	//std::cout <<"break!!! layer "<< layer << std::endl;
			break;
		}
	}
    } 

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_csc_sh |= 1;
      else etrk_[1].has_csc_sh |= 2;

      if (odd) etrk_[1].nlayers_csc_sh_odd = nlayers;
      else etrk_[1].nlayers_csc_sh_even = nlayers;

      if (odd) gp_sh_odd[1] = gp;
      else gp_sh_even[1] = gp;
      if (odd) gv_sh_odd[1] = ym;
      else gv_sh_even[1] = ym;


      if (odd){
      	etrk_[1].pt_sh_odd = ym.perp();
      	etrk_[1].pteta_sh_odd = ym.eta();
      	etrk_[1].ptphi_sh_odd = ym.phi();
      }else {
      	etrk_[1].pt_sh_even = ym.perp();
      	etrk_[1].pteta_sh_even = ym.eta();
      	etrk_[1].ptphi_sh_even = ym.phi();
      }
    if (odd){
    	if (match_sh.hitsInDetId(id1.rawId()).size()>0){
		etrk_[1].eta_layer1_sh_odd = gp1.eta();
		etrk_[1].phi_layer1_sh_odd = gp1.phi();
		etrk_[1].perp_layer1_sh_odd = gp1.perp();
		}
    	if (match_sh.hitsInDetId(id6.rawId()).size()>0){
		etrk_[1].eta_layer6_sh_odd = gp6.eta();
		etrk_[1].phi_layer6_sh_odd = gp6.phi();
		etrk_[1].perp_layer6_sh_odd = gp6.perp();
		}
    }else{
    	if (match_sh.hitsInDetId(id1.rawId()).size()>0){
		etrk_[1].eta_layer1_sh_even = gp1.eta();
		etrk_[1].phi_layer1_sh_even = gp1.phi();
		etrk_[1].perp_layer1_sh_even = gp1.perp();
		}
    	if (match_sh.hitsInDetId(id6.rawId()).size()>0){
		etrk_[1].eta_layer6_sh_even = gp6.eta();
		etrk_[1].phi_layer6_sh_even = gp6.phi();
		etrk_[1].perp_layer6_sh_even = gp6.perp();
		}
    }
      etrk_[1].bending_sh = match_sh.LocalBendingInChamber(d);
    }

  }

  // CSC strip digis
  for(auto d: match_cd.chamberIdsStrip(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const int nlayers(match_cd.nLayersWithStripInChamber(d));
    if (nlayers < minNHitsChamberCSCStripDigi_) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_csc_strips |= 1;
    else etrk_[st].has_csc_strips |= 2;

    if (odd) etrk_[st].nlayers_st_dg_odd = nlayers;
    else etrk_[st].nlayers_st_dg_even = nlayers;
    
    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_csc_strips |= 1;
      else etrk_[1].has_csc_strips |= 2;

      if (odd) etrk_[1].nlayers_st_dg_odd = nlayers;
      else etrk_[1].nlayers_st_dg_even = nlayers;
    }  
  }

  // CSC wire digis
  for(auto d: match_cd.chamberIdsWire(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const int nlayers(match_cd.nLayersWithWireInChamber(d));
    if (nlayers < minNHitsChamberCSCWireDigi_) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_csc_wires |= 1;
    else etrk_[st].has_csc_wires |= 2;

    if (odd) etrk_[st].nlayers_wg_dg_odd = nlayers;
    else etrk_[st].nlayers_wg_dg_even = nlayers;

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_csc_wires |= 1;
      else etrk_[1].has_csc_wires |= 2;

      if (odd) etrk_[1].nlayers_wg_dg_odd = nlayers;
      else etrk_[1].nlayers_wg_dg_even = nlayers;
    }  
  }

  // CSC CLCTs
  for(auto d: match_lct.chamberIdsCLCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    auto clct = match_lct.clctInChamber(d);

    if (odd) etrk_[st].halfstrip_odd = digi_channel(clct);
    else etrk_[st].halfstrip_even = digi_channel(clct);

    if (odd) etrk_[st].quality_clct_odd = digi_quality(clct);
    else etrk_[st].quality_clct_even = digi_quality(clct);

    if (odd) etrk_[st].has_clct |= 1;
    else etrk_[st].has_clct |= 2;

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].halfstrip_odd = digi_channel(clct);
      else etrk_[1].halfstrip_even = digi_channel(clct);

      if (odd) etrk_[1].quality_clct_odd = digi_quality(clct);
      else etrk_[1].quality_clct_even = digi_quality(clct);
      
      if (odd) etrk_[1].has_clct |= 1;
      else etrk_[1].has_clct |= 2;
    }  
  }

  // CSC ALCTs
  for(auto d: match_lct.chamberIdsALCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    auto alct = match_lct.alctInChamber(d);

    if (odd) etrk_[st].wiregroup_odd = digi_channel(alct);
    else etrk_[st].wiregroup_even = digi_channel(alct);

    if (odd) etrk_[st].quality_alct_odd = digi_quality(alct);
    else etrk_[st].quality_alct_even = digi_quality(alct);

    if (odd) etrk_[st].has_alct |= 1;
    else etrk_[st].has_alct |= 2;

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].wiregroup_odd = digi_channel(alct);
      else etrk_[1].wiregroup_even = digi_channel(alct);

      if (odd) etrk_[1].quality_alct_odd = digi_quality(alct);
      else etrk_[1].quality_alct_even = digi_quality(alct);
      
      if (odd) etrk_[1].has_alct |= 1;
      else etrk_[1].has_alct |= 2;      
    }
  }

  // holders for track's LCTs
  Digi lct_odd[12];
  Digi lct_even[12];
  GlobalPoint gp_lct_odd[12];
  GlobalPoint gp_lct_even[12];
  for (auto s: stations_to_use_)
  {
    lct_odd[s] = make_digi();
    lct_even[s] = make_digi();

    // case ME11
    if (s==2 or s==3){
      lct_odd[1] = make_digi();
      lct_even[1] = make_digi();
    }
  }

  // LCT stubs
  for(auto d: match_lct.chamberIdsLCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_lct |= 1;
    else etrk_[st].has_lct |= 2;

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_lct |= 1;
      else etrk_[1].has_lct |= 2;
    }
    
    auto lct = match_lct.lctInChamber(d);
    const int bend(LCT_BEND_PATTERN[digi_pattern(lct)]);
    auto gp = match_lct.digiPosition(lct);

    if (odd)
    {
      lct_odd[st] = lct;
      gp_lct_odd[st] = gp;
      etrk_[st].bend_lct_odd = bend;
      etrk_[st].phi_lct_odd = gp.phi();
      etrk_[st].eta_lct_odd = gp.eta();
      etrk_[st].perp_lct_odd = gp.perp();
      etrk_[st].dphi_lct_odd = digi_dphi(lct);
      etrk_[st].bx_lct_odd = digi_bx(lct);
      etrk_[st].hs_lct_odd = digi_channel(lct);
      etrk_[st].wg_lct_odd = digi_wg(lct);
      etrk_[st].chamber_odd |= 2;
      etrk_[st].quality_odd = digi_quality(lct);
      etrk_[st].passdphi_odd = match_lct.passDPhicut(id, chargesign, digi_dphi(lct), pt);
    }
    else
    {
      lct_even[st] = lct;
      gp_lct_even[st] = gp;
      etrk_[st].bend_lct_even = bend;
      etrk_[st].phi_lct_even = gp.phi();
      etrk_[st].eta_lct_even = gp.eta();
      etrk_[st].perp_lct_even = gp.perp();
      etrk_[st].dphi_lct_even = digi_dphi(lct);
      etrk_[st].bx_lct_even = digi_bx(lct);
      etrk_[st].hs_lct_even = digi_channel(lct);
      etrk_[st].wg_lct_even = digi_wg(lct);
      etrk_[st].chamber_even |= 2;
      etrk_[st].quality_even = digi_quality(lct);
      etrk_[st].passdphi_even = match_lct.passDPhicut(id, chargesign, digi_dphi(lct), pt);
    }

    // case ME11
    if (st==2 or st==3){
      if (odd)
      {
        lct_odd[1] = lct;
        gp_lct_odd[1] = gp;
        etrk_[1].bend_lct_odd = bend;
        etrk_[1].phi_lct_odd = gp.phi();
        etrk_[1].eta_lct_odd = gp.eta();
        etrk_[1].perp_lct_even = gp.perp();
        etrk_[1].dphi_lct_odd = digi_dphi(lct);
        etrk_[1].bx_lct_odd = digi_bx(lct);
        etrk_[1].hs_lct_odd = digi_channel(lct);
        etrk_[1].wg_lct_odd = digi_wg(lct);
        etrk_[1].chamber_odd |= 2;
        etrk_[1].quality_odd = digi_quality(lct);
        etrk_[1].passdphi_odd = match_lct.passDPhicut(id, chargesign, digi_dphi(lct), pt);
      }
      else
      {
        lct_even[1] = lct;
        gp_lct_even[1] = gp;
        etrk_[1].bend_lct_even = bend;
        etrk_[1].phi_lct_even = gp.phi();
        etrk_[1].eta_lct_even = gp.eta();
        etrk_[1].perp_lct_even = gp.perp();
        etrk_[1].dphi_lct_even = digi_dphi(lct);
        etrk_[1].bx_lct_even = digi_bx(lct);
        etrk_[1].hs_lct_even = digi_channel(lct);
        etrk_[1].wg_lct_even = digi_wg(lct);
        etrk_[1].chamber_even |= 2;
        etrk_[1].quality_even = digi_quality(lct);
        etrk_[1].passdphi_even = match_lct.passDPhicut(id, chargesign, digi_dphi(lct), pt);
      }

    }
    //find digis that are associatd to LCTs
      std::vector<GlobalPoint> gps;
      std::vector<float> phis;
      std::vector<float> xs;
      std::vector<float> ys;
      std::vector<float> zs;
      std::vector<float> ephis;
      std::vector<float> exs;
      std::vector<float> ezs;
      std::vector<float> status;
      
      auto bestMatchingLCT = match_lct.bestCscLctInChamber(d);
      match_lct.positionsOfComparatorInLCT(d, bestMatchingLCT, gps);
      if (gps.size()>=3){
	  for (auto gp: gps){
	        if (gp.z() > 0)
          		zs.push_back(gp.z());
		else zs.push_back(-gp.z());
		xs.push_back(gp.x());
		ys.push_back(gp.y());
          	ezs.push_back(0);
          	phis.push_back(gp.phi());
		ephis.push_back(gemvalidation::cscHalfStripWidth(id)/sqrt(12));
		float R=0.0;
		if (id.ring() == 1 or id.ring() == 4) R=200;//cm
		if (id.ring() == 2) R=400;//cm
		exs.push_back(gemvalidation::cscHalfStripWidth(id)/sqrt(12)*R);
	  }
      }else {
      	if (verbose_) std::cout <<" the size of gloabl points in this chamber is less than 3 "<< std::endl;
      }
     

      float alpha = 0., beta = 0.;
      calculateAlphaBeta(zs, phis, ezs, ephis, status,
                         alpha, beta);
      float alphax = 0., betax = 0.;
      calculateAlphaBeta(zs, xs, ezs, xs, status,
                         alphax, betax);
      float alphay = 0., betay = 0.;
      calculateAlphaBeta(zs, ys, ezs, xs, status,
                         alphay, betay);
      //ME11. even layer1 581.98; odd layer1 611.38, d=2.2
      //ME12. even layer1 678.706, layer6 691.406; odd layer1 706.106, layer6 718.806, d=2.54
      //ME21. even layer1 809.506; odd layer1 834.306 , d=2.54
      //ME22, even layer1 809.506, layer6 822.206; odd layer1 834.306, layer6 847.006, d=2.54
      //Z(layern) = Z(layer1)+d*(n-1)
      if(odd and gps.size()){
      	if (id.station()==1 and (id.ring()==1 or id.ring()==4)){
      		etrk_[st].phi_layer1_fit_odd = alpha+beta*611.38;
      		etrk_[st].phi_layer3_fit_odd = alpha+beta*(611.38+2.2*2);
		etrk_[st].phi_layer6_fit_odd = alpha+beta*(611.38+2.2*5);
		etrk_[st].phiM_fitxyz_odd = atan(betay/betax);
      		etrk_[1].phi_layer1_fit_odd = alpha+beta*611.38;
      		etrk_[1].phi_layer3_fit_odd = alpha+beta*(611.38+2.2*2);
		etrk_[1].phi_layer6_fit_odd = alpha+beta*(611.38+2.2*5);
		etrk_[1].phiM_fitxyz_odd = atan(betay/betax);
	}	
      	if (id.station()==1 and id.ring()==2){
      		etrk_[st].phi_layer1_fit_odd = alpha+beta*706.106;
      		etrk_[st].phi_layer3_fit_odd = alpha+beta*(706.106+2.54*2);
		etrk_[st].phi_layer6_fit_odd = alpha+beta*718.806;
		etrk_[st].phiM_fitxyz_odd = atan(betay/betax);
	}
      	if (id.station()==2){
      		etrk_[st].phi_layer1_fit_odd = alpha+beta*834.306;
      		etrk_[st].phi_layer3_fit_odd = alpha+beta*(834.306+2.54*2);
		etrk_[st].phi_layer6_fit_odd = alpha+beta*847.006;
		etrk_[st].phiM_fitxyz_odd = atan(betay/betax);
	}
	//std::cout <<"cscid "<< id <<" alpha "<< alpha <<" beta "<< beta <<" phi layer1 "<< etrk_[st].phi_layer1_fit_odd<<" layer6 "<<etrk_[st].phi_layer6_fit_odd<< std::endl;
      }else if (gps.size()){
      	if (id.station()==1 and (id.ring()==1 or id.ring(4))){
      		etrk_[st].phi_layer1_fit_even = alpha+beta*581.98;
      		etrk_[st].phi_layer3_fit_even = alpha+beta*(581.98+2.2*2);
		etrk_[st].phi_layer6_fit_even = alpha+beta*(581.98+2.2*5);
		etrk_[st].phiM_fitxyz_even = atan(betay/betax);
      		etrk_[1].phi_layer1_fit_even = alpha+beta*581.98;
      		etrk_[1].phi_layer3_fit_even = alpha+beta*(581.98+2.2*2);
		etrk_[1].phi_layer6_fit_even = alpha+beta*(581.98+2.2*5);
		etrk_[1].phiM_fitxyz_even = atan(betay/betax);
	}
      	if (id.station()==1 and id.ring()==2){
      		etrk_[st].phi_layer1_fit_even = alpha+beta*678.706;
      		etrk_[st].phi_layer3_fit_even = alpha+beta*(678.706+2.54*2);
		etrk_[st].phi_layer6_fit_even = alpha+beta*691.406;
		etrk_[st].phiM_fitxyz_even = atan(betay/betax);
	}
      	if (id.station()==2){
      		etrk_[st].phi_layer1_fit_even = alpha+beta*809.506;
      		etrk_[st].phi_layer3_fit_even = alpha+beta*(809.506+2.54*2);
		etrk_[st].phi_layer6_fit_even = alpha+beta*822.206;
		etrk_[st].phiM_fitxyz_even = atan(betay/betax);
	}
	//std::cout <<"cscid "<< id <<" alpha "<< alpha <<" beta "<< beta <<" phi layer1 "<< etrk_[st].phi_layer1_fit_even<<" layer6 "<<etrk_[st].phi_layer6_fit_even<< std::endl;
      }
      
  }
   

   //for GEMs in station1, it will be also filled in ME11
  // GEM simhits in superchamber
  for(auto d: match_sh.superChamberIdsGEM())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (match_sh.hitsInSuperChamber(d).size() > 0)
    {
      if (odd) etrk_[st].has_gem_sh |= 1;
      else     etrk_[st].has_gem_sh |= 2;
 
      for (int layer=1; layer<3; layer++){
	GEMDetId id_tmp(id.region(), id.ring(), id.station(), layer, id.chamber(), 0);
	GlobalPoint keygp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id_tmp.rawId()));
	GlobalVector keygv = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(id_tmp.rawId()));
	float centralperp_gemsh = match_sh.simHitsGEMCentralPosition(match_sh.hitsInChamber(id_tmp.rawId()));match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id_tmp.rawId()));
	if(match_sh.hitsInChamber(id_tmp).size()==0) continue;
	//std::cout <<" GEM Id "<< id <<" gp.eta "<< sh_gp.eta() <<" gp.phi "<< sh_gp.phi() << std::endl;
        //std::cout <<" GEM Id "<< id_tmp <<" eta "<< keygp.eta() <<" phi "<< keygp.phi()<<" perp "<< keygp.perp() <<" x " << keygp.x() <<" y "<< keygp.y() << " z "<< keygp.z()<< std::endl;
      	if (odd) etrk_[st].eta_gemsh_odd = keygp.eta();
     	else     etrk_[st].eta_gemsh_even = keygp.eta();
      	if (odd) etrk_[st].phi_gemsh_odd = keygp.phi();
      	else     etrk_[st].phi_gemsh_even = keygp.phi();
      	if (odd) etrk_[st].ptphi_gemsh_odd = keygv.phi();
      	else     etrk_[st].ptphi_gemsh_even = keygv.phi();
      	if (odd) etrk_[st].perp_gemsh_odd = keygp.perp();
      	else     etrk_[st].perp_gemsh_even = keygp.perp();
	if (odd) etrk_[st].centralperp_gemsh_odd = centralperp_gemsh;
	else     etrk_[st].centralperp_gemsh_even = centralperp_gemsh;
	float phi_m = -9;
	float phi_m2 = -9;
	int st_csc = (st==3 and (etrk_[st].has_csc_sh&3)==0)? 1:st;
      	if (odd and etrk_[st_csc].phi_cscsh_odd>-9){
	    etrk_[st].dphi_sh_odd = deltaPhi(etrk_[st_csc].phi_cscsh_odd, keygp.phi());// here should we use st or st_csc?
    	    etrk_[st].dphipositionpt_gemsh_odd = deltaPhi(keygp.phi(),etrk_[st_csc].ptphi_sh_odd);
	    phi_m =  PhiMomentum_Radius(etrk_[st].dphi_sh_odd, etrk_[st].phi_gemsh_odd, etrk_[st_csc].perp_cscsh_odd, etrk_[st].perp_gemsh_odd);
	    phi_m2 =  PhiMomentum_Radius(etrk_[st].dphi_sh_odd, etrk_[st].phi_gemsh_odd, etrk_[st_csc].perp_cscsh_odd, etrk_[st].centralperp_gemsh_odd);
	    etrk_[st].phiM_gemcsc_odd = phi_m;
	    etrk_[st].phiM_gemcsc_central_odd = phi_m2;
	    //std::cout <<"dphi(M,P) "<< deltaPhi(etrk_[st].phi_gemsh_odd, gv_sh_odd[st_csc].phi())<<" momentum phi from simhits "<< gv_sh_odd[st_csc].phi() <<" phi from gem-csc "<< phi_m <<" central "<< phi_m2<< std::endl;
	}
      	else if (not(odd) and etrk_[st_csc].phi_cscsh_even>-9) {
	    etrk_[st].dphi_sh_even = deltaPhi(etrk_[st_csc].phi_cscsh_even, keygp.phi());
    	    etrk_[st].dphipositionpt_gemsh_even = deltaPhi(keygp.phi(),etrk_[st_csc].ptphi_sh_even);
	    phi_m =  PhiMomentum_Radius(etrk_[st].dphi_sh_even, etrk_[st].phi_gemsh_even, etrk_[st_csc].perp_cscsh_even, etrk_[st].perp_gemsh_even);
	    phi_m2 =  PhiMomentum_Radius(etrk_[st].dphi_sh_even, etrk_[st].phi_gemsh_even, etrk_[st_csc].perp_cscsh_even, etrk_[st].centralperp_gemsh_even);
	    etrk_[st].phiM_gemcsc_even = phi_m;
	    etrk_[st].phiM_gemcsc_central_even = phi_m2;
	    //std::cout <<"dphi(M,P) "<< deltaPhi(etrk_[st].phi_gemsh_even, gv_sh_even[st_csc].phi())<<" momentum phi from simhits "<< gv_sh_even[st_csc].phi() <<" phi from gem-csc "<< phi_m <<" central "<< phi_m2<< std::endl;
	}
	if (st==2 or st==3){
      		if (odd) etrk_[1].eta_gemsh_odd = keygp.eta();
      		else     etrk_[1].eta_gemsh_even = keygp.eta();
      		if (odd) etrk_[1].phi_gemsh_odd = keygp.phi();
      		else     etrk_[1].phi_gemsh_even = keygp.phi();
      		if (odd) etrk_[1].ptphi_gemsh_odd = keygv.phi();
      		else     etrk_[1].ptphi_gemsh_even = keygv.phi();
      		if (odd) etrk_[1].perp_gemsh_odd = keygp.perp();
      		else     etrk_[1].perp_gemsh_even = keygp.perp();
		if (odd) etrk_[1].centralperp_gemsh_odd = centralperp_gemsh;
		else     etrk_[1].centralperp_gemsh_even = centralperp_gemsh;
      		if (odd and etrk_[1].phi_cscsh_odd>-9) etrk_[1].dphi_sh_odd = deltaPhi(etrk_[1].phi_cscsh_odd,keygp.phi());
      		else if (etrk_[1].phi_cscsh_even>-9)     etrk_[1].dphi_sh_even = deltaPhi(etrk_[1].phi_cscsh_even,keygp.phi());
		if (odd) etrk_[1].phiM_gemcsc_odd = phi_m;
		else etrk_[1].phiM_gemcsc_even = phi_m;
		if (odd) etrk_[1].phiM_gemcsc_central_odd = phi_m2;
		else etrk_[1].phiM_gemcsc_central_even = phi_m2;
		if (odd) etrk_[1].dphipositionpt_gemsh_odd = deltaPhi(keygp.phi(),etrk_[1].ptphi_sh_odd);
		else etrk_[1].dphipositionpt_gemsh_even = deltaPhi(keygp.phi(),etrk_[1].ptphi_sh_even);
	}
	//std::cout <<"dphi(GEM, CSC) at sim "<<(odd? etrk_[st].dphi_sh_odd:etrk_[st].dphi_sh_even) <<" phi_csc "<< (odd? etrk_[st_csc].phi_cscsh_odd:etrk_[st_csc].phi_cscsh_even) <<" phi_gp "<< keygp.phi()<< std::endl;
	if (id_tmp.layer()==1) break;
  		
      }

      const float mean_strip(match_sh.simHitsMeanStrip(match_sh.hitsInSuperChamber(d)));
      if (odd) etrk_[st].strip_gemsh_odd = mean_strip;
      else     etrk_[st].strip_gemsh_even = mean_strip;
    }

    if (match_sh.nLayersWithHitsInSuperChamber(d) > 1)
    {
      if (odd) etrk_[st].has_gem_sh2 |= 1;
      else     etrk_[st].has_gem_sh2 |= 2;
    }
    //ME11 Case
    if (st==2 or st==3)
    {
      if (odd) etrk_[1].has_gem_sh |= 1;
      else     etrk_[1].has_gem_sh |= 2;

      const float mean_strip(match_sh.simHitsMeanStrip(match_sh.hitsInSuperChamber(d)));
      if (odd) etrk_[1].strip_gemsh_odd = mean_strip;
      else     etrk_[1].strip_gemsh_even = mean_strip;
      
    if (match_sh.nLayersWithHitsInSuperChamber(d) > 1)
    {
      if (odd) etrk_[1].has_gem_sh2 |= 1;
      else etrk_[1].has_gem_sh2 |= 2;

    }
  }//end of ME11 case

  }

  // placeholders for best mtching pads
  GlobalPoint best_pad_odd[12];
  GlobalPoint best_pad_even[12];

  // GEM digis and pads in superchambers
  for(auto d: match_gd.superChamberIdsDigi())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (match_gd.nLayersWithDigisInSuperChamber(d) > 1)
    {
      if (odd) etrk_[st].has_gem_dg2 |= 1;
      else     etrk_[st].has_gem_dg2 |= 2;
    }

    auto digis = match_gd.digisInSuperChamber(d);
    const int median_strip(match_gd.median(digis));
    if (odd && digis.size() > 0)
    {
      etrk_[st].has_gem_dg |= 1;
      etrk_[st].strip_gemdg_odd = median_strip;
    }
    else if (digis.size() > 0)
    {
      etrk_[st].has_gem_dg |= 2;
      etrk_[st].strip_gemdg_even = median_strip;
    }

    if (match_gd.nLayersWithPadsInSuperChamber(d) > 1)
    {
      if (odd) etrk_[st].has_gem_pad2 |= 1;
      else     etrk_[st].has_gem_pad2 |= 2;
    }
    for (int layer=1; layer<3; layer++){
      GEMDetId id_tmp(id.region(), id.ring(), id.station(), layer, id.chamber(), 0);
      auto pads = match_gd.padsInChamber(id_tmp.rawId());
      if(pads.size() == 0) continue;
      if (odd)
      {
        etrk_[st].has_gem_pad |= 1;
        etrk_[st].chamber_odd |= 1;
        etrk_[st].pad_odd = digi_channel(pads.at(0));
        etrk_[st].hsfromgem_odd = match_gd.extrapolateHsfromGEMPad( d, digi_channel(pads.at(0)));
        if (is_valid(lct_odd[st]))
        {
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_odd[st]);
        	best_pad_odd[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_odd = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].phi_pad_odd = best_pad_odd[st].phi();
        	etrk_[st].eta_pad_odd = best_pad_odd[st].eta();
        	etrk_[st].dphi_pad_odd = deltaPhi(etrk_[st].phi_lct_odd, etrk_[st].phi_pad_odd);
        	etrk_[st].deta_pad_odd = etrk_[st].eta_lct_odd - etrk_[st].eta_pad_odd;
      	}
    	}
      else
      {
      	etrk_[st].has_gem_pad |= 2;
        etrk_[st].chamber_even |= 1;
        etrk_[st].pad_even = digi_channel(pads.at(0));
        etrk_[st].hsfromgem_even = match_gd.extrapolateHsfromGEMPad( d, digi_channel(pads.at(0)));
        if (is_valid(lct_even[st]))
        {
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_even[st]);
        	best_pad_even[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_even = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].phi_pad_even = best_pad_even[st].phi();
        	etrk_[st].eta_pad_even = best_pad_even[st].eta();
        	etrk_[st].dphi_pad_even = deltaPhi(etrk_[st].phi_lct_even, etrk_[st].phi_pad_even);
        	etrk_[st].deta_pad_even = etrk_[st].eta_lct_even - etrk_[st].eta_pad_even;
      	}
    	}
      if (id_tmp.layer()==1) break;
    }
  }

  //ME11Case
  for(auto d: match_gd.superChamberIdsDigi())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();

    const int stations(detIdToMEStation(MEStation,id.ring()));
    int st;
    if (stations==2 or stations==3) st=1;
    else continue;
    
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (match_gd.nLayersWithDigisInSuperChamber(d) > 1)
    {
      if (odd) etrk_[st].has_gem_dg2 |= 1;
      else     etrk_[st].has_gem_dg2 |= 2;
    }

    auto digis = match_gd.digisInSuperChamber(d);
    const int median_strip(match_gd.median(digis));
    if (odd && digis.size() > 0)
    {
      etrk_[st].has_gem_dg |= 1;
      etrk_[st].strip_gemdg_odd = median_strip;
    }
    else if (digis.size() > 0)
    {
      etrk_[st].has_gem_dg |= 2;
      etrk_[st].strip_gemdg_even = median_strip;
    }

    if (match_gd.nLayersWithPadsInSuperChamber(d) > 1)
    {
      if (odd) etrk_[st].has_gem_pad2 |= 1;
      else     etrk_[st].has_gem_pad2 |= 2;
    }

    for (int layer=1; layer<3; layer++){
      GEMDetId id_tmp(id.region(), id.ring(), id.station(), layer, id.chamber(), 0);
      auto pads = match_gd.padsInChamber(id_tmp.rawId());
      if(pads.size() == 0) continue;
      
      if (odd)
      {
      	etrk_[st].has_gem_pad |= 1;
      	etrk_[st].chamber_odd |= 1;
      	etrk_[st].pad_odd = digi_channel(pads.at(0));
      	if (is_valid(lct_odd[st]))
      	{
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_odd[st]);
        	best_pad_odd[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_odd = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].phi_pad_odd = best_pad_odd[st].phi();
        	etrk_[st].eta_pad_odd = best_pad_odd[st].eta();
        	etrk_[st].dphi_pad_odd = deltaPhi(etrk_[st].phi_lct_odd, etrk_[st].phi_pad_odd);
        	etrk_[st].deta_pad_odd = etrk_[st].eta_lct_odd - etrk_[st].eta_pad_odd;
      	}
    	}
      else
      {
      	etrk_[st].has_gem_pad |= 2;
      	etrk_[st].chamber_even |= 1;
      	etrk_[st].pad_even = digi_channel(pads.at(0));
      	if (is_valid(lct_even[st]))
      	{
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_even[st]);
        	best_pad_even[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_even = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].phi_pad_even = best_pad_even[st].phi();
        	etrk_[st].eta_pad_even = best_pad_even[st].eta();
        	etrk_[st].dphi_pad_even = deltaPhi(etrk_[st].phi_lct_even, etrk_[st].phi_pad_even);
        	etrk_[st].deta_pad_even = etrk_[st].eta_lct_even - etrk_[st].eta_pad_even;
      	}
    	}
      if (id_tmp.layer()==1) break;
    }
   }

  for(auto d: match_gd.superChamberIdsCoPad())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_gem_copad |= 1;
    else     etrk_[st].has_gem_copad |= 2;
    
    auto copads = match_gd.coPadsInSuperChamber(d);
    if (copads.size() == 0) continue;
    if (odd) etrk_[st].Copad_odd = digi_channel(copads.at(0));
    else etrk_[st].Copad_even = digi_channel(copads.at(0));

    if (st==2 or st==3)
    {
    if (odd) etrk_[1].has_gem_copad |= 1;
    else     etrk_[1].has_gem_copad |= 2;
    
    auto copads = match_gd.coPadsInSuperChamber(d);
    if (copads.size() == 0) continue;
    if (odd) etrk_[1].Copad_odd = digi_channel(copads.at(0));
    else etrk_[1].Copad_even = digi_channel(copads.at(0));
    }
  }
 
  // placeholders for best mtching rpcstrips
  GlobalPoint best_rpcstrip_odd[12];
  GlobalPoint best_rpcstrip_even[12];

  if (false) for (auto d: match_sh.chamberIdsRPC())
  {
    RPCDetId id(d);
    const int st(detIdToMEStation(id.station(), id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int cscchamber = CSCTriggerNumbering::chamberFromTriggerLabels(id.sector(), 0, id.station(), id.subsector());
    cscchamber = (cscchamber+16)%18+1; 
    if ( (match_sh.hitsInChamber(d)).size() >0 )
    {
      bool odd(cscchamber%2 == 1);
      if (odd)   etrk_[st].has_rpc_sh |= 1;
      else etrk_[st].has_rpc_sh |=2;  
    }	
  }

  if (false) for (auto d: match_rd.detIds())
  {
    RPCDetId id(d);
    const int st(detIdToMEStation(id.station(), id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    //meanstrip in rpc 
    auto rpcdigis = match_rd.digisInDetId(id); 
    const int rpc_medianstrip(match_rd.median(rpcdigis));
    const int cscchamber = CSCTriggerNumbering::chamberFromTriggerLabels(id.sector(), 0, id.station(), id.subsector());
    //std::cout <<"rpc detid " << id << " csc chamebr:"<< cscchamber << std::endl;
    const bool odd(cscchamber%2 == 1);
    if (odd)
    {
      etrk_[st].has_rpc_dg |= 1;
//       etrk_[st].chamber_odd |= 3;
      etrk_[st].strip_rpcdg_odd = rpc_medianstrip;
      etrk_[st].hsfromrpc_odd = match_rd.extrapolateHsfromRPC( d, rpc_medianstrip);
      if (is_valid(lct_odd[st]))
      {
        auto rpc_dg_and_gp = match_gd.digiInRPCClosestToCSC(rpcdigis, gp_lct_odd[st]);
        best_rpcstrip_odd[st] = rpc_dg_and_gp.second;
        etrk_[st].bx_rpcstrip_odd = digi_bx(rpc_dg_and_gp.first);
        etrk_[st].phi_rpcstrip_odd = best_rpcstrip_odd[st].phi();
        etrk_[st].eta_rpcstrip_odd = best_rpcstrip_odd[st].eta();
        etrk_[st].dphi_rpcstrip_odd = deltaPhi(etrk_[st].phi_lct_odd, etrk_[st].phi_rpcstrip_odd);
        etrk_[st].deta_rpcstrip_odd = etrk_[st].eta_lct_odd - etrk_[st].eta_rpcstrip_odd;
      }
    }
    else
    {
      etrk_[st].has_rpc_dg |= 2;
//       etrk_[st].chamber_even |= 3;
      etrk_[st].strip_rpcdg_even = rpc_medianstrip;
      etrk_[st].hsfromrpc_even = match_rd.extrapolateHsfromRPC( d, rpc_medianstrip);
      if (is_valid(lct_even[st]))
      {
        auto rpc_dg_and_gp = match_gd.digiInRPCClosestToCSC(rpcdigis, gp_lct_even[st]);
        best_rpcstrip_even[st] = rpc_dg_and_gp.second;
        etrk_[st].bx_rpcstrip_even = digi_bx(rpc_dg_and_gp.first);
        etrk_[st].phi_rpcstrip_even = best_rpcstrip_even[st].phi();
        etrk_[st].eta_rpcstrip_even = best_rpcstrip_even[st].eta();
        etrk_[st].dphi_rpcstrip_even = deltaPhi(etrk_[st].phi_lct_even, etrk_[st].phi_rpcstrip_even);
        etrk_[st].deta_rpcstrip_even = etrk_[st].eta_lct_even - etrk_[st].eta_rpcstrip_even;
      }
    }
  }



  GlobalPoint gp1,gp2, gp3;
  GlobalVector gv1,gv2;
  int npar=-1;
  float csc_bending_angle_12=-9;
  float csc_bending_angle12_gemcsc = -9;
  float csc_bending_angle12_gemcsc_central = -9;
  float csc_bending_angle12_xfactor = -9;
  float csc_bending_angle12_xfactor_smear1 = -9;
  float csc_bending_angle12_xfactor_smear2 = -9;
  float csc_bending_angle12_xfactor_smear3 = -9;
  float csc_bending_angle12_xfactor_smear0 = -9;
  float csc_bending_angle12_xfactor_L1_1 = -9;
  float csc_bending_angle12_xfactor_L1_2 = -9;
  //float dphiME11 = -9;
  //float dphiME21 = -9;
  //float dphiME11_gauss = -9;
  //float dphiME21_gauss = -9;
  float gp1_phi_gauss = -9;
  float gp2_phi_gauss = -9;
  float	gem1_phi_gauss = -9;
  float gem2_phi_gauss = -9;

//x_factor = (R(ME21)/R(ME11)-1)/D(ME11, ME21)
//X: st1, X=D(GEM,CSC)*x_factor, st2: D(GEM,CSC)*x_factor/(D(ME11,ME21)*x+1)
  float xfactor = 0.0;
  float xfactor_st1 = 0.0;
  float xfactor_st2 = 0.0;
  //float xfactor_smear = 0.0;
  float phiM_st1_sh = 0.0;
  float phiM_st2_sh = 0.0;
  float phiM_st1_smear0 = 0.0;
  float phiM_st2_smear0 = 0.0;
  float phiM_st1_smear1 = 0.0;
  float phiM_st2_smear1 = 0.0;
  float phiM_st1_smear2 = 0.0;
  float phiM_st2_smear2 = 0.0;
  float phiM_st1_smear3 = 0.0;
  float phiM_st2_smear3 = 0.0;
  float phiM_st1_L1_1 = 0.0;
  float phiM_st2_L1_1 = 0.0;
  float phiM_st1_L1_2 = 0.0;
  float phiM_st2_L1_2 = 0.0;


   
  if (etrk_[1].has_csc_sh>0 and etrk_[6].has_csc_sh>0){
     if ((etrk_[1].has_csc_sh&1)>0 and (etrk_[6].has_csc_sh&2)>0){
        gp1=gp_sh_odd[1];
        gp2=gp_sh_even[6];
        gv1=gv_sh_odd[1];
        gv2=gv_sh_even[6];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
        gem1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_gemsh_odd, .00021);
        gem2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_gemsh_even, .00014);
        gp1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_cscsh_odd, .00026);
        gp2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_cscsh_even, .00045);
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(814.586-615.33);
	//xfactor_smear = CLHEP::RandGauss::shoot(xfactor,0.00005);
	//std::cout <<"gp1 perp "<< gp1.perp() <<" gp2 perp "<< gp2.perp() <<"xfactor "<< xfactor <<" smear "<< xfactor_smear << std::endl;
	xfactor_st1 = (615.33-567.82)*xfactor;
	xfactor_st2 = (814.586-794.29)*xfactor/(xfactor*(814.586-615.33)+1);
	npar=0;
	if ((etrk_[8].has_csc_sh&2)>0)
	    gp3=gp_sh_even[8];
	if (etrk_[1].dphi_sh_odd>-9 and etrk_[6].dphi_sh_even>-9)
	{
	    //dphiME11 = etrk_[1].dphi_sh_odd;
	    //dphiME21 = etrk_[6].dphi_sh_even;
    	    //dphiME11_gauss = CLHEP::RandGauss::shoot(dphiME11, .0005);
     	    //dphiME21_gauss = CLHEP::RandGauss::shoot(dphiME21, .0005);
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[1].phi_cscsh_odd, etrk_[1].phi_gemsh_odd, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear0 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, etrk_[1].phi_gemsh_odd, xfactor_st1);
	    phiM_st2_smear0 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear1 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear2 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear2 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, gem2_phi_gauss, xfactor_st2);

	    phiM_st1_L1_1 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_odd, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_2 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_odd, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_2 = PhiMomentum_Xfactor_V2(etrk_[6].phi_layer3_fit_even, gem2_phi_gauss, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear0 = deltaPhi(phiM_st1_smear0,phiM_st2_smear0);
		csc_bending_angle12_xfactor_smear1 = deltaPhi(phiM_st1_smear1,phiM_st2_smear1);
		csc_bending_angle12_xfactor_smear2 = deltaPhi(phiM_st1_smear2,phiM_st2_smear2);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
	    }
	    if (phiM_st1_L1_2>-9 and phiM_st2_L1_2>-9){
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(phiM_st1_L1_1,phiM_st2_L1_1);
		csc_bending_angle12_xfactor_L1_2 = deltaPhi(phiM_st1_L1_2,phiM_st2_L1_2);
	    }

	    if (etrk_[1].phiM_gemcsc_odd>-9 and etrk_[6].phiM_gemcsc_even>-9)
     	    	csc_bending_angle12_gemcsc = deltaPhi(etrk_[1].phiM_gemcsc_odd, etrk_[6].phiM_gemcsc_even);
	    if (etrk_[1].phiM_gemcsc_central_odd>-9 and etrk_[6].phiM_gemcsc_central_even>-9)
     	    	csc_bending_angle12_gemcsc_central = deltaPhi(etrk_[1].phiM_gemcsc_central_odd, etrk_[6].phiM_gemcsc_central_even);
	}
     }else if ((etrk_[1].has_csc_sh&1)>0 and (etrk_[6].has_csc_sh&1)>0 ){ 
        gp1=gp_sh_odd[1];
        gp2=gp_sh_odd[6];
        gv1=gv_sh_odd[1];
        gv2=gv_sh_odd[6];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
        gem1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_gemsh_odd, .00021);
        gem2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_gemsh_odd, .00014);
        gp1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_cscsh_odd, .00026);
        gp2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_cscsh_odd, .00045);
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(839.386-615.33);
        //xfactor_smear = CLHEP::RandGauss::shoot(xfactor,0.00005);
	//std::cout <<"gp1 perp "<< gp1.perp() <<" gp2 perp "<< gp2.perp() <<"xfactor "<< xfactor <<" smear "<< xfactor_smear << std::endl;
	xfactor_st1 = (615.33-567.82)*xfactor;
	xfactor_st2 = (839.386-796.792)*xfactor/(xfactor*(839.386-615.33)+1);
	npar=1;
	if ((etrk_[8].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[8];
	if (etrk_[1].dphi_sh_odd>-9 and etrk_[6].dphi_sh_odd>-9)
	{
	    //dphiME11 = etrk_[1].dphi_sh_odd;
	    //dphiME21 = etrk_[6].dphi_sh_odd;
    	    //dphiME11_gauss = CLHEP::RandGauss::shoot(dphiME11, .0005);
     	    //dphiME21_gauss = CLHEP::RandGauss::shoot(dphiME21, .0005);
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[1].phi_cscsh_odd, etrk_[1].phi_gemsh_odd, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear0 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, etrk_[1].phi_gemsh_odd, xfactor_st1);
	    phiM_st2_smear0 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear1 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear2 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear2 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_1 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_odd, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_2 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_odd, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_2 = PhiMomentum_Xfactor_V2(etrk_[6].phi_layer3_fit_odd, gem2_phi_gauss, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear0 = deltaPhi(phiM_st1_smear0,phiM_st2_smear0);
		csc_bending_angle12_xfactor_smear1 = deltaPhi(phiM_st1_smear1,phiM_st2_smear1);
		csc_bending_angle12_xfactor_smear2 = deltaPhi(phiM_st1_smear2,phiM_st2_smear2);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
	    }
	    if (phiM_st1_L1_2>-9 and phiM_st2_L1_2>-9){
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(phiM_st1_L1_1,phiM_st2_L1_1);
		csc_bending_angle12_xfactor_L1_2 = deltaPhi(phiM_st1_L1_2,phiM_st2_L1_2);
	    }

	    if (etrk_[1].phiM_gemcsc_odd>-9 and etrk_[6].phiM_gemcsc_odd>-9)
     	    	csc_bending_angle12_gemcsc = deltaPhi(etrk_[1].phiM_gemcsc_odd, etrk_[6].phiM_gemcsc_odd);
	    if (etrk_[1].phiM_gemcsc_central_odd>-9 and etrk_[6].phiM_gemcsc_central_odd>-9)
     	    	csc_bending_angle12_gemcsc_central = deltaPhi(etrk_[1].phiM_gemcsc_central_odd, etrk_[6].phiM_gemcsc_central_odd);
	}
     }else if ((etrk_[1].has_csc_sh&2)>0 and (etrk_[6].has_csc_sh&2)>0 ){ 
        gp1=gp_sh_even[1];
        gp2=gp_sh_even[6];
        gv1=gv_sh_even[1];
        gv2=gv_sh_even[6];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
        gem1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_gemsh_even, .00021);
        gem2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_gemsh_even, .00014);
        gp1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_cscsh_even, .00026);
        gp2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_cscsh_even, .00045);
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(814.586-586.03);
        //xfactor_smear = CLHEP::RandGauss::shoot(xfactor,0.00005);
	//std::cout <<"gp1 perp "<< gp1.perp() <<" gp2 perp "<< gp2.perp() <<"xfactor "<< xfactor <<" smear "<< xfactor_smear << std::endl;
	xfactor_st1 = (586.03-565.32)*xfactor;
	xfactor_st2 = (814.586-794.29)*xfactor/(xfactor*(814.586-586.03)+1);
	npar=2;
	if ((etrk_[8].has_csc_sh&2)>0)
	    gp3=gp_sh_even[8];
	if (etrk_[1].dphi_sh_even>-9 and etrk_[6].dphi_sh_even>-9)
	{
	    //dphiME11 = etrk_[1].dphi_sh_even;
	    //dphiME21 = etrk_[6].dphi_sh_even;
    	    //dphiME11_gauss = CLHEP::RandGauss::shoot(dphiME11, .0005);
     	    //dphiME21_gauss = CLHEP::RandGauss::shoot(dphiME21, .0005);
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[1].phi_cscsh_even, etrk_[1].phi_gemsh_even, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear0 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, etrk_[1].phi_gemsh_even, xfactor_st1);
	    phiM_st2_smear0 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear1 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear2 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear2 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, etrk_[6].phi_gemsh_even, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_1 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_even, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_even, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_2 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_even, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_2 = PhiMomentum_Xfactor_V2(etrk_[6].phi_layer3_fit_even, gem2_phi_gauss, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear0 = deltaPhi(phiM_st1_smear0,phiM_st2_smear0);
		csc_bending_angle12_xfactor_smear1 = deltaPhi(phiM_st1_smear1,phiM_st2_smear1);
		csc_bending_angle12_xfactor_smear2 = deltaPhi(phiM_st1_smear2,phiM_st2_smear2);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
	    }
	    if (phiM_st1_L1_2>-9 and phiM_st2_L1_2>-9){
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(phiM_st1_L1_1,phiM_st2_L1_1);
		csc_bending_angle12_xfactor_L1_2 = deltaPhi(phiM_st1_L1_2,phiM_st2_L1_2);
	    }

	    if (etrk_[1].phiM_gemcsc_even>-9 and etrk_[6].phiM_gemcsc_even>-9)
     	    	csc_bending_angle12_gemcsc = deltaPhi(etrk_[1].phiM_gemcsc_even, etrk_[6].phiM_gemcsc_even);
	    if (etrk_[1].phiM_gemcsc_central_even>-9 and etrk_[6].phiM_gemcsc_central_even>-9)
     	    	csc_bending_angle12_gemcsc_central = deltaPhi(etrk_[1].phiM_gemcsc_central_even, etrk_[6].phiM_gemcsc_central_even);
	}
     }else if ((etrk_[1].has_csc_sh&2)>0 and (etrk_[6].has_csc_sh&1)>0 ){ 
        gp1=gp_sh_even[1];
        gp2=gp_sh_odd[6];
        gv1=gv_sh_even[1];
        gv2=gv_sh_odd[6];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
        gem1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_gemsh_even, .00021);
        gem2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_gemsh_odd, .00014);
        gp1_phi_gauss = CLHEP::RandGauss::shoot(etrk_[1].phi_cscsh_even, .00026);
        gp2_phi_gauss = CLHEP::RandGauss::shoot(etrk_[6].phi_cscsh_odd, .00045);
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(839.386-586.03);
        //xfactor_smear = CLHEP::RandGauss::shoot(xfactor,0.00005);
	xfactor_st1 = (586.03-565.32)*xfactor;
	xfactor_st2 = (839.386-796.792)*xfactor/(xfactor*(839.386-586.03)+1);
	//std::cout <<"gp1 perp "<< gp1.perp() <<" gp2 perp "<< gp2.perp() <<"xfactor "<< xfactor <<" smear "<< xfactor_smear <<" xfactor_st1 "<< xfactor_st1 <<" xfactor_st2 "<< xfactor_st2 << std::endl;
	npar=3;
	if ((etrk_[8].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[8];
	if (etrk_[1].dphi_sh_even>-9 and etrk_[6].dphi_sh_odd>-9)
	{
	    //dphiME11 = etrk_[1].dphi_sh_even;
	    //dphiME21 = etrk_[6].dphi_sh_odd;
    	    //dphiME11_gauss = CLHEP::RandGauss::shoot(dphiME11, .0005);//0.00055
     	    //dphiME21_gauss = CLHEP::RandGauss::shoot(dphiME21, .0005);//0.00095
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[1].phi_cscsh_even, etrk_[1].phi_gemsh_even, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear0 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, etrk_[1].phi_gemsh_even, xfactor_st1);
	    phiM_st2_smear0 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear1 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear2 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear2 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, etrk_[6].phi_gemsh_odd, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(gp1_phi_gauss, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(gp2_phi_gauss, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_1 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_even, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_1 = PhiMomentum_Xfactor_V2(etrk_[6].phi_cscsh_odd, gem2_phi_gauss, xfactor_st2);
	    phiM_st1_L1_2 = PhiMomentum_Xfactor_V2(etrk_[1].phi_layer3_fit_even, gem1_phi_gauss, xfactor_st1);
	    phiM_st2_L1_2 = PhiMomentum_Xfactor_V2(etrk_[6].phi_layer3_fit_odd, gem2_phi_gauss, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear0 = deltaPhi(phiM_st1_smear0,phiM_st2_smear0);
		csc_bending_angle12_xfactor_smear1 = deltaPhi(phiM_st1_smear1,phiM_st2_smear1);
		csc_bending_angle12_xfactor_smear2 = deltaPhi(phiM_st1_smear2,phiM_st2_smear2);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
	    }
	    if (phiM_st1_L1_2>-9 and phiM_st2_L1_2>-9){
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(phiM_st1_L1_1,phiM_st2_L1_1);
		csc_bending_angle12_xfactor_L1_2 = deltaPhi(phiM_st1_L1_2,phiM_st2_L1_2);
	    }

	    if (etrk_[1].phiM_gemcsc_even>-9 and etrk_[6].phiM_gemcsc_odd>-9)
     	    	csc_bending_angle12_gemcsc = deltaPhi(etrk_[1].phiM_gemcsc_even, etrk_[6].phiM_gemcsc_odd);
	    if (etrk_[1].phiM_gemcsc_central_even>-9 and etrk_[6].phiM_gemcsc_central_odd>-9)
     	    	csc_bending_angle12_gemcsc_central = deltaPhi(etrk_[1].phiM_gemcsc_central_even, etrk_[6].phiM_gemcsc_central_odd);
	}
     }
     //std::cout <<"high eta region npar "<< npar<<" eta_st2 "<< gp2.eta() <<" csc bending angle12 "<< csc_bending_angle_12<<" xfactor1 "<< xfactor_st1<<" xfactor2 "<< xfactor_st2 <<" csc_bending_angle12_xfactor "<< csc_bending_angle12_xfactor << std::endl;
  } 

     //low eta region , ME12, ME22
  if (etrk_[4].has_csc_sh>0 and etrk_[7].has_csc_sh>0){
     if ((etrk_[4].has_csc_sh&1)>0 and (etrk_[7].has_csc_sh&2)>0){
        gp1=gp_sh_odd[4];
        gp2=gp_sh_even[7];
        gv1=gv_sh_odd[4];
        gv2=gv_sh_even[7];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(814.586-711.186);
	xfactor_st1 = (12.7)*xfactor;
	xfactor_st2 = (12.7)*xfactor/(xfactor*(814.586-711.186)+1);
	npar=0;
	if (etrk_[4].phi_layer1_sh_odd>-9  and etrk_[4].phi_layer6_sh_odd>-9 and etrk_[7].phi_layer1_sh_even>-9  and etrk_[7].phi_layer6_sh_even>-9){
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_sh_odd, etrk_[4].phi_layer1_sh_odd, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_sh_even, etrk_[7].phi_layer1_sh_even, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_fit_odd, etrk_[4].phi_layer1_fit_odd, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_fit_even, etrk_[7].phi_layer1_fit_even, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(etrk_[4].phiM_fitxyz_odd, etrk_[7].phiM_fitxyz_even);
		}
	}
     }else if ((etrk_[4].has_csc_sh&1)>0 and (etrk_[7].has_csc_sh&1)>0){
        gp1=gp_sh_odd[4];
        gp2=gp_sh_odd[7];
        gv1=gv_sh_odd[4];
        gv2=gv_sh_odd[7];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(839.386-711.186);
	xfactor_st1 = (12.7)*xfactor;
	xfactor_st2 = (12.7)*xfactor/(xfactor*(839.386-711.186)+1);
	npar=1;
	if (etrk_[4].phi_layer1_sh_odd>-9  and etrk_[4].phi_layer6_sh_odd>-9 and etrk_[7].phi_layer1_sh_odd>-9  and etrk_[7].phi_layer6_sh_odd>-9){
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_sh_odd, etrk_[4].phi_layer1_sh_odd, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_sh_odd, etrk_[7].phi_layer1_sh_odd, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_fit_odd, etrk_[4].phi_layer1_fit_odd, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_fit_odd, etrk_[7].phi_layer1_fit_odd, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(etrk_[4].phiM_fitxyz_odd, etrk_[7].phiM_fitxyz_odd);
		}
	}
     }else if ((etrk_[4].has_csc_sh&2)>0 and (etrk_[7].has_csc_sh&2)>0){
        gp1=gp_sh_even[4];
        gp2=gp_sh_even[7];
        gv1=gv_sh_even[4];
        gv2=gv_sh_even[7];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(814.586-683.786);
	xfactor_st1 = (12.7)*xfactor;
	xfactor_st2 = (12.7)*xfactor/(xfactor*(814.586-683.786)+1);
	npar=2;
	if (etrk_[4].phi_layer1_sh_even>-9  and etrk_[4].phi_layer6_sh_even>-9 and etrk_[7].phi_layer1_sh_even>-9  and etrk_[7].phi_layer6_sh_even>-9){
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_sh_even, etrk_[4].phi_layer1_sh_even, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_sh_even, etrk_[7].phi_layer1_sh_even, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_fit_even, etrk_[4].phi_layer1_fit_even, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_fit_even, etrk_[7].phi_layer1_fit_even, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(etrk_[4].phiM_fitxyz_even, etrk_[7].phiM_fitxyz_even);
		}
	}
     }else if ((etrk_[4].has_csc_sh&2)>0 and (etrk_[7].has_csc_sh&1)>0){
        gp1=gp_sh_even[4];
        gp2=gp_sh_odd[7];
        gv1=gv_sh_even[4];
        gv2=gv_sh_odd[7];
	csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
	xfactor = (gp2.perp()/gp1.perp()-1.0)/(839.386-683.786);
	xfactor_st1 = (12.7)*xfactor;
	xfactor_st2 = (12.7)*xfactor/(xfactor*(839.386-683.786)+1);
	npar=3;
	if (etrk_[4].phi_layer1_sh_even>-9  and etrk_[4].phi_layer6_sh_even>-9 and etrk_[7].phi_layer1_sh_odd>-9  and etrk_[7].phi_layer6_sh_odd>-9){
	    phiM_st1_sh = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_sh_even, etrk_[4].phi_layer1_sh_even, xfactor_st1);
	    phiM_st2_sh = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_sh_odd, etrk_[7].phi_layer1_sh_odd, xfactor_st2);
	    phiM_st1_smear3 = PhiMomentum_Xfactor_V2(etrk_[4].phi_layer6_fit_even, etrk_[4].phi_layer1_fit_even, xfactor_st1);
	    phiM_st2_smear3 = PhiMomentum_Xfactor_V2(etrk_[7].phi_layer6_fit_odd, etrk_[7].phi_layer1_fit_odd, xfactor_st2);
	    if (phiM_st1_sh>-9 and phiM_st2_sh>-9){
		csc_bending_angle12_xfactor = deltaPhi(phiM_st1_sh,phiM_st2_sh);
		csc_bending_angle12_xfactor_smear3 = deltaPhi(phiM_st1_smear3,phiM_st2_smear3);
		csc_bending_angle12_xfactor_L1_1 = deltaPhi(etrk_[4].phiM_fitxyz_even, etrk_[7].phiM_fitxyz_odd);
		}
	}
     }
     //std::cout <<"low eta region npar "<< npar <<" eta_st2 "<< gp2.eta() <<" csc bending angle12 "<< csc_bending_angle_12<<" xfactor1 "<< xfactor_st1<<" xfactor2 "<< xfactor_st2 <<" phiM_st1_sh "<< phiM_st1_sh <<" phiM_st2_sh "<< phiM_st1_sh <<" csc_bending_angle12_xfactor "<< csc_bending_angle12_xfactor <<" phiM_st1_smear3 "<< phiM_st1_smear3 <<" phiM_st2_smear3 "<< phiM_st2_smear3<<" csc_bending_angle12_xfactor_smear3 "<< csc_bending_angle12_xfactor_smear3 << std::endl;
  }



  if (fabs(csc_bending_angle12_xfactor)<4) {
     if ((etrk_[8].has_csc_sh || etrk_[9].has_csc_sh) || (etrk_[10].has_csc_sh || etrk_[11].has_csc_sh))
     	etrk_[0].hasSt3orSt4_sh=true; 
     etrk_[0].npar = npar;
     etrk_[0].npar_jose =((npar%2==0 and (etrk_[8].has_csc_sh&2)>0)or(npar%2==1 and (etrk_[8].has_csc_sh&1)>0))? npar:-1;
     etrk_[0].pt_direction_sh=Ptassign_Direction(csc_bending_angle_12, gp2.eta(), npar);  
     etrk_[0].ptphi_diff_sh = csc_bending_angle_12;
     etrk_[0].csc_bending_angle12_gemcsc = csc_bending_angle12_gemcsc;
     etrk_[0].csc_bending_angle12_xfactor = csc_bending_angle12_xfactor;
     etrk_[0].csc_bending_angle12_xfactor_smear0 = csc_bending_angle12_xfactor_smear0;
     etrk_[0].csc_bending_angle12_xfactor_smear1 = csc_bending_angle12_xfactor_smear1;
     etrk_[0].csc_bending_angle12_xfactor_smear2 = csc_bending_angle12_xfactor_smear2;
     etrk_[0].csc_bending_angle12_xfactor_smear3 = csc_bending_angle12_xfactor_smear3;
     etrk_[0].csc_bending_angle12_xfactor_L1_1 = csc_bending_angle12_xfactor_L1_1;
     etrk_[0].csc_bending_angle12_xfactor_L1_2 = csc_bending_angle12_xfactor_L1_2;
     etrk_[0].eta_st1_sh = gp1.eta();
     etrk_[0].eta_st2_sh = gp2.eta();
     
     //double gv1_phi_gauss = CLHEP::RandGauss::shoot(gv1.phi(), .00055);
     //double gv2_phi_gauss = CLHEP::RandGauss::shoot(gv2.phi(), .00095);
     GlobalPoint gp1_smeared( GlobalPoint::Cylindrical(gp1.perp(), gp1_phi_gauss, gp1.z()) );
     GlobalPoint gp2_smeared( GlobalPoint::Cylindrical(gp2.perp(), gp2_phi_gauss, gp2.z()) );

     etrk_[0].pt_direction_gemcsc=Ptassign_Direction(csc_bending_angle12_gemcsc, gp2.eta(), npar);  
     etrk_[0].pt_direction_gemcsc_central=Ptassign_Direction(csc_bending_angle12_gemcsc_central, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor=Ptassign_Direction(csc_bending_angle12_xfactor, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_smear0=Ptassign_Direction(csc_bending_angle12_xfactor_smear0, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_smear1=Ptassign_Direction(csc_bending_angle12_xfactor_smear1, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_smear2=Ptassign_Direction(csc_bending_angle12_xfactor_smear2, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_smear3=Ptassign_Direction(csc_bending_angle12_xfactor_smear3, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_L1_1=Ptassign_Direction(csc_bending_angle12_xfactor_L1_1, gp2.eta(), npar);  
     etrk_[0].pt_direction_xfactor_L1_2=Ptassign_Direction(csc_bending_angle12_xfactor_L1_2, gp2.eta(), npar);  
     if (verbose_)
     	std::cout <<" pt direction sh "<< etrk_[0].pt_direction_sh <<" pt_direction_gemcsc "<<etrk_[0].pt_direction_gemcsc<<"  pt_direction_gemcsc_central "<<  etrk_[0].pt_direction_gemcsc_central <<" pt_direction_xfactor "<< etrk_[0].pt_direction_xfactor <<" xfactor smear0 "<< etrk_[0].pt_direction_xfactor_smear0<<" xfactor smear1 "<< etrk_[0].pt_direction_xfactor_smear1<<" xfactro smear2 "<< etrk_[0].pt_direction_xfactor_smear2<<" L1 1 "<< etrk_[0].pt_direction_xfactor_L1_1 <<" L1 2 "<< etrk_[0].pt_direction_xfactor_L1_2 << std::endl;

     if (etrk_[8].has_csc_sh){
        double gp3_phi_gauss = CLHEP::RandGauss::shoot(gp3.phi(), .00075);
        GlobalPoint gp3_smeared( GlobalPoint::Cylindrical(gp3.perp(), gp3_phi_gauss, gp3.z()) );
	etrk_[0].pt_position_sh=Ptassign_Position_gp(gp1, gp2, gp3, gp2.eta(), npar); //t.momentum().eta() 
	etrk_[0].pt_position_smeared=Ptassign_Position_gp(gp1_smeared, gp2_smeared, gp3_smeared, gp2_smeared.eta(), npar); //t.momentum().eta() 
	if (npar >= 2 and etrk_[0].pt_position_sh>15 and   etrk_[0].pt_position_smeared<10 and verbose_)
	   std::cout <<" npar "<< npar <<" simpt "<< etrk_[0].pt <<" simeta "<< etrk_[0].eta <<" pt_position_sh "<< etrk_[0].pt_position_sh <<" pt_position_smeared "<< etrk_[0].pt_position_smeared <<" gp1 phi "<< gp1.phi() <<" after smearing "<< gp1_smeared.phi() <<" gp2 phi "<< gp2.phi() <<" eta "<< gp2.eta() <<" after smearing "<< gp2_smeared.phi() <<" eta "<< gp2_smeared.eta() <<" gp3 phi "<< gp3.phi() <<" after smearing "<< gp3_smeared.phi() << std::endl;
     }
  } 
  

  if (etrk_[1].has_lct>0 and etrk_[6].has_lct>0 and etrk_[8].has_lct>0){
     int npar=-1;
     GlobalPoint gp1,gp2, gp3;
     if ((etrk_[1].has_lct&1)>0 and (etrk_[6].has_lct&2)>0 and (etrk_[8].has_lct&2)>0){
        gp1=gp_lct_odd[1];
        gp2=gp_lct_even[6];
        gp3=gp_lct_even[8];
	npar=0;
     }else if ((etrk_[1].has_lct&1)>0 and (etrk_[6].has_lct&1)>0 and (etrk_[8].has_lct&1)>0){ 
        gp1=gp_lct_odd[1];
        gp2=gp_lct_odd[6];
        gp3=gp_lct_odd[8];
	npar=1;
     }else if ((etrk_[1].has_lct&2)>0 and (etrk_[6].has_lct&2)>0 and (etrk_[8].has_lct&2)>0){ 
        gp1=gp_lct_even[1];
        gp2=gp_lct_even[6];
        gp3=gp_lct_even[8];
	npar=2;
     }else if ((etrk_[1].has_lct&2)>0 and (etrk_[6].has_lct&1)>0 and (etrk_[8].has_lct&1)>0){ 
        gp1=gp_lct_even[1];
        gp2=gp_lct_odd[6];
        gp3=gp_lct_odd[8];
	npar=3;
     }
     etrk_[0].hasSt1St2St3=true; 
     etrk_[0].pt_position=Ptassign_Position_gp(gp1, gp2, gp3, etrk_[0].eta, npar);  
  
  } 
  
 //general propagation 
  auto propagate_odd_gp(match_track.simTrackPropagateGPs_odd());
  auto propagate_even_gp(match_track.simTrackPropagateGPs_even());
  auto propagate_interstat_odd(match_track.interStatPropagation_odd());
  auto propagate_interstat_even(match_track.interStatPropagation_even());
  for (auto s: stations_to_use_)
  {
    	auto cscdet(cscStationsCo_.at(s));
	if (cscdet.first<1 or cscdet.first>4)  continue;//not in st[1,4]
	int st = cscdet.first;
        auto odd(propagate_odd_gp.at(st-1));
//	std::cout <<" station = "<< cscdet.first <<"  ring = " << cscdet.second << std::endl;
	//take odd chamber as default one
        if (st==1)  {etrk_[s].eta_propagated_ME1 = odd.first; etrk_[s].phi_propagated_ME1 = odd.second;}
        if (st==2)  {etrk_[s].eta_propagated_ME2 = odd.first; etrk_[s].phi_propagated_ME2 = odd.second;}
        if (st==3)  {etrk_[s].eta_propagated_ME3 = odd.first; etrk_[s].phi_propagated_ME3 = odd.second;}
        if (st==4)  {etrk_[s].eta_propagated_ME4 = odd.first; etrk_[s].phi_propagated_ME4 = odd.second;}
        if (st==2 && !isnan(propagate_interstat_odd[12].eta()))  
	              {etrk_[s].eta_interStat12 = propagate_interstat_odd[12].eta(); 
	               etrk_[s].phi_interStat12 = propagate_interstat_odd[12].phi();}
	if (st==3 && !isnan(propagate_interstat_odd[23].eta()))  
	              {etrk_[s].eta_interStat23 = propagate_interstat_odd[23].eta(); 
	               etrk_[s].phi_interStat23 = propagate_interstat_odd[23].phi();}
	if (st==3 && !isnan(propagate_interstat_odd[13].eta()))  
	              {etrk_[s].eta_interStat13 = propagate_interstat_odd[13].eta();
	               etrk_[s].phi_interStat13 = propagate_interstat_odd[13].phi();}
  }
  if (match_track.tfTracks().size()) {
    etrk_[0].has_tfTrack = 1;
    TFTrack* besttrack = match_track.bestTFTrack();
    etrk_[0].trackpt = besttrack->pt();
    etrk_[0].tracketa = besttrack->eta();
    etrk_[0].trackphi = besttrack->phi();
  //  quality_packed;
   etrk_[0].pt_packed = besttrack->ptPacked();
   etrk_[0].eta_packed = besttrack->etaPacked();
   etrk_[0].phi_packed = besttrack->phiPacked();
   etrk_[0].quality_packed = besttrack->qPacked();
 // rank = 0;
   etrk_[0].deltaphi12 = besttrack->dPhi12();
   etrk_[0].deltaphi23 = besttrack->dPhi23();
   etrk_[0].hasME1 = besttrack->hasStubEndcap(1);
   etrk_[0].hasME2 = besttrack->hasStubEndcap(2);
   etrk_[0].nstubs = besttrack->nStubs();
   etrk_[0].deltaR = besttrack->dr();
   etrk_[0].chargesign = besttrack->chargesign();
   unsigned int lct1 = 999;
   auto me1b(besttrack->digiInME(1,1));
   auto me1a(besttrack->digiInME(1,4));
   if (me1a != 999) lct1 = me1a;
   if (me1b != 999) lct1 = me1b;

   if (lct1 < (besttrack->getTriggerDigis()).size()) 
   {
       auto id_me1((besttrack->getTriggerDigisIds()).at(lct1));
       if (id_me1.station() != 1 ) std::cout <<"Error!  CSCDetid should be in station1 " << id_me1 << std::endl;
  //     std::cout <<" CSCDetid in station1 " << id_me1 << std::endl;
       if (id_me1.chamber()%2 == 1)  etrk_[0].chamberME1 |= 1;
       if (id_me1.chamber()%2 == 0)  etrk_[0].chamberME1 |= 2;
       etrk_[0].ME1_ring = id_me1.ring();
       etrk_[0].passGE11 = besttrack->passDPhicutTFTrack(1,bendingcutPt_);
       etrk_[0].passGE11_pt5 = besttrack->passDPhicutTFTrack(1, 5);
       etrk_[0].passGE11_pt7 = besttrack->passDPhicutTFTrack(1, 7);
       etrk_[0].passGE11_pt10 = besttrack->passDPhicutTFTrack(1, 10);
       etrk_[0].passGE11_pt15 = besttrack->passDPhicutTFTrack(1, 15);
       etrk_[0].passGE11_pt20 = besttrack->passDPhicutTFTrack(1, 20);
       etrk_[0].passGE11_pt30 = besttrack->passDPhicutTFTrack(1, 30);
       etrk_[0].passGE11_pt40 = besttrack->passDPhicutTFTrack(1, 40);
       etrk_[0].dphiGE11 = ((besttrack->getTriggerDigis()).at(lct1))->getGEMDPhi();
       etrk_[0].ME1_hs = ((besttrack->getTriggerDigis()).at(lct1))->getStrip();
       etrk_[0].ME1_wg = ((besttrack->getTriggerDigis()).at(lct1))->getKeyWG();
       etrk_[0].passGE11_simpt = match_lct.passDPhicut(id_me1, etrk_[0].chargesign, etrk_[0].dphiGE11, pt);
       //std::cout <<" pass dphicut ?? " <<(etrk_[0].passGE11 ? "  Yes ":" No") << std::endl;
       //if (fabs(etrk_[0].dphiGE11)>1 and fabs(etrk_[0].dphiGE11)<99) std::cout <<" dphiGE11 " << etrk_[0].dphiGE11  << std::endl;
       //if (!etrk_[0].passGE11_simpt and etrk_[0].passGE11 and id_me1.ring()==1) std::cout <<"simpt dphicut failed,st "<< id_me1.station()<<(id_me1.chamber()%2==1 ? " odd": " even") <<" dphiGE11 " << etrk_[0].dphiGE11 << " simpt "<<pt <<" trackpt "<<etrk_[0].trackpt << std::endl; 
       //if (etrk_[0].passGE11_simpt and !etrk_[0].passGE11 and  id_me1.ring()==1) std::cout <<"trackpt dphicut failed,st "<< id_me1.station()<<(id_me1.chamber()%2==1 ? " odd": " even") <<" dphiGE11 " << etrk_[0].dphiGE11 << " simpt "<<pt <<" trackpt "<<etrk_[0].trackpt << std::endl; 

   }

   unsigned int lct2 = besttrack->digiInME(2,1);

   if (lct2 < (besttrack->getTriggerDigis()).size()) 
   {
       auto id_me2((besttrack->getTriggerDigisIds()).at(lct2));
       if (id_me2.station() != 2) std::cout <<"Error!  CSCDetid should be in station2 " << id_me2 << std::endl;
   //    std::cout <<" CSCDetid in station2 ring1 " << id_me2 << std::endl;
       if (id_me2.chamber()%2 == 1)  etrk_[0].chamberME2 |= 1;
       if (id_me2.chamber()%2 == 0)  etrk_[0].chamberME2 |= 2;
       etrk_[0].ME2_ring = id_me2.ring();
       etrk_[0].passGE21 = besttrack->passDPhicutTFTrack(2, bendingcutPt_);
       etrk_[0].passGE21_pt5 = besttrack->passDPhicutTFTrack(2, 5);
       etrk_[0].passGE21_pt7 = besttrack->passDPhicutTFTrack(2, 7);
       etrk_[0].passGE21_pt10 = besttrack->passDPhicutTFTrack(2, 10);
       etrk_[0].passGE21_pt15 = besttrack->passDPhicutTFTrack(2, 15);
       etrk_[0].passGE21_pt20 = besttrack->passDPhicutTFTrack(2, 20);
       etrk_[0].passGE21_pt30 = besttrack->passDPhicutTFTrack(2, 30);
       etrk_[0].passGE21_pt40 = besttrack->passDPhicutTFTrack(2, 40);
       etrk_[0].dphiGE21 = ((besttrack->getTriggerDigis()).at(lct2))->getGEMDPhi();
       etrk_[0].ME2_hs = ((besttrack->getTriggerDigis()).at(lct2))->getStrip();
       etrk_[0].ME2_wg = ((besttrack->getTriggerDigis()).at(lct2))->getKeyWG();
       etrk_[0].passGE21_simpt = match_lct.passDPhicut(id_me2, etrk_[0].chargesign, etrk_[0].dphiGE21, pt);
       //std::cout <<" pass dphicut ?? " <<(etrk_[0].passGE21 ? "  Yes ":" No") << std::endl;
       //if (fabs(etrk_[0].dphiGE21)>1 and fabs(etrk_[0].dphiGE21)<99) std::cout <<" dphiGE21 " << etrk_[0].dphiGE21  << std::endl;
       //if (!etrk_[0].passGE21_simpt and etrk_[0].passGE21 and id_me2.ring()==1) std::cout <<"simpt dphicut failed,st "<<id_me2.station()<<(id_me2.chamber()%2==1 ? " odd": " even")  << " dphiGE21 " << etrk_[0].dphiGE21 << " simpt "<<pt <<" trackpt "<<etrk_[0].trackpt << std::endl; 
       //if (etrk_[0].passGE21_simpt and !etrk_[0].passGE21 and id_me2.ring()==1) std::cout <<"trackpt dphicut failed,st "<<id_me2.station() <<(id_me2.chamber()%2==1 ? " odd": " even") <<" dphiGE21 " << etrk_[0].dphiGE21 << " simpt "<<pt <<" trackpt "<<etrk_[0].trackpt << std::endl; 

   }
    auto triggerDigiIds(besttrack->getTriggerDigisIds()); 
    auto triggerDigis(besttrack->getTriggerDigis()); 

    auto triggerDigiEtaPhi(besttrack->getTriggerEtaPhis());
    if (triggerDigiIds.size() == triggerDigiEtaPhi.size() && triggerDigis.size() == triggerDigiIds.size())
     {
        bool stub_Good_ME[4] = {1,1,1,1};

	for (unsigned int i=0; i<triggerDigiIds.size(); i++)
	{
	  auto id(triggerDigiIds.at(i));
	  auto etaphi(triggerDigiEtaPhi.at(i));

	   //std::cout << " CSCDetId " << id << " to fill Tree, Phi " << etaphi.second << std::endl;
	  int st = id.station();
	  bool IsOdd(id.chamber()%2==1);
	  if (IsOdd)
	  {
          auto odd(propagate_odd_gp.at(st-1));
	 // std::cout <<"  propagated position in odd chamber eta:"  << odd.first << "  phi:" << odd.second << std::endl;
          if (st==1)  {etrk_[0].eta_propagated_ME1 = odd.first; etrk_[0].phi_propagated_ME1 = odd.second;}
          if (st==2)  {etrk_[0].eta_propagated_ME2 = odd.first; etrk_[0].phi_propagated_ME2 = odd.second;}
          if (st==3)  {etrk_[0].eta_propagated_ME3 = odd.first; etrk_[0].phi_propagated_ME3 = odd.second;}
          if (st==4)  {etrk_[0].eta_propagated_ME4 = odd.first; etrk_[0].phi_propagated_ME4 = odd.second;}
	  
	  if (st==2 && !isnan(propagate_interstat_odd[12].eta()))  
	              {etrk_[0].eta_interStat12 = propagate_interstat_odd[12].eta(); 
	               etrk_[0].phi_interStat12 = propagate_interstat_odd[12].phi();}
	  if (st==3 && !isnan(propagate_interstat_odd[23].eta()))  
	              {etrk_[0].eta_interStat23 = propagate_interstat_odd[23].eta(); 
	               etrk_[0].phi_interStat23 = propagate_interstat_odd[23].phi();}
	  if (st==3 && !isnan(propagate_interstat_odd[13].eta()))  
	              {etrk_[0].eta_interStat13 = propagate_interstat_odd[13].eta();
	               etrk_[0].phi_interStat13 = propagate_interstat_odd[13].phi();}
           }
	  else {
          auto even(propagate_even_gp.at(st-1));
	  //std::cout <<"  propagated position in even chamber eta:"  << even.first << "  phi:" << even.second << std::endl;
          if (st==1)  {etrk_[0].eta_propagated_ME1 = even.first; etrk_[0].phi_propagated_ME1 = even.second;}
          if (st==2)  {etrk_[0].eta_propagated_ME2 = even.first; etrk_[0].phi_propagated_ME2 = even.second;}
          if (st==3)  {etrk_[0].eta_propagated_ME3 = even.first; etrk_[0].phi_propagated_ME3 = even.second;}
          if (st==4)  {etrk_[0].eta_propagated_ME4 = even.first; etrk_[0].phi_propagated_ME4 = even.second;}
	  
	  if (st==2 && !isnan(propagate_interstat_even[12].eta()))  
	              {etrk_[0].eta_interStat12 = propagate_interstat_even[12].eta(); 
	               etrk_[0].phi_interStat12 = propagate_interstat_even[12].phi();}
	  if (st==3 && !isnan(propagate_interstat_even[23].eta()))  
	              {etrk_[0].eta_interStat23 = propagate_interstat_even[23].eta(); 
	               etrk_[0].phi_interStat23 = propagate_interstat_even[23].phi();}
	  if (st==3 && !isnan(propagate_interstat_even[13].eta()))  
	              {etrk_[0].eta_interStat13 = propagate_interstat_even[13].eta();
	               etrk_[0].phi_interStat13 = propagate_interstat_even[13].phi();}
		       

	  }
          if (st==1)  {etrk_[0].eta_ME1_TF = etaphi.first; etrk_[0].phi_ME1_TF = etaphi.second;
	               stub_Good_ME[0] = match_lct.checkStubInChamber(id,*triggerDigis.at(i));}
          if (st==2)  {etrk_[0].eta_ME2_TF = etaphi.first; etrk_[0].phi_ME2_TF = etaphi.second;
	               stub_Good_ME[1] = match_lct.checkStubInChamber(id,*triggerDigis.at(i));}
          if (st==3)  {etrk_[0].eta_ME3_TF = etaphi.first; etrk_[0].phi_ME3_TF = etaphi.second;
	               stub_Good_ME[2] = match_lct.checkStubInChamber(id,*triggerDigis.at(i));}
          if (st==4)  {etrk_[0].eta_ME4_TF = etaphi.first; etrk_[0].phi_ME4_TF = etaphi.second;
	               stub_Good_ME[3] = match_lct.checkStubInChamber(id,*triggerDigis.at(i));}

	  
	  //if ( match_lct.checkStubInChamber(id,*triggerDigis.at(i))) std::cout << "stub in TF can be matched to simtrack" << std::endl;
	  //else std::cout << "stub in TF can NOT be matched to simtrack" << std::endl;
	}
         etrk_[0].allstubs_matched_TF = (stub_Good_ME[0] and stub_Good_ME[1] and stub_Good_ME[2] and stub_Good_ME[3]);
	
	 /*if (!stub_Good_ME[0]) std::cout << "In station1 stub can not be matched to simTrack" << std::endl; 
	 if (!stub_Good_ME[1]) std::cout << "In station2 stub can not be matched to simTrack" << std::endl; 
	 if (!stub_Good_ME[2]) std::cout << "In station3 stub can not be matched to simTrack" << std::endl; 
	 if (!stub_Good_ME[3]) std::cout << "In station4 stub can not be matched to simTrack" << std::endl; 
	 //for debug
	if (!etrk_[0].allstubs_matched_TF && abs(etrk_[0].eta)>1.65 && abs(etrk_[0].eta)<1.85)
	 for (unsigned int i=0; i<triggerDigiIds.size(); i++)
	 {
	  auto id(triggerDigiIds.at(i));
          std::cout << "stub in TF DetId " << triggerDigiIds.at(i) << "  " << *triggerDigis.at(i) << std::endl;
	  std::cout << "matched stub in this Detid " << std::endl;
	  for (auto p : match_lct.lctsInChamber(id.rawId()))  std::cout << "  " << p << std::endl;
	  if ( match_lct.checkStubInChamber(id,*triggerDigis.at(i))) std::cout << "stub in TF can be matched to simtrack" << std::endl;
	  else std::cout << "stub in TF can NOT be matched to simtrack" << std::endl;

	   }*/
         // check simhit in each station, station1->bit1, station2->bit2
	 if (etrk_[1].has_csc_sh>0 or etrk_[4].has_csc_sh>0 or etrk_[5].has_csc_sh>0) etrk_[0].has_csc_sh |= 1;
	 //std::cout << "simhits in station1 " << (std::bitset<8>)etrk_[0].has_csc_sh  << std::endl;
	 if (etrk_[6].has_csc_sh>0 or etrk_[7].has_csc_sh>0) etrk_[0].has_csc_sh |= 2;
	 //std::cout << "simhits in station12 " << (std::bitset<8>)etrk_[0].has_csc_sh  << std::endl;
	 if (etrk_[8].has_csc_sh>0 or etrk_[9].has_csc_sh>0) etrk_[0].has_csc_sh |= 4;
	 //std::cout << "simhits in station123 " << (std::bitset<8>)etrk_[0].has_csc_sh  << std::endl;
	 if (etrk_[10].has_csc_sh>0 or etrk_[11].has_csc_sh>0) etrk_[0].has_csc_sh |= 8;
	 //std::cout << "simhits in each station1234 " << (std::bitset<8>)etrk_[0].has_csc_sh  << std::endl;
     }//end if 

    if (triggerDigiEtaPhi.size()>1)
    {
         auto etaphi1(triggerDigiEtaPhi.at(0));
	 auto etaphi2(triggerDigiEtaPhi.at(1));
	 etrk_[0].lctdphi12 = etaphi1.second-etaphi2.second;
    
    }
   /*std::cout<<"check csc detids" << std::endl;
   for (auto id : besttrack->getTriggerDigisIds()) 
   {
      std::cout<<" stub id " << id << std::endl;
   }   
   std::cout<< std::endl;*/
     //    std::cout << "SimTrack has matched CSCTF track" << std::endl;
  }
  
  if (match_track.tfCands().size()) {
    etrk_[0].has_tfCand = 1;
    std::cout << "SimTrack has matched CSCTF Cand" << std::endl;
  }
  
  if (match_track.gmtRegCands().size()) {
    etrk_[0].has_gmtRegCand = 1;
    std::cout << "SimTrack has GMTRegCand" << std::endl;
  }

  if (match_track.gmtCands().size()) {
    etrk_[0].has_gmtCand = 1;
    std::cout << "SimTrack has GMTCand" << std::endl;
  }

  // L1Extra
  auto l1Extras(match_l1_gmt.getMatchedL1ExtraMuonParticles());
  if (l1Extras.size()) {
    etrk_[0].has_l1Extra = 1;

    auto l1Extra(l1Extras[0].first);
    etrk_[0].l1Extra_pt = l1Extra.pt();
    etrk_[0].l1Extra_eta = l1Extra.eta();
    etrk_[0].l1Extra_phi = l1Extra.phi();
    etrk_[0].l1Extra_dR = l1Extras[0].second;
    if (verbose_) {
      std::cout << "Number of matched L1Extras: " << l1Extras.size() << std::endl;
      std::cout << "l1Extra_pt " << etrk_[0].l1Extra_pt << std::endl;
      std::cout << "l1Extra_eta " << etrk_[0].l1Extra_eta << std::endl;
      std::cout << "l1Extra_phi " << etrk_[0].l1Extra_phi << std::endl;
      std::cout << "l1Extra_dR " << etrk_[0].l1Extra_dR << std::endl;
    }
  }

  // RecoTrackExtra
  auto recoTrackExtras(match_hlt_track.getMatchedRecoTrackExtras());
  if (recoTrackExtras.size()) {
    if (verbose_) std::cout << "Number of matched RecoTrackExtras: " << recoTrackExtras.size() << std::endl;
    etrk_[0].has_recoTrackExtra = 1;

    auto recoTrackExtra(recoTrackExtras[0]);
    etrk_[0].recoTrackExtra_pt_inner = recoTrackExtra.innerMomentum().Rho();
    etrk_[0].recoTrackExtra_eta_inner = recoTrackExtra.innerPosition().eta();
    etrk_[0].recoTrackExtra_phi_inner = recoTrackExtra.innerPosition().phi();

    etrk_[0].recoTrackExtra_pt_outer = recoTrackExtra.outerMomentum().Rho();
    etrk_[0].recoTrackExtra_eta_outer = recoTrackExtra.outerPosition().eta();
    etrk_[0].recoTrackExtra_phi_outer = recoTrackExtra.outerPosition().phi();
  }

  // RecoTrack
  auto recoTracks(match_hlt_track.getMatchedRecoTracks());
  if (match_hlt_track.getMatchedRecoTracks().size()) {
    if (verbose_) std::cout << "Number of matched RecoTracks: " << recoTracks.size() << std::endl;
    etrk_[0].has_recoTrack = 1;

    auto recoTrack(recoTracks[0]);
    etrk_[0].recoTrack_pt_outer = recoTrack.outerPt();
    etrk_[0].recoTrack_eta_outer = recoTrack.outerEta();
    etrk_[0].recoTrack_phi_outer = recoTrack.outerPhi();
  }

  // RecoChargedCandidate
  auto recoChargedCandidates(match_hlt_track.getMatchedRecoChargedCandidates());
  if (recoChargedCandidates.size()) {
    if (verbose_) std::cout << "Number of matched RecoChargedCandidates: " << recoChargedCandidates.size() << std::endl;
    etrk_[0].has_recoChargedCandidate = 1;

    auto recoChargedCandidate(recoChargedCandidates[0]);
    etrk_[0].recoChargedCandidate_pt = recoChargedCandidate.pt();
    etrk_[0].recoChargedCandidate_eta = recoChargedCandidate.eta();
    etrk_[0].recoChargedCandidate_phi = recoChargedCandidate.phi();
    etrk_[0].recoChargedCandidate_nValidDTHits = (recoChargedCandidate.track().get())->hitPattern().numberOfValidMuonDTHits();
    etrk_[0].recoChargedCandidate_nValidCSCHits = (recoChargedCandidate.track().get())->hitPattern().numberOfValidMuonCSCHits();
    etrk_[0].recoChargedCandidate_nValidRPCHits = (recoChargedCandidate.track().get())->hitPattern().numberOfValidMuonRPCHits();
    if (verbose_) {
      std::cout << "recoChargedCandidate_pt " << etrk_[0].recoChargedCandidate_pt << std::endl;
      std::cout << "recoChargedCandidate_eta " << etrk_[0].recoChargedCandidate_eta << std::endl;
      std::cout << "recoChargedCandidate_phi " << etrk_[0].recoChargedCandidate_phi << std::endl;
      std::cout << "nValidHits:" 
		<< " DT " << etrk_[0].recoChargedCandidate_nValidDTHits 
		<< " CSC " << etrk_[0].recoChargedCandidate_nValidCSCHits
		<< " RPC " << etrk_[0].recoChargedCandidate_nValidRPCHits << std::endl;
    }
  }

  for (auto s: stations_to_use_)
  {
    tree_eff_[s]->Fill();
  }
}



void GEMCSCAnalyzer::analyzeTrackChamberDeltas(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  const SimTrack &t = match_sh.trk();

  if (verbose_ > 1) // ---- SimHitMatcher debug printouts
  {
    cout<<"** GEM SimHits **"<<endl;
    cout<<"n_sh_ids "<<match_sh.detIdsGEM().size()<<endl;
    cout<<"n_sh_ids_copad "<<match_sh.detIdsGEMCoincidences().size()<<endl;
    auto gem_sh_sch_ids = match_sh.superChamberIdsGEM();
    cout<<"n_sh_ids_sch "<<gem_sh_sch_ids.size()<<endl;
    cout<<"n_sh_ids_cosch "<<match_sh.superChamberIdsGEMCoincidences().size()<<endl;
    cout<<"n_sh_pad "<<match_sh.nPadsWithHits()<<endl;
    cout<<"n_sh_copad "<<match_sh.nCoincidencePadsWithHits()<<endl;
    for (auto id: gem_sh_sch_ids)
    {
      auto gem_simhits = match_sh.hitsInSuperChamber(id);
      auto gem_simhits_gp = match_sh.simHitsMeanPosition(gem_simhits);
      cout<<"shtrk "<<trk_no<<": "<<t.momentum().eta()<<" "<<t.momentum().phi()<<" "<<t.vertIndex()
          <<" | "<<gem_simhits.size()<<" "<<gem_simhits_gp.phi()<<endl;
    }

    const int nsch(match_sh.superChamberIdsGEM().size());
    auto gem_sh_ids = match_sh.detIdsGEM();
    for(auto d: gem_sh_ids)
    {
      GEMDetId id(d);
      auto strips = match_sh.hitStripsInDetId(d);
      for(auto s: strips)
      {
        cout<<"sch_strip "<<nsch<<" "<<s<<" "<<id.roll()<<" "<<id.chamber()<<" "<<strips.size()<<endl;
        //if (nsch > 1)cout<<"many_sch_strip "<<s<<" "<<id.roll()<<" "<<id.chamber()<<endl;
        //if (nsch == 1)cout<<"1_sch_strip "<<s<<" "<<id.roll()<<endl;
      }
    }

    cout<<"** CSC SimHits **"<<endl;
    cout<<"n_csh_ids "<<match_sh.detIdsCSC().size()<<endl;
    auto csc_csh_ch_ids = match_sh.chamberIdsCSC();
    cout<<"n_csh_ids_ch "<<csc_csh_ch_ids.size()<<endl;
    cout<<"n_csh_coch "<<match_sh.nCoincidenceCSCChambers(minNHitsChamberCSCSimHit_)<<endl;
    for (auto id: csc_csh_ch_ids)
    {
      auto csc_simhits = match_sh.hitsInChamber(id);
      auto csc_simhits_gp = match_sh.simHitsMeanPosition(csc_simhits);
      cout<<"cshtrk "<<trk_no<<": "<<t.momentum().eta()<<" "<<t.momentum().phi()
          <<" | "<<csc_simhits.size()<<" "<<csc_simhits_gp.phi()<<endl;
    }

    const int ncch(match_sh.chamberIdsCSC().size());
    auto csc_sh_ids = match_sh.detIdsCSC();
    for(auto d: csc_sh_ids)
    {
      CSCDetId id(d);
      auto strips = match_sh.hitStripsInDetId(d);
      for(auto s: strips)
      {
        cout<<"cscch_strip "<<ncch<<" "<<s<<" "<<id.chamber()<<" "<<strips.size()<<endl;
      }
    }
  }

  if (verbose_ > 1) // ---- GEMDigiMatcher debug printouts
  {
    cout<<"** GEM Digis **"<<endl;
    cout<<"n_gd_ids "<<match_gd.detIdsDigi().size()<<endl;
    auto gem_gd_sch_ids = match_gd.superChamberIdsDigi();
    cout<<"n_gd_ids_sch "<<gem_gd_sch_ids.size()<<endl;
    cout<<"n_gd_ids_cosch "<<match_gd.superChamberIdsCoPad().size()<<endl;
    cout<<"n_gd_pad "<<match_gd.nPads()<<endl;
    cout<<"n_gd_copad "<<match_gd.nCoPads()<<endl;
    for (auto id: gem_gd_sch_ids)
    {
      auto gem_digis = match_gd.digisInSuperChamber(id);
      auto gem_digis_gp = match_gd.digisMeanPosition(gem_digis);
      cout<<"gdtrk "<<trk_no<<": "<<t.momentum().eta()<<" "<<t.momentum().phi()<<" "<<t.vertIndex()
          <<" | "<<gem_digis.size()<<" "<<gem_digis_gp.phi()<<endl;
    }
  }

  if (verbose_ > 1) // ---- CSCDigiMatcher debug printouts
  {
    cout<<"** CSC Digis **"<<endl;
    cout<<"n_sd_ids "<<match_cd.detIdsStrip().size()<<endl;
    auto csc_sd_ch_ids = match_cd.chamberIdsStrip();
    cout<<"n_sd_ids_ch "<<csc_sd_ch_ids.size()<<endl;
    //cout<<"n_sd_lay "<<cdm.nLayersWithStripInChamber(id)<<endl;
    cout<<"n_sd_coch "<<match_cd.nCoincidenceStripChambers()<<endl;
    for (auto id: csc_sd_ch_ids)
    {
      auto csc_digis = match_cd.stripDigisInChamber(id);
      auto csc_digis_gp = match_cd.digisMeanPosition(csc_digis);
      cout<<"sdtrk "<<trk_no<<": "<<t.momentum().eta()<<" "<<t.momentum().phi()
          <<" | "<<csc_digis.size()<<" "<<csc_digis_gp.phi()<<endl;
    }

    cout<<"n_wd_ids "<<match_cd.detIdsWire().size()<<endl;
    auto csc_wd_ch_ids = match_cd.chamberIdsWire();
    cout<<"n_wd_ids_ch "<<csc_wd_ch_ids.size()<<endl;
    //cout<<"n_wd_lay "<<cdm.nLayersWithStripInChamber(id)<<endl;
    cout<<"n_wd_coch "<<match_cd.nCoincidenceWireChambers()<<endl;
  }

  // debug possible mismatch in number of pads from digis and simhits
  if (verbose_ > 0 && match_gd.nPads() != match_sh.nPadsWithHits())
  {
    cout<<"mismatch "<<match_sh.nPadsWithHits()<<" "<<match_gd.nPads()<<endl;
    auto gdids = match_gd.detIdsDigi();
    for (auto d: gdids)
    {
      auto pad_ns = match_gd.padNumbersInDetId(d);
      cout<<"gd "<<GEMDetId(d)<<" ";
      copy(pad_ns.begin(), pad_ns.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }
    auto shids = match_sh.detIdsGEM();
    for (auto d: shids)
    {
      auto pad_ns = match_sh.hitPadsInDetId(d);
      cout<<"sh "<<GEMDetId(d)<<" ";
      copy(pad_ns.begin(), pad_ns.end(), ostream_iterator<int>(cout, " "));
      cout<<endl;
    }
  }

  // fill the information for delta-tree
  // only for tracks with enough hit layers in CSC and at least a pad in GEM
  if ( match_gd.nPads() > 0 &&
       match_cd.nCoincidenceStripChambers(minNHitsChamberCSCStripDigi_) > 0 &&
       match_cd.nCoincidenceWireChambers(minNHitsChamberCSCWireDigi_) > 0 )
  {
    dtrk_.pt = t.momentum().pt();
    dtrk_.phi = t.momentum().phi();
    dtrk_.eta = t.momentum().eta();
    dtrk_.charge = t.charge();

    auto csc_sd_ch_ids = match_cd.chamberIdsStrip();
    auto gem_d_sch_ids = match_gd.superChamberIdsDigi();
    if (verbose_) cout<<"will match csc & gem  "<<csc_sd_ch_ids.size()<<" "<<gem_d_sch_ids.size()<<endl;
    for (auto csc_d: csc_sd_ch_ids)
    {
      CSCDetId csc_id(csc_d);

      // require CSC chamber to have at least 4 layers with comparator digis
      if (match_cd.nLayersWithStripInChamber(csc_d) < minNHitsChamberCSCStripDigi_) continue;

      bool is_odd = csc_id.chamber() & 1;
      int region = (csc_id.endcap() == 1) ? 1 : -1;

      auto csc_sh = match_sh.hitsInChamber(csc_d);
      GlobalPoint csc_sh_gp = match_sh.simHitsMeanPosition(csc_sh);

      // CSC trigger strips and wire digis
      auto csc_sd = match_cd.stripDigisInChamber(csc_d);
      auto csc_wd = match_cd.wireDigisInChamber(csc_d);

      GlobalPoint csc_dg_gp = match_cd.digisCSCMedianPosition(csc_sd, csc_wd);

      //GlobalPoint csc_sd_gp = match_cd.digisMeanPosition(csc_sd);
      //cout<<"test csc_dg_gp  "<<csc_sd_gp<<" "<<csc_dg_gp<<" "<<csc_sd_gp.phi() - csc_dg_gp.phi()<<endl;

      if ( std::abs(csc_dg_gp.z()) < 0.001 ) { cout<<"bad csc_dg_gp"<<endl; continue; }

      auto lct_digi = match_lct.lctInChamber(csc_d);
      GlobalPoint csc_lct_gp;
      if (is_valid(lct_digi))
      {
        csc_lct_gp = match_lct.digiPosition(lct_digi);
      }


      // match with signal in GEM in corresponding superchamber
      for(auto gem_d: gem_d_sch_ids)
      {
        GEMDetId gem_id(gem_d);

        // gotta be the same endcap
        if (gem_id.region() != region) continue;
        // gotta be the same chamber#
        if (gem_id.chamber() != csc_id.chamber()) continue;

        auto gem_sh = match_sh.hitsInSuperChamber(gem_d);
        GlobalPoint gem_sh_gp = match_sh.simHitsMeanPosition(gem_sh);

        auto gem_dg = match_gd.digisInSuperChamber(gem_d);
        //GlobalPoint gem_dg_gp = match_gd.digisMeanPosition(gem_dg);
        auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(gem_dg, csc_dg_gp);
        //auto best_gem_dg = gem_dg_and_gp.first;
        GlobalPoint gem_dg_gp = gem_dg_and_gp.second;

        auto gem_pads = match_gd.padsInSuperChamber(gem_d);
        //GlobalPoint gem_pads_gp = match_gd.digisMeanPosition(gem_pads);
        auto gem_pad_and_gp = match_gd.digiInGEMClosestToCSC(gem_pads, csc_dg_gp);
        auto best_gem_pad = gem_pad_and_gp.first;
        GlobalPoint gem_pad_gp = gem_pad_and_gp.second;

        if (gem_sh.size() == 0 || gem_dg.size() == 0 || gem_pads.size() == 0) continue;

        /*
        float avg_roll = 0.;
        for (auto& d: gem_pads )
        {
          GEMDetId id(digi_id(d));
          avg_roll += id.roll();
        }
        avg_roll = avg_roll/gem_pads.size();
        */
        GEMDetId id_of_best_gem(digi_id(best_gem_pad));

        dtrk_.odd = is_odd;
        dtrk_.chamber = csc_id.chamber();
        dtrk_.endcap = csc_id.endcap();
        dtrk_.roll = id_of_best_gem.roll();
        dtrk_.csc_sh_phi = csc_sh_gp.phi();
        dtrk_.csc_dg_phi = csc_dg_gp.phi();
        dtrk_.gem_sh_phi = gem_sh_gp.phi();
        dtrk_.gem_dg_phi = gem_dg_gp.phi();
        dtrk_.gem_pad_phi = gem_pad_gp.phi();
        dtrk_.dphi_sh = deltaPhi(csc_sh_gp.phi(), gem_sh_gp.phi());
        dtrk_.dphi_dg = deltaPhi(csc_dg_gp.phi(), gem_dg_gp.phi());
        dtrk_.dphi_pad = deltaPhi(csc_dg_gp.phi(), gem_pad_gp.phi());
        dtrk_.csc_sh_eta = csc_sh_gp.eta();
        dtrk_.csc_dg_eta = csc_dg_gp.eta();
        dtrk_.gem_sh_eta = gem_sh_gp.eta();
        dtrk_.gem_dg_eta = gem_dg_gp.eta();
        dtrk_.gem_pad_eta = gem_pad_gp.eta();
        dtrk_.deta_sh = csc_sh_gp.eta() - gem_sh_gp.eta();
        dtrk_.deta_dg = csc_dg_gp.eta() - gem_dg_gp.eta();
        dtrk_.deta_pad = csc_dg_gp.eta() - gem_pad_gp.eta();
        dtrk_.bend = -99;
        dtrk_.csc_lct_phi = -99.;
        dtrk_.dphi_lct_pad = -99.;
        dtrk_.csc_lct_eta = -99.;
        dtrk_.deta_lct_pad = -99.;
        if (std::abs(csc_lct_gp.z()) > 0.001)
        {
          dtrk_.bend = LCT_BEND_PATTERN[digi_pattern(lct_digi)];
          dtrk_.csc_lct_phi = csc_lct_gp.phi();
          dtrk_.dphi_lct_pad = deltaPhi(csc_lct_gp.phi(), gem_pad_gp.phi());
          dtrk_.csc_lct_eta = csc_lct_gp.eta();
          dtrk_.deta_lct_pad = csc_lct_gp.eta() - gem_pad_gp.eta();
        }

        tree_delta_->Fill();

        /*
        if (csc_id.endcap()==1)
        {
          auto best_gem_dg = gem_dg_and_gp.first;
          GEMDetId id_of_best_dg(digi_id(best_gem_dg));
          cout<<"funny_deta "<<gem_dg_gp.eta() - gem_pad_gp.eta()<<" "
              <<digi_channel(best_gem_pad)<<" "<<digi_channel(best_gem_dg)<<" "
              <<id_of_best_gem.roll()<<" "<<id_of_best_dg.roll()<<" "
              <<id_of_best_gem.layer()<<" "<<id_of_best_dg.layer()<<" "
              <<match_gd.nLayersWithDigisInSuperChamber(gem_d)<<endl;
        }*/

        if (verbose_ > 1) // debug printout for the stuff in delta-tree
        {
          cout<<"got match "<<csc_id<<"  "<<gem_id<<endl;
          cout<<"matchdphis "<<is_odd<<" "<<csc_id.chamber()<<" "
              <<csc_sh_gp.phi()<<" "<<csc_dg_gp.phi()<<" "<<gem_sh_gp.phi()<<" "<<gem_dg_gp.phi()<<" "<<gem_pad_gp.phi()<<" "
              <<dtrk_.dphi_sh<<" "<<dtrk_.dphi_dg<<" "<<dtrk_.dphi_pad<<"   "
              <<csc_sh_gp.eta()<<" "<<csc_dg_gp.eta()<<" "<<gem_sh_gp.eta()<<" "<<gem_dg_gp.eta()<<" "<<gem_pad_gp.eta()<<" "
              <<dtrk_.deta_sh<<" "<<dtrk_.deta_dg<<" "<<dtrk_.deta_pad<<endl;
        }
      }
    }
  }
}


void GEMCSCAnalyzer::bookSimTracksDeltaTree()
{
  edm::Service< TFileService > fs;
  tree_delta_ = fs->make<TTree>("trk_delta", "trk_delta");
  tree_delta_->Branch("odd", &dtrk_.odd);
  tree_delta_->Branch("charge", &dtrk_.charge);
  tree_delta_->Branch("chamber", &dtrk_.chamber);
  tree_delta_->Branch("endcap", &dtrk_.endcap);
  tree_delta_->Branch("roll", &dtrk_.roll);
  tree_delta_->Branch("bend", &dtrk_.bend);
  tree_delta_->Branch("pt", &dtrk_.pt);
  tree_delta_->Branch("eta", &dtrk_.eta);
  tree_delta_->Branch("phi", &dtrk_.phi);
  tree_delta_->Branch("csc_sh_phi", &dtrk_.csc_sh_phi);
  tree_delta_->Branch("csc_dg_phi", &dtrk_.csc_dg_phi);
  tree_delta_->Branch("gem_sh_phi", &dtrk_.gem_sh_phi);
  tree_delta_->Branch("gem_dg_phi", &dtrk_.gem_dg_phi);
  tree_delta_->Branch("gem_pad_phi", &dtrk_.gem_pad_phi);
  tree_delta_->Branch("dphi_sh", &dtrk_.dphi_sh);
  tree_delta_->Branch("dphi_dg", &dtrk_.dphi_dg);
  tree_delta_->Branch("dphi_pad", &dtrk_.dphi_pad);
  tree_delta_->Branch("csc_sh_eta", &dtrk_.csc_sh_eta);
  tree_delta_->Branch("csc_dg_eta", &dtrk_.csc_dg_eta);
  tree_delta_->Branch("gem_sh_eta", &dtrk_.gem_sh_eta);
  tree_delta_->Branch("gem_dg_eta", &dtrk_.gem_dg_eta);
  tree_delta_->Branch("gem_pad_eta", &dtrk_.gem_pad_eta);
  tree_delta_->Branch("deta_sh", &dtrk_.deta_sh);
  tree_delta_->Branch("deta_dg", &dtrk_.deta_dg);
  tree_delta_->Branch("deta_pad", &dtrk_.deta_pad);
  tree_delta_->Branch("csc_lct_phi", &dtrk_.csc_lct_phi);
  tree_delta_->Branch("dphi_lct_pad", &dtrk_.dphi_lct_pad);
  tree_delta_->Branch("csc_lct_eta", &dtrk_.csc_lct_eta);
  tree_delta_->Branch("deta_lct_pad", &dtrk_.deta_lct_pad);
  //tree_delta_->Branch("", &dtrk_.);
  tree_delta_->Branch("dphi_gem_sh_csc_sh", &dtrk_.dphi_gem_sh_csc_sh);
  tree_delta_->Branch("dphi_gem_dg_csc_dg", &dtrk_.dphi_gem_dg_csc_dg);
  tree_delta_->Branch("dphi_gem_pad_csc_lct", &dtrk_.dphi_gem_pad_csc_lct);
  tree_delta_->Branch("dphi_gem_rh_csc_seg", &dtrk_.dphi_gem_rh_csc_seg);
}


 void GEMCSCAnalyzer::printout(SimTrackMatchManager& match, int trk_no, const char msg[300])
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const RPCDigiMatcher& match_rd = match.rpcDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  const TrackMatcher& match_track = match.tracks();
  const SimTrack &t = match_sh.trk();

  
  std::cout << "======================== matching information ========================= " << std::endl;
  std::cout << msg << std::endl;
  std::cout << "  pt:"<<t.momentum().pt()
            << "  phi:"<<t.momentum().phi()
            << "  eta:"<<t.momentum().eta()
            << "  chage:"<<t.charge() << std::endl;
  
  std::cout << "######matching simhit to simtrack " << std::endl;
  for (auto d: match_sh.chamberIdsCSC(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers = match_sh.nLayersWithHitsInSuperChamber(d);
    const auto& hits = match_sh.hitsInChamber(d);
    auto gp = match_sh.simHitsMeanPosition(hits);
    float mean_strip = match_sh.simHitsMeanStrip(hits);
    std::cout << "CSC Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers<<" global eta:"<<gp.eta()<<" mean strip:"<<mean_strip<<endl;
  }     
  
  if (false) for (auto d: match_sh.chamberIdsRPC())
  {
    RPCDetId id(d);
    const int st(detIdToMEStation(id.station(), id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers = match_sh.nLayersWithHitsInSuperChamber(d);
    const auto& hits = match_sh.hitsInChamber(d);
    auto gp = match_sh.simHitsMeanPosition(hits);
    float mean_strip = match_sh.simHitsMeanStrip(hits);
    std::cout << "RPC Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers<<" global eta:"<<gp.eta()<<" mean strip:"<<mean_strip<<endl;
    int cscchamber = CSCTriggerNumbering::chamberFromTriggerLabels(id.sector(), 0, id.station(), id.subsector());
    std::cout <<"rpc detid " << id << " csc chamebr:"<< cscchamber << std::endl;
  }     

  for(auto d: match_sh.superChamberIdsGEM())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    int nlayers = match_sh.nLayersWithHitsInSuperChamber(d);
    auto gp = match_sh.simHitsMeanPosition(match_sh.hitsInSuperChamber(d));
    float mean_strip = match_sh.simHitsMeanStrip(match_sh.hitsInSuperChamber(d));
    std::cout << "GEM Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers<<" global eta:"<<gp.eta()<<" mean strip:"<<mean_strip<<endl;

  }
  
  std::cout << "######matching Cathode Digi to simtrack " << std::endl;
  for (auto d: match_cd.chamberIdsStrip(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers = match_cd.nLayersWithStripInChamber(d);
    std::cout <<"CSC Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers<<std::endl;
    auto strips = match_cd.stripDigisInChamber(d);
    // std::cout <<"strips:"  ;
    for ( auto p : strips )
      std::cout << p << std::endl;
  }

  std::cout << "######matching Anode Digi to simtrack " << std::endl;
  for (auto d: match_cd.chamberIdsWire(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers = match_cd.nLayersWithWireInChamber(d);
    std::cout <<"CSC Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers<<std::endl;
    auto wires = match_cd.wireDigisInChamber(d);
    //  std::cout <<"WireGroups:"  ;
    for ( auto p : wires)
      std::cout << p <<std::endl; 
  }
  
  std::cout << "######matching GEM Digi to simtrack " << std::endl;
  for(auto d: match_gd.superChamberIdsDigi())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
    
    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    
    int nlayers = match_gd.nLayersWithDigisInSuperChamber(d);
    auto digis = match_gd.digisInSuperChamber(d);
    int median_strip = match_gd.median(digis);
    int hs = match_gd.extrapolateHsfromGEMStrip( d, median_strip);
    std::cout <<"GEM Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers
              <<" Medianstrip in Digi:" <<median_strip<<" hs:" << hs<<std::endl;
    // std::cout <<"GEM Pads:"  ;
    auto pads = match_gd.padsInSuperChamber(d);
    for ( auto p=pads.begin(); p != pads.end(); p++)
      std::cout << "  "<< *p <<std::endl; 
  }
  
  std::cout << "######matching Copad to simtrack " << std::endl;
  for (auto d: match_gd.superChamberIdsCoPad())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
    
    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    
    std::cout <<"Copad GEM Chamber: "<<d<<" "<<id<<std::endl;
    auto Copads = match_gd.coPadsInSuperChamber(d);
    // std::cout <<"GEM Copads:"  ;
    for ( auto p=Copads.begin(); p != Copads.end(); p++)
      {  std::cout << "  "<< *p ; }
    std::cout << std::endl;
  }

  
  std::cout << "######matching RPC Digi to simtrack " << std::endl;
  for (auto d: match_rd.detIds())
  {
    RPCDetId id(d);
    const int st(detIdToMEStation(id.station(), id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    
    auto rpcdigis = match_rd.digisInDetId(d); 
    int medianstrip(match_rd.median(rpcdigis));
    int hs = match_rd.extrapolateHsfromRPC( d, medianstrip);
    std::cout<< "RPC chamber: "<<d<<" "<<id<<" median strip:" << medianstrip <<" hs:" << hs<<std::endl; 
    for (auto p : rpcdigis)
    	std::cout << p << std::endl;
   
  }

  std::cout << "######matching LCT to Simtrack " << std::endl;
  for(auto d: match_lct.chamberIdsAllLCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    auto lcts = match_lct.allLCTsInChamber(d);
    for (auto p : lcts)    
       std::cout<<id<< p <<std::endl;
    std::cout << "-------matched lcts-------" << std::endl;
    auto lcts_matched = match_lct.lctsInChamber(d);
    for (auto q : lcts_matched)    
       std::cout<<id<< q <<std::endl;
    std::cout << "-------    end     -------" << std::endl;
    
  }


  std::cout << "######  matching Tracks to Simtrack " << std::endl;
  if (match_track.tfTracks().size()) {
     TFTrack* besttrack = match_track.bestTFTrack();
     std::cout << "       Best TFTrack                  " << std::endl;
     besttrack->print();
	 /*for (unsigned int i=0; i<triggerDigiIds.size(); i++)
	 {
	  auto id(triggerDigiIds.at(i));
          std::cout << "stub in TF DetId " << triggerDigiIds.at(i) << "  " << *triggerDigis.at(i) << std::endl;
	  std::cout << "matched stub in this Detid " << std::endl;
	  for (auto p : match_lct.lctsInChamber(id.rawId()))  std::cout << "  " << p << std::endl;
	  if ( match_lct.checkStubInChamber(id,*triggerDigis.at(i))) std::cout << "stub in TF can be matched to simtrack" << std::endl;
	  else std::cout << "stub in TF can NOT be matched to simtrack" << std::endl;

	   }*/
     std::cout << " propagated information " << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME1 << " phi " << etrk_[0].phi_propagated_ME1 << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME2 << " phi " << etrk_[0].phi_propagated_ME2 << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME3 << " phi " << etrk_[0].phi_propagated_ME3 << std::endl;
     std::cout << " propagated phi in  ME1 " << etrk_[0].phi_propagated_ME1 <<" stub phi in ME1 " <<etrk_[0].phi_ME1_TF << std::endl; 
     std::cout << " propagated phi in  ME2 " << etrk_[0].phi_interStat12 <<" stub phi in ME2 " << etrk_[0].phi_ME2_TF << std::endl; 
     std::cout << " propagated phi in  ME3 " << etrk_[0].phi_interStat23 <<" stub phi in ME3 " << etrk_[0].phi_ME3_TF << std::endl; 


  }
  else std::cout << "NO matched TFtracks"  << std::endl;


  std::cout << "==========================  end of printing ========================\n\n" << std::endl;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GEMCSCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GEMCSCAnalyzer);

