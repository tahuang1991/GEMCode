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
#include "GEMCode/GEMValidation/interface/PtassignmentHelper.h"
#include "GEMCode/GEMValidation/interface/DisplacedMuonTriggerPtassignment.h"
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

static const int NumOfTrees = 13;
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
  Float_t genGdMu_lxy;

  Float_t beamSpot_x;
  Float_t beamSpot_y;
  Float_t beamSpot_z;


  Char_t chamber_ME1_csc_sh;//bit1:odd, bit2:even
  Char_t chamber_ME2_csc_sh;
  Int_t nlayers_csc_sh_odd;
  Int_t nlayers_csc_sh_even;
  Int_t chamber_lct_odd, chamber_dg_odd, chamber_sh_odd; //
  Int_t chamber_lct_even, chamber_dg_even, chamber_sh_even; //
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
  Int_t nlayers_wg_dg_odd;
  Int_t nlayers_st_dg_odd;
  Int_t nlayers_wg_dg_even;
  Int_t nlayers_st_dg_even;

  Char_t has_rechits;
  Int_t nlayers_rechits_odd;
  Int_t nlayers_rechits_even;

  Char_t has_clct; // bit1: in odd, bit2: even
  Char_t has_alct; // bit1: in odd, bit2: even
  Char_t has_lct; // bit1: in odd, bit2: even

  Float_t phi_layer1_fit_even, phi_layer3_fit_even, phi_layer6_fit_even, phi_layer1_fit_odd, phi_layer3_fit_odd, phi_layer6_fit_odd;
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
  Float_t phi_me0_st2_odd;
  Float_t phi_me0_st2_even;
  Float_t eta_me0_st2_odd;
  Float_t eta_me0_st2_even;
  Float_t dphi_lct_odd; // dphi stored as data member in LCT
  Float_t dphi_lct_even;
  Float_t chi2_lct_odd, chi2_lct_even;
  Float_t timeErr_lct_odd, timeErr_lct_even;
  Int_t nHits_lct_odd, nHits_lct_even;
  Float_t dR_sim_lct_odd, dR_sim_lct_even;
  Float_t me0_st1_dphi_odd, me0_st1_dphi_even;
  Bool_t me0_st1_isEven_odd, me0_st1_isEven_even;
  Float_t perp_lct_odd;
  Float_t perp_lct_even;
  Bool_t passdphi_odd;
  Bool_t passdphi_even;
  Float_t fitperp_lct_odd;
  Float_t fitperp_lct_even;

  Int_t wiregroup_odd;
  Int_t wiregroup_even;
  Int_t halfstrip_odd;
  Int_t halfstrip_even;

  Int_t quality_clct_odd;
  Int_t quality_clct_even;
  Int_t quality_alct_odd;
  Int_t quality_alct_even;

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
  Float_t dphi_pad_fit_odd;
  Float_t dphi_pad_fit_even;

  Float_t phi_1strippad_odd;
  Float_t phi_1strippad_even;
  Float_t phi_2strippad_odd;
  Float_t phi_2strippad_even;
  Float_t phi_4strippad_odd;
  Float_t phi_4strippad_even;
  Float_t phi_8strippad_odd;
  Float_t phi_8strippad_even;
  Float_t dphi_1strippad_odd;
  Float_t dphi_1strippad_even;
  Float_t dphi_2strippad_odd;
  Float_t dphi_2strippad_even;
  Float_t dphi_4strippad_odd;
  Float_t dphi_4strippad_even;
  Float_t dphi_8strippad_odd;
  Float_t dphi_8strippad_even;
  Float_t dphi_1strippad_fit_odd;
  Float_t dphi_1strippad_fit_even;
  Float_t dphi_2strippad_fit_odd;
  Float_t dphi_2strippad_fit_even;
  Float_t dphi_4strippad_fit_odd;
  Float_t dphi_4strippad_fit_even;
  Float_t dphi_8strippad_fit_odd;
  Float_t dphi_8strippad_fit_even;

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
  /*
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
  */

  Bool_t allstubs_matched_TF;
  //L1Mu
  Float_t bestdRGmtCand;
  Float_t L1Mu_pt, L1Mu_eta, L1Mu_phi, L1Mu_quality, L1Mu_bx;
  UInt_t L1Mu_charge;
  Float_t L1Mu_me0_eta, L1Mu_me0_phi, L1Mu_me0_dPhi, L1Mu_me0_mindPhi1, L1Mu_me0_mindPhi2, L1Mu_me0_dR, L1Mu_me0_st1_dphi, L1Mu_me0_st2_eta, L1Mu_me0_st2_phi;
  Bool_t L1Mu_me0_st1_isEven;

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
  Int_t meRing;
  Float_t phiM_gemcsc_even,phiM_gemcsc_odd;
  Float_t phiM_st1_sh, phiM_st2_sh;
  Float_t phiM_st12_sh, phiM_st23_sh;
  Float_t phiM_fitxyz_even,phiM_fitxyz_odd;
  Float_t phiM_st1_L1, phiM_st2_L1;
  Float_t phiM_st12_L1, phiM_st23_L1;
  Float_t phiM_ME11only_L1, phiM_ME21only_L1;
  Float_t dphi_dir_ME11only_st2_L1, dphi_dir_st1_ME21only_L1, dphi_dir_ME11only_ME21only_L1;
  Float_t pt_sh_even,pt_sh_odd,ptphi_sh_even,ptphi_sh_odd,pteta_sh_even,pteta_sh_odd;
  Float_t ptphi_gemsh_odd, ptphi_gemsh_even;
  Float_t dphi_dir_st1_st2_sh, dphi_dir_st1_st2_L1, dphi_dir_st1_st2_L1_csc;
  Float_t dphi_dir_st1_st12_sh,dphi_dir_st2_st23_sh, dphi_dir_st12_st23_sh;
  Float_t dphi_dir_st1_st12_L1,dphi_dir_st2_st23_L1,dphi_dir_st1_st12_L1_csc,dphi_dir_st2_st23_L1_csc, dphi_dir_st12_st23_L1;
  Float_t hybrid_pt, hybrid_noge21_pt, position_pt, direction_ge21_pt, direction_noge21_pt;
  //Float_t ptphi_diff_sh_11,ptphi_diff_sh_12,ptphi_diff_sh_21,ptphi_diff_sh_22;
  Float_t ptphi_diff_sh;
  Float_t deltay12_fit, deltay23_fit;
  Float_t deltay12_lct, deltay23_lct;
  Float_t deltay12_sim, deltay23_sim;
  Float_t deltay12_test, deltay23_test, deltay123_test;
  Float_t eta_st2_sh, eta_st2_L1, eta_st1_L1;
  Bool_t isMe0Region;
  Int_t npar_L1;
  Float_t pt_position_sh, pt_direction_sh;
  Float_t pt_position, pt_direction_gemcsc, pt_direction_gemcsc_central, pt_direction_xfactor;
  Float_t pt_position_fit;
  Int_t npar;
  Int_t npar_lct;
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
  genGdMu_lxy = 0.0;

  beamSpot_x =0.0;
  beamSpot_y =0.0;
  beamSpot_z =0.0;


  chamber_ME1_csc_sh=0;
  chamber_ME2_csc_sh=0;
  chamber_sh_odd = -1;
  chamber_sh_even = -1;
  nlayers_csc_sh_odd = -1;
  nlayers_csc_sh_even = -1;
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


  has_csc_strips = 0;
  has_csc_wires = 0;
  has_rechits = 0;
  nlayers_wg_dg_odd = -1;
  nlayers_st_dg_odd = -1;
  nlayers_wg_dg_even = -1;
  nlayers_st_dg_even = -1;
  nlayers_rechits_odd = -1;
  nlayers_rechits_even = -1;


  phi_layer1_fit_even =-9.0;
  phi_layer1_fit_odd =-9.0;
  phi_layer3_fit_even =-9.0;
  phi_layer3_fit_odd =-9.0;
  phi_layer6_fit_even =-9.0;
  phi_layer6_fit_odd =-9.0;
  has_csc_sh = 0;
  chamber_dg_odd = -1;
  chamber_dg_even = -1;
  has_alct = 0;
  has_clct = 0;
  has_lct = 0;
  chamber_lct_odd = -1;
  chamber_lct_even = -1;
  bend_lct_odd = -9;
  bend_lct_even = -9;
  dphi_lct_odd = -9;
  dphi_lct_even = -9;
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
  phi_me0_st2_odd = -9.;
  phi_me0_st2_even = -9.;
  eta_me0_st2_odd = -9.;
  eta_me0_st2_even = -9.;
  chi2_lct_odd = -99999;
  chi2_lct_even = -99999;
  timeErr_lct_odd = -9999;
  timeErr_lct_even = -9999;
  nHits_lct_odd = 0;
  nHits_lct_even = 0;
  dR_sim_lct_odd = 10;
  dR_sim_lct_even = 10;
  passdphi_odd = false;
  passdphi_even = false;
  perp_lct_odd = -1;
  perp_lct_even = -1;
  fitperp_lct_odd = -1;
  fitperp_lct_even = -1;
  me0_st1_dphi_odd = -9;
  me0_st1_dphi_even = -9;
  me0_st1_isEven_odd = false;
  me0_st1_isEven_even = false;



  wiregroup_odd = -1;
  wiregroup_even =-1;
  halfstrip_odd =-1;
  halfstrip_even = -1;
  quality_clct_odd = -1;
  quality_clct_even = -1;
  quality_alct_odd = -1;
  quality_alct_even = -1;
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
  dphi_pad_fit_odd = -9.;
  dphi_pad_fit_even = -9.;
  deta_pad_odd = -9.;
  deta_pad_even = -9.;

  phi_1strippad_odd = -9.0;
  phi_1strippad_even = -9.0;
  phi_2strippad_odd = -9.0;
  phi_2strippad_even = -9.0;
  phi_4strippad_odd = -9.0;
  phi_4strippad_even = -9.0;
  phi_8strippad_odd = -9.0;
  phi_8strippad_even = -9.0;
  dphi_1strippad_odd = -9.0;
  dphi_1strippad_even = -9.0;
  dphi_2strippad_odd = -9.0;
  dphi_2strippad_even = -9.0;
  dphi_4strippad_odd = -9.0;
  dphi_4strippad_even = -9.0;
  dphi_8strippad_odd = -9.0;
  dphi_8strippad_even = -9.0;
  dphi_1strippad_fit_odd = -9.0;
  dphi_1strippad_fit_even = -9.0;
  dphi_2strippad_fit_odd = -9.0;
  dphi_2strippad_fit_even = -9.0;
  dphi_4strippad_fit_odd = -9.0;
  dphi_4strippad_fit_even = -9.0;
  dphi_8strippad_fit_odd = -9.0;
  dphi_8strippad_fit_even = -9.0;

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
  trackpt = -1;
  tracketa = -9;
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
  /*
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
  */

  allstubs_matched_TF = false;

  bestdRGmtCand = 99;
  L1Mu_pt = -99;
  L1Mu_eta = -99;
  L1Mu_phi = -99;
  L1Mu_quality = -99;
  L1Mu_bx = -99;
  L1Mu_charge = -99;
  L1Mu_me0_eta = -9;
  L1Mu_me0_phi = -9;
  L1Mu_me0_st2_eta = -9;
  L1Mu_me0_st2_phi = -9;
  L1Mu_me0_dPhi = -9;
  L1Mu_me0_dR = 10;
  L1Mu_me0_mindPhi1 = 10;
  L1Mu_me0_mindPhi2 = 10;
  L1Mu_me0_st1_dphi = -9;
  L1Mu_me0_st1_isEven = false;


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
  npar_lct = -1;
  meRing = -1;
  npar_L1 = -1;
  eta_st2_L1 = -9;
  eta_st1_L1 = -9;
  eta_st2_sh = -9;
  isMe0Region = false;
  ptphi_diff_sh = -9;
  phiM_st1_sh = -9;
  phiM_st2_sh = -9;
  phiM_st12_sh = -9;
  phiM_st23_sh = -9;
  phiM_st1_L1 = -9;
  phiM_st2_L1 = -9;
  phiM_st12_L1 = -9;
  phiM_st23_L1 = -9;
  phiM_ME11only_L1 = -9;
  phiM_ME21only_L1 = -9;
  dphi_dir_st1_st2_sh = -9;
  dphi_dir_st1_st2_L1 = -9;
  dphi_dir_ME11only_st2_L1 = -9;
  dphi_dir_st1_ME21only_L1 = -9;
  dphi_dir_ME11only_ME21only_L1 = -9;
  dphi_dir_st1_st2_L1_csc = -9;
  dphi_dir_st1_st12_sh = -9;
  dphi_dir_st2_st23_sh = -9;
  dphi_dir_st12_st23_sh = -9;
  dphi_dir_st1_st12_L1 = -9;
  dphi_dir_st2_st23_L1 = -9;
  dphi_dir_st1_st12_L1_csc = -9;
  dphi_dir_st2_st23_L1_csc = -9;
  dphi_dir_st12_st23_L1 = -9;
  hybrid_pt = 0.0;
  hybrid_noge21_pt = 0.0;
  position_pt = 0.0;
  direction_noge21_pt = 0.0;
  direction_ge21_pt = 0.0;
  deltay12_sim = -99;
  deltay23_sim = -99;
  deltay12_fit = -99;
  deltay23_fit = -99;
  deltay12_lct = -99;
  deltay23_lct = -99;
  deltay12_test = -99;
  deltay23_test = -99;
  deltay123_test = -99;
  pt_position_sh=-1;
  pt_position=-1;
  pt_position_fit =-1;
  pt_direction_gemcsc = -1;
  pt_direction_sh=-1;
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
  t->Branch("genGdMu_lxy", &genGdMu_lxy, "genGdMu_lxy/F");
  t->Branch("genGdMu_dR", &genGdMu_dR, "genGdMu_dR/F");


  t->Branch("chamber_ME1_csc_sh", &chamber_ME1_csc_sh);
  t->Branch("chamber_ME2_csc_sh", &chamber_ME2_csc_sh);
  t->Branch("chamber_sh_odd", &chamber_sh_odd);
  t->Branch("chamber_sh_even", &chamber_sh_even);
  t->Branch("nlayers_csc_sh_odd", &nlayers_csc_sh_odd);
  t->Branch("nlayers_csc_sh_even", &nlayers_csc_sh_even);
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

  t->Branch("chamber_dg_odd", &chamber_dg_odd);
  t->Branch("chamber_dg_even", &chamber_dg_even);
  t->Branch("has_csc_sh", &has_csc_sh);
  t->Branch("has_csc_strips", &has_csc_strips);
  t->Branch("has_csc_wires", &has_csc_wires);
  t->Branch("has_rechits", &has_rechits);
  t->Branch("nlayers_wg_dg_odd", &nlayers_wg_dg_odd);
  t->Branch("nlayers_wg_dg_even", &nlayers_wg_dg_even);
  t->Branch("nlayers_st_dg_odd", &nlayers_st_dg_odd);
  t->Branch("nlayers_st_dg_even", &nlayers_st_dg_even);
  t->Branch("nlayers_rechits_odd", &nlayers_rechits_odd);
  t->Branch("nlayers_rechits_even", &nlayers_rechits_even);


  t->Branch("has_clct", &has_clct);
  t->Branch("has_alct", &has_alct);
  t->Branch("has_lct", &has_lct);
  t->Branch("chamber_lct_odd", &chamber_lct_odd);
  t->Branch("chamber_lct_even", &chamber_lct_even);
  t->Branch("phi_layer1_fit_even", &phi_layer1_fit_even);
  t->Branch("phi_layer1_fit_odd", &phi_layer1_fit_odd);
  t->Branch("phi_layer3_fit_even", &phi_layer3_fit_even);
  t->Branch("phi_layer3_fit_odd", &phi_layer3_fit_odd);
  t->Branch("phi_layer6_fit_even", &phi_layer6_fit_even);
  t->Branch("phi_layer6_fit_odd", &phi_layer6_fit_odd);
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
  t->Branch("phi_me0_st2_odd", &phi_me0_st2_odd);
  t->Branch("phi_me0_st2_even", &phi_me0_st2_even);
  t->Branch("eta_me0_st2_odd", &eta_me0_st2_odd);
  t->Branch("eta_me0_st2_even", &eta_me0_st2_even);
  t->Branch("fitperp_lct_odd", &fitperp_lct_odd);
  t->Branch("fitperp_lct_even", &fitperp_lct_even);
  t->Branch("dphi_lct_odd", &dphi_lct_odd);
  t->Branch("dphi_lct_even", &dphi_lct_even);
  t->Branch("chi2_lct_odd", &chi2_lct_odd);
  t->Branch("chi2_lct_even", &chi2_lct_even);
  t->Branch("timeErr_lct_odd", &timeErr_lct_odd);
  t->Branch("timeErr_lct_even", &timeErr_lct_even);
  t->Branch("nHits_lct_odd", &nHits_lct_odd);
  t->Branch("nHits_lct_even", &nHits_lct_even);
  t->Branch("dR_sim_lct_odd", &dR_sim_lct_odd);
  t->Branch("dR_sim_lct_even", &dR_sim_lct_even);
  t->Branch("me0_st1_dphi_odd", &me0_st1_dphi_odd);
  t->Branch("me0_st1_dphi_even", &me0_st1_dphi_even);
  t->Branch("me0_st1_isEven_odd", &me0_st1_isEven_odd);
  t->Branch("me0_st1_isEven_even", &me0_st1_isEven_even);
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
  t->Branch("dphi_pad_fit_odd", &dphi_pad_fit_odd);
  t->Branch("dphi_pad_fit_even", &dphi_pad_fit_even);
  t->Branch("deta_pad_odd", &deta_pad_odd);
  t->Branch("deta_pad_even", &deta_pad_even);

  t->Branch("phi_1strippad_odd", &phi_1strippad_odd);
  t->Branch("phi_1strippad_even", &phi_1strippad_even);
  t->Branch("phi_2strippad_odd", &phi_2strippad_odd);
  t->Branch("phi_2strippad_even", &phi_2strippad_even);
  t->Branch("phi_4strippad_odd", &phi_4strippad_odd);
  t->Branch("phi_4strippad_even", &phi_4strippad_even);
  t->Branch("phi_8strippad_odd", &phi_8strippad_odd);
  t->Branch("phi_8strippad_even", &phi_8strippad_even);
  t->Branch("dphi_1strippad_odd", &dphi_1strippad_odd);
  t->Branch("dphi_1strippad_even", &dphi_1strippad_even);
  t->Branch("dphi_2strippad_odd", &dphi_2strippad_odd);
  t->Branch("dphi_2strippad_even", &dphi_2strippad_even);
  t->Branch("dphi_4strippad_odd", &dphi_4strippad_odd);
  t->Branch("dphi_4strippad_even", &dphi_4strippad_even);
  t->Branch("dphi_8strippad_odd", &dphi_8strippad_odd);
  t->Branch("dphi_8strippad_even", &dphi_8strippad_even);
  t->Branch("dphi_1strippad_fit_odd", &dphi_1strippad_fit_odd);
  t->Branch("dphi_1strippad_fit_even", &dphi_1strippad_fit_even);
  t->Branch("dphi_2strippad_fit_odd", &dphi_2strippad_fit_odd);
  t->Branch("dphi_2strippad_fit_even", &dphi_2strippad_fit_even);
  t->Branch("dphi_4strippad_fit_odd", &dphi_4strippad_fit_odd);
  t->Branch("dphi_4strippad_fit_even", &dphi_4strippad_fit_even);
  t->Branch("dphi_8strippad_fit_odd", &dphi_8strippad_fit_odd);
  t->Branch("dphi_8strippad_fit_even", &dphi_8strippad_fit_even);

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

  /*
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
  */

  t->Branch("allstubs_matched_TF",&allstubs_matched_TF);

  t->Branch("bestdRGmtCand", &bestdRGmtCand);
  t->Branch("L1Mu_pt", &L1Mu_pt);
  t->Branch("L1Mu_eta", &L1Mu_eta);
  t->Branch("L1Mu_phi", &L1Mu_phi);
  t->Branch("L1Mu_quality", &L1Mu_quality);
  t->Branch("L1Mu_bx", &L1Mu_bx);
  t->Branch("L1Mu_charge", &L1Mu_charge);
  t->Branch("L1Mu_me0_eta", &L1Mu_me0_eta);
  t->Branch("L1Mu_me0_phi", &L1Mu_me0_phi);
  t->Branch("L1Mu_me0_st2_eta", &L1Mu_me0_st2_eta);
  t->Branch("L1Mu_me0_st2_phi", &L1Mu_me0_st2_phi);
  t->Branch("L1Mu_me0_dPhi", &L1Mu_me0_dPhi);
  t->Branch("L1Mu_me0_dR", &L1Mu_me0_dR);
  t->Branch("L1Mu_me0_mindPhi1", &L1Mu_me0_mindPhi1);
  t->Branch("L1Mu_me0_mindPhi2", &L1Mu_me0_mindPhi2);
  t->Branch("L1Mu_me0_st1_dphi", &L1Mu_me0_st1_dphi);
  t->Branch("L1Mu_me0_st1_isEven", &L1Mu_me0_st1_isEven);

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


  t->Branch("meRing", &meRing);
  t->Branch("npar", &npar);
  t->Branch("npar_L1", &npar_L1);
  t->Branch("eta_st2_L1", &eta_st2_L1);
  t->Branch("eta_st1_L1", &eta_st1_L1);
  t->Branch("eta_st2_sh", &eta_st2_sh);
  t->Branch("isMe0Region", &isMe0Region);
  t->Branch("ptphi_gemsh_even", &ptphi_gemsh_even);
  t->Branch("ptphi_gemsh_odd", &ptphi_gemsh_odd);
  t->Branch("ptphi_diff_sh", &ptphi_diff_sh);
  t->Branch("phiM_st1_sh", &phiM_st1_sh);
  t->Branch("phiM_st2_sh", &phiM_st2_sh);
  t->Branch("phiM_st12_sh", &phiM_st12_sh);
  t->Branch("phiM_st23_sh", &phiM_st23_sh);
  t->Branch("phiM_st12_L1", &phiM_st12_L1);
  t->Branch("phiM_st23_L1", &phiM_st23_L1);
  t->Branch("phiM_st1_L1", &phiM_st1_L1);
  t->Branch("phiM_st2_L1", &phiM_st2_L1);
  t->Branch("phiM_ME11only_L1", &phiM_ME11only_L1);
  t->Branch("phiM_ME21only_L1", &phiM_ME21only_L1);
  t->Branch("dphi_dir_st1_st2_sh", &dphi_dir_st1_st2_sh);
  t->Branch("dphi_dir_st1_st2_L1", &dphi_dir_st1_st2_L1);
  t->Branch("dphi_dir_st1_st2_L1_csc", &dphi_dir_st1_st2_L1_csc);
  t->Branch("dphi_dir_ME11only_st2_L1", &dphi_dir_ME11only_st2_L1);
  t->Branch("dphi_dir_st1_ME21only_L1", &dphi_dir_st1_ME21only_L1);
  t->Branch("dphi_dir_ME11only_ME21only_L1", &dphi_dir_ME11only_ME21only_L1);
  t->Branch("dphi_dir_st1_st12_sh", &dphi_dir_st1_st12_sh);
  t->Branch("dphi_dir_st12_st23_sh", &dphi_dir_st12_st23_sh);
  t->Branch("dphi_dir_st2_st23_sh", &dphi_dir_st2_st23_sh);
  t->Branch("dphi_dir_st1_st12_L1", &dphi_dir_st1_st12_L1);
  t->Branch("dphi_dir_st2_st23_L1", &dphi_dir_st2_st23_L1);
  t->Branch("dphi_dir_st12_st23_L1", &dphi_dir_st12_st23_L1);
  t->Branch("dphi_dir_st1_st12_L1_csc", &dphi_dir_st1_st12_L1_csc);
  t->Branch("dphi_dir_st2_st23_L1_csc", &dphi_dir_st2_st23_L1_csc);
  t->Branch("position_pt", &position_pt);
  t->Branch("direction_noge21_pt", &direction_noge21_pt);
  t->Branch("direction_ge21_pt", &direction_ge21_pt);
  t->Branch("hybrid_pt", &hybrid_pt);
  t->Branch("hybrid_noge21_pt", &hybrid_noge21_pt);
  t->Branch("deltay12_sim", &deltay12_sim);
  t->Branch("deltay23_sim", &deltay23_sim);
  t->Branch("deltay12_fit", &deltay12_fit);
  t->Branch("deltay23_fit", &deltay23_fit);
  t->Branch("deltay12_lct", &deltay12_lct);
  t->Branch("deltay23_lct", &deltay23_lct);
  t->Branch("deltay12_test", &deltay12_test);
  t->Branch("deltay23_test", &deltay23_test);
  t->Branch("deltay123_test", &deltay123_test);
  t->Branch("npar_lct", &npar_lct);
  t->Branch("pt_position_sh", &pt_position_sh);
  t->Branch("pt_position", &pt_position);
  t->Branch("pt_position_fit", &pt_position_fit);
  t->Branch("pt_direction_sh", &pt_direction_sh);
  t->Branch("pt_direction_gemcsc", &pt_direction_gemcsc);
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

  enum {ME0=12};


private:

  void bookSimTracksDeltaTree();

  void analyzeTrackChamberDeltas(SimTrackMatchManager& match, int trk_no);
  void analyzeTrackEff(SimTrackMatchManager& match, int trk_no);
  void printout(SimTrackMatchManager& match, int trk_no, const char msg[300]);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);

  edm::ParameterSet cfg_;
  edm::ParameterSet displacedMuPt_cfg_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleInput_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertexInput_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackInput_;

  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> me0SimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> dtSimHitInput_;

  edm::EDGetTokenT<GEMDigiCollection> gemDigiInput_;
  edm::EDGetTokenT<GEMPadDigiCollection> gemPadDigiInput_;
  edm::EDGetTokenT<GEMCoPadDigiCollection> gemCoPadDigiInput_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHitInput_;

  edm::EDGetTokenT<ME0DigiPreRecoCollection> me0DigiInput_;
  edm::EDGetTokenT<ME0RecHitCollection> me0RecHitInput_;
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentInput_;

  edm::EDGetTokenT<CSCComparatorDigiCollection> cscComparatorDigiInput_;
  edm::EDGetTokenT<CSCWireDigiCollection> cscWireDigiInput_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clctInputs_;
  edm::EDGetTokenT<CSCALCTDigiCollection> alctInputs_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctInputs_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> mplctInputs_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHit2DInput_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentInput_;

  edm::EDGetTokenT<DTDigiCollection> dtDigiInput_;
  edm::EDGetTokenT<DTLocalTriggerCollection> dtStubInput_;
  edm::EDGetTokenT<DTRecHitCollection> dtRecHit1DPairInput_;
  edm::EDGetTokenT<DTRecSegment2DCollection> dtRecSegment2DInput_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtRecSegment4DInput_;

  edm::EDGetTokenT<RPCDigiCollection> rpcDigiInput_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitInput_;

  //edm::EDGetTokenT<L1CSCTrackCollection> cscTfTrackInputLabel_;
  //edm::EDGetTokenT<L1MuRegionalCandCollection> cscTfCandInputLabel_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> emtfTrackInputLabel_;
  edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>> regMuonCandInputLabel_;
  edm::EDGetTokenT<BXVector<l1t::Muon>> gmtInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> dtTfCandInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> rpcfTfCandInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> rpcbTfCandInputLabel_;

  edm::EDGetTokenT<L1MuRegionalCandCollection> gmtRegCandCSCInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> gmtRegCandDTInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> gmtRegCandRPCfInputLabel_;
  edm::EDGetTokenT<L1MuRegionalCandCollection> gmtRegCandRPCbInputLabel_;
  edm::EDGetTokenT<L1MuGMTCandCollection> gmtCandInputLabel_;
  edm::EDGetTokenT<l1extra::L1MuonParticleCollection> l1ExtraMuonInputLabel_;

  edm::EDGetTokenT<reco::TrackExtraCollection> recoTrackExtraInputLabel_;
  edm::EDGetTokenT<reco::TrackCollection> recoTrackInputLabel_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> recoChargedCandidateInputLabel_;


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

  TTree *tree_eff_[NumOfTrees];
  TTree *tree_delta_;

  MyTrackEff  etrk_[NumOfTrees];
  MyTrackChamberDelta dtrk_;

  int minNHitsChamberCSCSimHit_;
  int minNHitsChamberME0SimHit_;
  int minNHitsChamberCSCWireDigi_;
  int minNHitsChamberCSCStripDigi_;
  int minNHitsChamberME0Digi_;
  int minNHitsChamberME0RecHit_;
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

  /*
  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrack.getParameter<int>("verbose");
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
  */

  displacedMuPt_cfg_ = cfg_.getParameter<edm::ParameterSet>("displacedMuPtAssignment");

  auto displacedGenMu = cfg_.getParameter<edm::ParameterSet>("displacedGenMu");
  genParticleInput_ = consumes<reco::GenParticleCollection>(displacedGenMu.getParameter<edm::InputTag>("validInputTags"));

  auto simVertex = cfg_.getParameter<edm::ParameterSet>("simVertex");
  simVertexInput_ = consumes<edm::SimVertexContainer>(simVertex.getParameter<edm::InputTag>("validInputTags"));

  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrack.getParameter<int>("verbose");
  simTrackInput_ = consumes<edm::SimTrackContainer>(simTrack.getParameter<edm::InputTag>("validInputTags"));
  simTrackMinPt_ = simTrack.getParameter<double>("minPt");
  simTrackMinEta_ = simTrack.getParameter<double>("minEta");
  simTrackMaxEta_ = simTrack.getParameter<double>("maxEta");
  simTrackOnlyMuon_ = simTrack.getParameter<bool>("onlyMuon");

  auto gemSimHit_ = cfg_.getParameter<edm::ParameterSet>("gemSimHit");
  gemSimHitInput_ = consumes<edm::PSimHitContainer>(gemSimHit_.getParameter<edm::InputTag>("validInputTags"));

  auto cscSimHit_= cfg_.getParameter<edm::ParameterSet>("cscSimHit");
  cscSimHitInput_ = consumes<edm::PSimHitContainer>(cscSimHit_.getParameter<edm::InputTag>("validInputTags"));
  minNHitsChamberCSCSimHit_ = cscSimHit_.getParameter<int>("minNHitsChamber");

  auto me0SimHit_ = cfg_.getParameter<edm::ParameterSet>("me0SimHit");
  me0SimHitInput_ = consumes<edm::PSimHitContainer>(me0SimHit_.getParameter<edm::InputTag>("validInputTags"));
  minNHitsChamberME0SimHit_ = me0SimHit_.getParameter<int>("minNHitsChamber");

  auto rpcSimHit_ = cfg_.getParameter<edm::ParameterSet>("rpcSimHit");
  rpcSimHitInput_ = consumes<edm::PSimHitContainer>(rpcSimHit_.getParameter<edm::InputTag>("validInputTags"));

  auto dtSimHit_ = cfg_.getParameter<edm::ParameterSet>("dtSimHit");
  dtSimHitInput_ = consumes<edm::PSimHitContainer>(dtSimHit_.getParameter<edm::InputTag>("validInputTags"));

  auto gemDigi_= cfg_.getParameter<edm::ParameterSet>("gemStripDigi");
  gemDigiInput_ = consumes<GEMDigiCollection>(gemDigi_.getParameter<edm::InputTag>("validInputTags"));

  auto gemPad_= cfg_.getParameter<edm::ParameterSet>("gemPadDigi");
  gemPadDigiInput_ = consumes<GEMPadDigiCollection>(gemPad_.getParameter<edm::InputTag>("validInputTags"));

  auto gemCoPad_= cfg_.getParameter<edm::ParameterSet>("gemCoPadDigi");
  gemCoPadDigiInput_ = consumes<GEMCoPadDigiCollection>(gemCoPad_.getParameter<edm::InputTag>("validInputTags"));

  auto gemRecHit_= cfg_.getParameter<edm::ParameterSet>("gemRecHit");
  gemRecHitInput_ = consumes<GEMRecHitCollection>(gemRecHit_.getParameter<edm::InputTag>("validInputTags"));

  auto me0Digi_= cfg_.getParameter<edm::ParameterSet>("me0DigiPreReco");
  me0DigiInput_ = consumes<ME0DigiPreRecoCollection>(me0Digi_.getParameter<edm::InputTag>("validInputTags"));
  minNHitsChamberME0Digi_ =  me0Digi_.getParameter<int>("minNHitsChamber");

  auto me0RecHit_ = cfg_.getParameter<edm::ParameterSet>("me0RecHit");
  me0RecHitInput_ = consumes<ME0RecHitCollection>(me0RecHit_.getParameter<edm::InputTag>("validInputTags"));
  minNHitsChamberME0RecHit_ =  me0RecHit_.getParameter<int>("minNHitsChamber");

  auto me0Segment_ = cfg_.getParameter<edm::ParameterSet>("me0Segment");
  me0SegmentInput_ = consumes<ME0SegmentCollection>(me0Segment_.getParameter<edm::InputTag>("validInputTags"));

  auto cscComparatorDigi = cfg_.getParameter<edm::ParameterSet>("cscStripDigi");
  minNHitsChamberCSCStripDigi_ = cscComparatorDigi.getParameter<int>("minNHitsChamber");
  cscComparatorDigiInput_ = consumes<CSCComparatorDigiCollection>(cscComparatorDigi.getParameter<edm::InputTag>("validInputTags"));

  auto cscWireDigi = cfg_.getParameter<edm::ParameterSet>("cscWireDigi");
  minNHitsChamberCSCWireDigi_ = cscWireDigi.getParameter<int>("minNHitsChamber");
  cscWireDigiInput_ = consumes<CSCWireDigiCollection>(cscWireDigi.getParameter<edm::InputTag>("validInputTags"));

  auto cscCLCT = cfg_.getParameter<edm::ParameterSet>("cscCLCT");
  minNHitsChamberCLCT_ = cscCLCT.getParameter<int>("minNHitsChamber");
  clctInputs_ = consumes<CSCCLCTDigiCollection>(cscCLCT.getParameter<edm::InputTag>("validInputTags"));

  auto cscALCT = cfg_.getParameter<edm::ParameterSet>("cscALCT");
  minNHitsChamberALCT_ = cscALCT.getParameter<int>("minNHitsChamber");
  alctInputs_ = consumes<CSCALCTDigiCollection>(cscALCT.getParameter<edm::InputTag>("validInputTags"));

  auto cscLCT = cfg_.getParameter<edm::ParameterSet>("cscLCT");
  minNHitsChamberLCT_ = cscLCT.getParameter<int>("minNHitsChamber");
  lctInputs_ = consumes<CSCCorrelatedLCTDigiCollection>(cscLCT.getParameter<edm::InputTag>("validInputTags"));

  auto cscMPLCT = cfg_.getParameter<edm::ParameterSet>("cscMPLCT");
  minNHitsChamberMPLCT_ = cscMPLCT.getParameter<int>("minNHitsChamber");
  mplctInputs_ = consumes<CSCCorrelatedLCTDigiCollection>(cscMPLCT.getParameter<edm::InputTag>("validInputTags"));

  auto cscRecHit2D = cfg_.getParameter<edm::ParameterSet>("cscRecHit");
  cscRecHit2DInput_ = consumes<CSCRecHit2DCollection>(cscRecHit2D.getParameter<edm::InputTag>("validInputTags"));

  auto cscSegment2D = cfg_.getParameter<edm::ParameterSet>("cscSegment");
  cscSegmentInput_ = consumes<CSCSegmentCollection>(cscSegment2D.getParameter<edm::InputTag>("validInputTags"));

  auto dtDigi_= cfg_.getParameter<edm::ParameterSet>("dtDigi");
  dtDigiInput_ = consumes<DTDigiCollection>(dtDigi_.getParameter<edm::InputTag>("validInputTags"));

  auto dtStub_= cfg_.getParameter<edm::ParameterSet>("dtLocalTrigger");
  dtStubInput_ = consumes<DTLocalTriggerCollection>(dtStub_.getParameter<edm::InputTag>("validInputTags"));

  auto dtRecHit1DPair = cfg_.getParameter<edm::ParameterSet>("dtRecHit");
  dtRecHit1DPairInput_ = consumes<DTRecHitCollection>(dtRecHit1DPair.getParameter<edm::InputTag>("validInputTags"));

  auto dtSegment2D = cfg_.getParameter<edm::ParameterSet>("dtRecSegment2D");
  dtRecSegment2DInput_ = consumes<DTRecSegment2DCollection>(dtSegment2D.getParameter<edm::InputTag>("validInputTags"));

  auto dtSegment4D = cfg_.getParameter<edm::ParameterSet>("dtRecSegment4D");
  dtRecSegment4DInput_ = consumes<DTRecSegment4DCollection>(dtSegment4D.getParameter<edm::InputTag>("validInputTags"));

  auto rpcDigi_= cfg_.getParameter<edm::ParameterSet>("rpcStripDigi");
  rpcDigiInput_ = consumes<RPCDigiCollection>(rpcDigi_.getParameter<edm::InputTag>("validInputTags"));

  auto rpcRecHit_= cfg_.getParameter<edm::ParameterSet>("rpcRecHit");
  rpcRecHitInput_ = consumes<RPCRecHitCollection>(rpcRecHit_.getParameter<edm::InputTag>("validInputTags"));

  //auto tfTrack = cfg_.getParameter<edm::ParameterSet>("cscTfTrack");
  //cscTfTrackInputLabel_ = consumes<L1CSCTrackCollection>(tfTrack.getParameter<edm::InputTag>("validInputTags"));
  auto emtfTrack = cfg_.getParameter<edm::ParameterSet>("upgradeEmtfTrack");
  emtfTrackInputLabel_ = consumes<l1t::EMTFTrackCollection>(emtfTrack.getParameter<edm::InputTag>("validInputTags"));

  //auto tfCand = cfg_.getParameter<edm::ParameterSet>("cscTfCand");
  //cscTfCandInputLabel_ = consumes<L1MuRegionalCandCollection>(tfCand.getParameter<edm::InputTag>("validInputTags"));
  auto upgradeemtfCand = cfg_.getParameter<edm::ParameterSet>("upgradeEmtfCand");
  regMuonCandInputLabel_ = consumes< BXVector<l1t::RegionalMuonCand> >(upgradeemtfCand.getParameter<edm::InputTag>("validInputTags"));

  auto upgradegmt = cfg_.getParameter<edm::ParameterSet>("upgradeGMT");
  gmtInputLabel_ = consumes< BXVector<l1t::Muon> >(upgradegmt.getParameter<edm::InputTag>("validInputTags"));

  auto dtTfCand = cfg_.getParameter<edm::ParameterSet>("dtTfCand");
  dtTfCandInputLabel_ = consumes<L1MuRegionalCandCollection>(dtTfCand.getParameter<edm::InputTag>("validInputTags"));

  auto rpcfTfCand = cfg_.getParameter<edm::ParameterSet>("rpcfTfCand");
  rpcfTfCandInputLabel_ = consumes<L1MuRegionalCandCollection>(rpcfTfCand.getParameter<edm::InputTag>("validInputTags"));

  auto rpcbTfCand = cfg_.getParameter<edm::ParameterSet>("rpcbTfCand");
  rpcbTfCandInputLabel_ = consumes<L1MuRegionalCandCollection>(rpcbTfCand.getParameter<edm::InputTag>("validInputTags"));

  auto gmtRegCandCSC = cfg_.getParameter<edm::ParameterSet>("gmtRegCandCSC");
  auto gmtRegCandDT = cfg_.getParameter<edm::ParameterSet>("gmtRegCandDT");
  auto gmtRegCandRPCf = cfg_.getParameter<edm::ParameterSet>("gmtRegCandRPCf");
  auto gmtRegCandRPCb = cfg_.getParameter<edm::ParameterSet>("gmtRegCandRPCb");
  auto gmtCand = cfg_.getParameter<edm::ParameterSet>("gmtCand");
  auto l1ExtraMuonParticle = cfg_.getParameter<edm::ParameterSet>("l1ExtraMuonParticle");

  gmtRegCandCSCInputLabel_ = consumes<L1MuRegionalCandCollection>(gmtRegCandCSC.getParameter<edm::InputTag>("validInputTags"));
  gmtRegCandDTInputLabel_ = consumes<L1MuRegionalCandCollection>(gmtRegCandDT.getParameter<edm::InputTag>("validInputTags"));
  gmtRegCandRPCfInputLabel_ = consumes<L1MuRegionalCandCollection>(gmtRegCandRPCf.getParameter<edm::InputTag>("validInputTags"));
  gmtRegCandRPCbInputLabel_ = consumes<L1MuRegionalCandCollection>(gmtRegCandRPCb.getParameter<edm::InputTag>("validInputTags"));
  gmtCandInputLabel_ = consumes<L1MuGMTCandCollection>(gmtCand.getParameter<edm::InputTag>("validInputTags"));
  l1ExtraMuonInputLabel_ = consumes<l1extra::L1MuonParticleCollection>(l1ExtraMuonParticle.getParameter<edm::InputTag>("validInputTags"));

  auto recoTrackExtra = cfg_.getParameter<edm::ParameterSet>("recoTrackExtra");
  recoTrackExtraInputLabel_ = consumes<reco::TrackExtraCollection>(recoTrackExtra.getParameter<edm::InputTag>("validInputTags"));

  auto recoTrack = cfg_.getParameter<edm::ParameterSet>("recoTrack");
  recoTrackInputLabel_ = consumes<reco::TrackCollection>(recoTrack.getParameter<edm::InputTag>("validInputTags"));

  auto recoChargedCandidate = cfg_.getParameter<edm::ParameterSet>("recoChargedCandidate");
  recoChargedCandidateInputLabel_ = consumes<reco::RecoChargedCandidateCollection>(recoChargedCandidate.getParameter<edm::InputTag>("validInputTags"));


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
  cscStationsCo_.push_back(std::make_pair(0,1));
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
  ev.getByToken(simTrackInput_, sim_tracks);
  const edm::SimTrackContainer & sim_track = *sim_tracks.product();

  edm::Handle<edm::SimVertexContainer> sim_vertices;
  ev.getByToken(simVertexInput_, sim_vertices);
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
    //SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);
    SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es,
                               genParticleInput_,
                               simVertexInput_,
                               simTrackInput_,
                               gemSimHitInput_,
                               cscSimHitInput_,
                               rpcSimHitInput_,
                               me0SimHitInput_,
                               dtSimHitInput_,
                               gemDigiInput_,
                               gemPadDigiInput_,
                               gemCoPadDigiInput_,
                               gemRecHitInput_,
                               me0DigiInput_,
                               me0RecHitInput_,
                               me0SegmentInput_,
                               cscComparatorDigiInput_,
                               cscWireDigiInput_,
                               clctInputs_,
                               alctInputs_,
                               lctInputs_,
                               mplctInputs_,
                               cscRecHit2DInput_,
                               cscSegmentInput_,
                               dtDigiInput_,
                               dtStubInput_,
                               dtRecHit1DPairInput_,
                               dtRecSegment2DInput_,
                               dtRecSegment4DInput_,
                               rpcDigiInput_,
                               rpcRecHitInput_,
                               //cscTfTrackInputLabel_,
                               //cscTfCandInputLabel_,
                               emtfTrackInputLabel_,
                               regMuonCandInputLabel_,
                               gmtInputLabel_,
                               dtTfCandInputLabel_,
                               rpcfTfCandInputLabel_,
                               rpcbTfCandInputLabel_,
                               gmtRegCandCSCInputLabel_,
                               gmtRegCandDTInputLabel_,
                               gmtRegCandRPCfInputLabel_,
                               gmtRegCandRPCbInputLabel_,
                               gmtCandInputLabel_,
                               l1ExtraMuonInputLabel_,
                               recoTrackExtraInputLabel_,
                               recoTrackInputLabel_,
                               recoChargedCandidateInputLabel_
			       );

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
  const ME0DigiMatcher& match_me0digi = match.me0Digis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  const ME0RecHitMatcher& match_me0rh = match.me0RecHits();
  //const L1TrackMatcher& match_track = match.l1Tracks();
  const UpgradeL1TrackMatcher& match_track = match.l1Tracks();
  //const L1GlobalMuonTriggerMatcher& match_l1_gmt = match.l1GMTCands();
  //const HLTTrackMatcher& match_hlt_track = match.hltTracks();
  const SimTrack &t = match_sh.trk();


  /*
  auto matchedDarkBoson(match_gen.getMatchedDarkBoson());
  auto matchedGENMuon(match_gen.getMatchedGENMuon());
  if (match_gen.checkRunOK()){
  	std::cout <<" matched dark boson mass "<< matchedDarkBoson->mass() << std::endl;
 	std::cout <<"matched muon pt "<< matchedGENMuon->pt()<<" phi "<< matchedGENMuon->phi()<<" eta "<< matchedGENMuon->eta() <<" dR "<< match_gen.matchedGenMudR()<<" dxy "<< match_gen.matchedGenMudxy() <<std::endl;
  }
  std::cout <<"Sim trk_no " << trk_no <<" eta "<< t.momentum().eta() << " phi "<< t.momentum().phi() << " pt "<< t.momentum().pt()<<" pz "<<  t.momentum().pz()<<std::endl;
  */

  if (verbose_) std::cout <<"GEMCSCAnalyzer step1 "<< std::endl;
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
    auto matchedDarkBoson(match_gen.getMatchedDarkBoson());
    auto matchedGENMuon(match_gen.getMatchedGENMuon());
    if (matchedDarkBoson){
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
    }
    if (matchedGENMuon){
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
	 etrk_[s].genGdMu_lxy = std::sqrt(etrk_[s].genGdMu_vx*etrk_[s].genGdMu_vx + etrk_[s].genGdMu_vy*etrk_[s].genGdMu_vy);
	 etrk_[s].genGdMu_dxy = match_gen.matchedGenMudxy();
	 etrk_[s].genGdMu_dR = match_gen.matchedGenMudR();
    }
  }
  }


  int chargesign = (t.charge()>0? 1:0);
  float pt = t.momentum().pt();
  // SimHits
  auto csc_simhits(match_sh.chamberIdsCSC(0));
  GlobalPoint gp_sh_odd[NumOfTrees];
  GlobalPoint gp_sh_even[NumOfTrees];
  GlobalVector gv_sh_odd[NumOfTrees];
  GlobalVector gv_sh_even[NumOfTrees];
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
    //std::cout <<"CSCid with simhits "<< id << " nlayer "<< nlayers << std::endl;
    if (nlayers < minNHitsChamberCSCSimHit_) continue;
    GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
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

    if (odd) etrk_[st].chamber_sh_odd = id.chamber();
    else etrk_[st].chamber_sh_even = id.chamber();
    if (odd) etrk_[st].nlayers_csc_sh_odd = nlayers;
    else etrk_[st].nlayers_csc_sh_even = nlayers;

    if (odd) gv_sh_odd[st] = ym;
    else gv_sh_even[st] = ym;
    GlobalPoint keygp(match_sh.simHitPositionKeyLayer(id));
    if (odd) gp_sh_odd[st] = keygp;
    else gp_sh_even[st] = keygp;
    if (odd) etrk_[st].eta_cscsh_odd = keygp.eta();
    else     etrk_[st].eta_cscsh_even = keygp.eta();
    if (odd) etrk_[st].phi_cscsh_odd = keygp.phi();
    else     etrk_[st].phi_cscsh_even = keygp.phi();
    if (odd) etrk_[st].perp_cscsh_odd = keygp.perp();
    else     etrk_[st].perp_cscsh_even = keygp.perp();
    if (odd) etrk_[st].dphipositionpt_cscsh_odd = deltaPhi(float(keygp.phi()),etrk_[st].ptphi_sh_odd);
    else     etrk_[st].dphipositionpt_cscsh_even = deltaPhi(float(keygp.phi()),etrk_[st].ptphi_sh_even);

    if (st==2 or st==3){
      	if (odd) etrk_[1].eta_cscsh_odd = keygp.eta();
      	else     etrk_[1].eta_cscsh_even = keygp.eta();
      	if (odd) etrk_[1].phi_cscsh_odd = keygp.phi();
      	else     etrk_[1].phi_cscsh_even = keygp.phi();
      	if (odd) etrk_[1].perp_cscsh_odd = keygp.perp();
      	else     etrk_[1].perp_cscsh_even = keygp.perp();
    	if (odd) etrk_[1].dphipositionpt_cscsh_odd = deltaPhi(float(keygp.phi()),etrk_[st].ptphi_sh_odd);
    	else     etrk_[1].dphipositionpt_cscsh_even = deltaPhi(float(keygp.phi()),etrk_[st].ptphi_sh_even);
	}

    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_csc_sh |= 1;
      else etrk_[1].has_csc_sh |= 2;

      if (odd) etrk_[1].nlayers_csc_sh_odd = nlayers;
      else etrk_[1].nlayers_csc_sh_even = nlayers;

      if (odd) etrk_[st].chamber_sh_odd = id.chamber();
      else etrk_[st].chamber_sh_even = id.chamber();
      if (odd) gp_sh_odd[1] = keygp;
      else gp_sh_even[1] = keygp;
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

  //ME0 simhits
  auto me0_simhits(match_sh.superChamberIdsME0());
  std::cout <<"me0 simthits , chamber id size "<< me0_simhits.size() << std::endl;
  for (auto d : me0_simhits){
    const ME0DetId id(d);
    int nlayers = match_sh.nLayersWithHitsInSuperChamber(d);
    std::cout <<"ME0 Detid "<< id <<" nlayer hits "<< nlayers << std::endl;
    if (nlayers < minNHitsChamberME0SimHit_) continue;
    bool odd(id.chamber()%2 == 1);
    //ME0DetId(int region, int layer,int chamber, int roll)
    GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInSuperChamber(d));
    GlobalPoint keygp(match_sh.simHitsMeanPosition(match_sh.hitsInSuperChamber(d)));
    //etrk_[ME0].bending_sh = match_sh.LocalBendingInChamber(d);
    const ME0DetId id1(id.region(), 1, id.chamber(), 0);
    const ME0DetId id6(id.region(), 6, id.chamber(), 0);
    GlobalPoint gp1 = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id1.rawId()));
    GlobalPoint gp6 = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id6.rawId()));
    if (odd)  etrk_[ME0].nlayers_csc_sh_odd = nlayers;
    else  etrk_[ME0].nlayers_csc_sh_even = nlayers;
    if (odd) etrk_[ME0].has_csc_sh |= 1;
    else etrk_[ME0].has_csc_sh |= 2;
    if (odd) etrk_[ME0].chamber_sh_odd = id.chamber();
    else etrk_[ME0].chamber_sh_even = id.chamber();
    if (odd) etrk_[ME0].eta_cscsh_odd = keygp.eta();
    else     etrk_[ME0].eta_cscsh_even = keygp.eta();
    if (odd) etrk_[ME0].phi_cscsh_odd = keygp.phi();
    else     etrk_[ME0].phi_cscsh_even = keygp.phi();
    if (odd) etrk_[ME0].perp_cscsh_odd = keygp.perp();
    if (odd) etrk_[ME0].perp_cscsh_even = keygp.perp();
    if (odd){
    	etrk_[ME0].pteta_sh_odd = ym.eta();
    	etrk_[ME0].ptphi_sh_odd = ym.phi();
    	etrk_[ME0].pt_sh_odd = ym.perp();
    }else{
    	etrk_[ME0].pteta_sh_even = ym.eta();
    	etrk_[ME0].ptphi_sh_even = ym.phi();
    	etrk_[ME0].pt_sh_even = ym.perp();
    }
    if (odd){
    	if (match_sh.hitsInChamber(id1.rawId()).size()>0){
		etrk_[ME0].eta_layer1_sh_odd = gp1.eta();
		etrk_[ME0].phi_layer1_sh_odd = gp1.phi();
		etrk_[ME0].perp_layer1_sh_odd = gp1.perp();
		//std::cout <<"layer1 id "<< id1 <<" phi "<< gp1.phi() <<" eta "<< gp1.eta() <<" x "<< gp1.x()<<" y "<< gp1.y()<<" z "<< gp1.z() <<" perp "<< gp1.perp() << std::endl;
		}
    	if (match_sh.hitsInChamber(id6.rawId()).size()>0){
		etrk_[ME0].eta_layer6_sh_odd = gp6.eta();
		etrk_[ME0].phi_layer6_sh_odd = gp6.phi();
		etrk_[ME0].perp_layer6_sh_odd = gp6.perp();
		//std::cout <<"layer6 id "<< id6 <<" phi "<< gp6.phi() <<" eta "<< gp6.eta() <<" x "<< gp6.x()<<" y "<< gp6.y()<<" z "<< gp6.z() <<" perp "<< gp6.perp() << std::endl;
		}
    }else{
    	if (match_sh.hitsInChamber(id1.rawId()).size()>0){
		etrk_[ME0].eta_layer1_sh_even = gp1.eta();
		etrk_[ME0].phi_layer1_sh_even = gp1.phi();
		etrk_[ME0].perp_layer1_sh_even = gp1.perp();
		//std::cout <<"layer1 id "<< id1 <<" phi "<< gp1.phi() <<" eta "<< gp1.eta() <<" x "<< gp1.x()<<" y "<< gp1.y()<<" z "<< gp1.z() <<" perp "<< gp1.perp() << std::endl;
		}
    	if (match_sh.hitsInChamber(id6.rawId()).size()>0){
		etrk_[ME0].eta_layer6_sh_even = gp6.eta();
		etrk_[ME0].phi_layer6_sh_even = gp6.phi();
		etrk_[ME0].perp_layer6_sh_even = gp6.perp();
		//std::cout <<"layer6 id "<< id6 <<" phi "<< gp6.phi() <<" eta "<< gp6.eta() <<" x "<< gp6.x()<<" y "<< gp6.y()<<" z "<< gp6.z() <<" perp "<< gp6.perp() << std::endl;
		}
    }

  }
  if (verbose_) std::cout <<"GEMCSCAnalyzer step2 "<< std::endl;
  // CSC strip digis
  for(auto d: match_cd.chamberIdsStrip(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const int nlayers(match_cd.nLayersWithStripInChamber(d));
    //std::cout <<"CSC strip digi, CSCid "<< id <<" nlayer "<< nlayers << std::endl;
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
    //std::cout <<"CSC wire digi, CSCid "<< id <<" nlayer "<< nlayers << std::endl;
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

    //std::cout <<"CSC CLCT, CSCid "<< id <<" clct "<< clct << std::endl;
    if (odd) etrk_[st].chamber_dg_odd = id.chamber();
    else etrk_[st].chamber_dg_even = id.chamber();
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

    //std::cout <<"CSC ALCT, CSCid "<< id <<" alct "<< alct << std::endl;
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

  if (verbose_) std::cout <<"GEMCSCAnalyzer step3 "<< std::endl;
  // holders for track's LCTs
  Digi lct_odd[12];
  Digi lct_even[12];
  GlobalPoint gp_lct_odd[12];
  GlobalPoint gp_lct_even[12];
  GlobalPoint gp_fit_odd[12];
  GlobalPoint gp_fit_even[12];
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

    //std::cout <<"CSC LCT, CSCid "<< id <<" lct "<< lct <<" gp eta "<< gp.eta()<<" phi "<< gp.phi() <<" perp "<< gp.perp()<<std::endl;
    if (csc_simhits.find(d) == csc_simhits.end())
	std::cout <<"failed to find simhits in CSCid "<< id <<" but find LCT" << lct << std::endl;
    if (odd)
    {
      lct_odd[st] = lct;
      gp_lct_odd[st] = gp;
      etrk_[st].bend_lct_odd = bend;
      etrk_[st].phi_lct_odd = gp.phi();
      etrk_[st].eta_lct_odd = gp.eta();
      etrk_[st].perp_lct_odd = gp.perp();
      if (fabs(etrk_[st].perp_lct_odd-etrk_[st].perp_cscsh_odd)>5.0 and etrk_[st].perp_cscsh_odd>10.0)
	  std::cout <<"CSCid "<< id <<" perp_cscsh_odd "<< etrk_[st].perp_cscsh_odd<<" perp_lct_odd "<< etrk_[st].perp_lct_odd<<" csc_phi "<< etrk_[st].phi_cscsh_odd<<" phi_lct_odd "<< etrk_[st].phi_lct_odd << std::endl;

      etrk_[st].dphi_lct_odd = digi_dphi(lct);
      etrk_[st].bx_lct_odd = digi_bx(lct);
      etrk_[st].hs_lct_odd = digi_channel(lct);
      etrk_[st].wg_lct_odd = digi_wg(lct);
      etrk_[st].chamber_lct_odd = id.chamber() ;
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
      if (fabs(etrk_[st].perp_lct_even-etrk_[st].perp_cscsh_even)>5.0 and etrk_[st].perp_cscsh_even>10.0)
	  std::cout <<"CSCid "<< id <<" perp_cscsh_even "<< etrk_[st].perp_cscsh_even <<" perp_lct_even "<< etrk_[st].perp_lct_even <<" csc_phi "<< etrk_[st].phi_cscsh_even <<" phi_lct_even "<< etrk_[st].phi_lct_even << std::endl;
      etrk_[st].dphi_lct_even = digi_dphi(lct);
      etrk_[st].bx_lct_even = digi_bx(lct);
      etrk_[st].hs_lct_even = digi_channel(lct);
      etrk_[st].wg_lct_even = digi_wg(lct);
      etrk_[st].chamber_lct_even = id.chamber();
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
        etrk_[1].perp_lct_odd = gp.perp();
        etrk_[1].dphi_lct_odd = digi_dphi(lct);
        etrk_[1].bx_lct_odd = digi_bx(lct);
        etrk_[1].hs_lct_odd = digi_channel(lct);
        etrk_[1].wg_lct_odd = digi_wg(lct);
        etrk_[1].chamber_lct_odd = id.chamber();
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
        etrk_[1].chamber_lct_even = id.chamber();
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
    GlobalPoint lctgp(match_lct.getGlobalPosition(id.rawId(), bestMatchingLCT));
    //std::cout <<"csc id "<< id <<" lctgp x "<< lctgp.x() <<"  y "<< lctgp.y() << " phi "<< lctgp.phi() <<" size of gps "<< gps.size()<< std::endl;
    if (gps.size()>=3){
      for (auto gp: gps){
	if (gp.z() > 0)
		zs.push_back(gp.z());
	else zs.push_back(-gp.z());
	xs.push_back(gp.x());
	ys.push_back(gp.y());
	ezs.push_back(0);
	float gpphi = gp.phi();
	if (phis.size()>0 and gpphi>0 and phis[0]<0 and  (gpphi-phis[0])>3.1416)
	    phis.push_back(gpphi-2*3.1415926);
	else if (phis.size()>0 and gpphi<0 and phis[0]>0 and (gpphi-phis[0])<-3.1416)
	    phis.push_back(gpphi+2*3.1415926);
	else
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


    float alpha = -99., beta = 0.;
    PtassignmentHelper::calculateAlphaBeta(zs, phis, ezs, ephis, status,
		     alpha, beta);
    //float alphax = -99., betax = 0.;
    //PtassignmentHelper::calculateAlphaBeta(zs, xs, ezs, exs, status,
    //		     alphax, betax);
    //float alphay = -99., betay = 0.;
    //PtassignmentHelper::calculateAlphaBeta(zs, ys, ezs, exs, status,
    //		     alphay, betay);
    if (phis.size() < 3 or fabs(alpha)>=99){
	std::cout <<"warning, falied to fit comparator digis, cscid "<< id <<",num of digis: "<< phis.size()<<" alpha "<< alpha <<" beta "<< beta << std::endl;
	alpha = lctgp.phi();
	beta = 0.0;
    }
    //ME11. even layer1 581.98; odd layer1 611.38, d=2.2
    //ME12. even layer1 678.706, layer6 691.406; odd layer1 706.106, layer6 718.806, d=2.54
    //ME21. even layer1 809.506; odd layer1 834.306 , d=2.54
    //ME22, even layer1 809.506, layer6 822.206; odd layer1 834.306, layer6 847.006, d=2.54
    //Z(layern) = Z(layer1)+d*(n-1)
    if(odd){
	if (id.station()==1 and (id.ring()==1 or id.ring()==4)){
	    etrk_[1].phi_layer1_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 1));
	    etrk_[1].phi_layer3_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 3));
	    etrk_[1].phi_layer6_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 6));
	    //etrk_[1].phiM_fitxyz_odd = atan(betay/betax);
	    gp_fit_odd[st] = GlobalPoint(GlobalPoint::Cylindrical(lctgp.perp(), alpha+beta*match_lct.zpositionOfLayer(d, 3), lctgp.z()));
	    gp_fit_odd[1] = GlobalPoint(GlobalPoint::Cylindrical(lctgp.perp(), alpha+beta*match_lct.zpositionOfLayer(d, 3), lctgp.z()));
	    //std::cout <<((etrk_[1].has_lct&1)>0?" odd ":" not odd ") <<std::endl;
	}
	etrk_[st].phi_layer1_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 1));
	etrk_[st].phi_layer3_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 3));
	etrk_[st].phi_layer6_fit_odd = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 6));
	//etrk_[st].phiM_fitxyz_odd = atan(betay/betax);
	gp_fit_odd[st] = GlobalPoint(GlobalPoint::Cylindrical(lctgp.perp(), alpha+beta*match_lct.zpositionOfLayer(d, 3), lctgp.z()));
        //std::cout <<"cscid "<< id <<" ring "<< id.ring() <<" st "<< st <<" alpha "<< alpha <<" beta "<< beta <<" phi layer1 "<< etrk_[st].phi_layer1_fit_odd<<" layer6 "<<etrk_[st].phi_layer6_fit_odd<<" gp x "<< gp_fit_odd[st].x()<<" y "<< gp_fit_odd[st].y() <<std::endl;
    }else{
	if (id.station()==1 and (id.ring()==1 or id.ring() ==4)){
	    etrk_[1].phi_layer1_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 1));
	    etrk_[1].phi_layer3_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 3));
	    etrk_[1].phi_layer6_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 6));
	    //etrk_[1].phiM_fitxyz_even = atan(betay/betax);
	    gp_fit_even[1] = GlobalPoint(GlobalPoint::Cylindrical(lctgp.perp(), alpha+beta*match_lct.zpositionOfLayer(d, 3), lctgp.z()));
	    //std::cout <<((etrk_[1].has_lct&2)>0?" even ":" not even ") <<std::endl;
	}
	etrk_[st].phi_layer1_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 1));
	etrk_[st].phi_layer3_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 3));
	etrk_[st].phi_layer6_fit_even = PtassignmentHelper::normalizePhi(alpha+beta*match_lct.zpositionOfLayer(d, 6));
	//etrk_[st].phiM_fitxyz_even = atan(betay/betax);
	gp_fit_even[st] = GlobalPoint(GlobalPoint::Cylindrical(lctgp.perp(), alpha+beta*match_lct.zpositionOfLayer(d, 3), lctgp.z()));
        //std::cout <<"cscid "<< id <<" ring "<< id.ring() <<" st "<< st <<" alpha "<< alpha <<" beta "<< beta <<" phi layer1 "<< etrk_[st].phi_layer1_fit_even<<" layer6 "<<etrk_[st].phi_layer6_fit_even<<" gp x "<< gp_fit_even[st].x()<<" y "<< gp_fit_even[st].y() <<std::endl;
    }
    if (gps.size()>=3 and id.ring()==1 and id.station()==2 and fabs(etrk_[st].phi_layer3_fit_even-etrk_[st].phi_lct_even)>0.3 and fabs(etrk_[st].phi_lct_even)>3 and fabs(etrk_[st].phi_lct_even)<4){
      std::cout <<"id "<< id <<" phi from fit "<< etrk_[st].phi_layer3_fit_even <<" phi from simhits "<< etrk_[st].phi_cscsh_even
	  <<" phi from lct "<< etrk_[st].phi_lct_even << " fitting alpha "<< alpha <<" beta "<<beta <<" zposition "<< match_lct.zpositionOfLayer(d, 3) << std::endl;
      int igp=0;
      for (auto phi: phis){
	  std::cout <<" igp "<< igp <<" phi "<< phi << std::endl;
	  igp++;
      }
    }
  }


  if (verbose_) std::cout <<"GEMCSCAnalyzer step4 "<< std::endl;
   //for GEMs in station1, it will be also filled in ME11
  // GEM simhits in superchamber
  GlobalPoint gp_gemsh_odd[12];
  GlobalPoint gp_gemsh_even[12];
  for(auto d: match_sh.superChamberIdsGEM())
  {
    GEMDetId id(d);
    int MEStation = id.station();

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
	if (odd) gp_gemsh_odd[st] = keygp;
	else gp_gemsh_even[st] = keygp;
	//std::cout <<"Matching GEMHits GEM Id "<< id <<" gp.eta "<< keygp.eta() <<" gp.phi "<< keygp.phi() << std::endl;
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
	int st_csc = (st==3 and (etrk_[st].has_csc_sh&3)==0)? 1:st;
      	if (odd and etrk_[st_csc].phi_cscsh_odd>-9){
	    etrk_[st].dphi_sh_odd = deltaPhi(etrk_[st_csc].phi_cscsh_odd, float(keygp.phi()));// here should we use st or st_csc?
    	    etrk_[st].dphipositionpt_gemsh_odd = deltaPhi(float(keygp.phi()),etrk_[st_csc].ptphi_sh_odd);
            phi_m =  PtassignmentHelper::PhiMomentum_Radius(etrk_[st].dphi_sh_odd, etrk_[st].phi_gemsh_odd, etrk_[st_csc].perp_cscsh_odd, etrk_[st].perp_gemsh_odd);
	    etrk_[st].phiM_gemcsc_odd = phi_m;
	    //std::cout <<"dphi(M,P) "<< deltaPhi(etrk_[st].phi_gemsh_odd, gv_sh_odd[st_csc].phi())<<" momentum phi from simhits "<< gv_sh_odd[st_csc].phi() <<" phi from gem-csc "<< phi_m << std::endl;
	}
      	else if (not(odd) and etrk_[st_csc].phi_cscsh_even>-9) {
	    etrk_[st].dphi_sh_even = deltaPhi(etrk_[st_csc].phi_cscsh_even, float(keygp.phi()));
    	    etrk_[st].dphipositionpt_gemsh_even = deltaPhi(float(keygp.phi()),etrk_[st_csc].ptphi_sh_even);
	    phi_m =  PtassignmentHelper::PhiMomentum_Radius(etrk_[st].dphi_sh_even, etrk_[st].phi_gemsh_even, etrk_[st_csc].perp_cscsh_even, etrk_[st].perp_gemsh_even);
	    etrk_[st].phiM_gemcsc_even = phi_m;
	    //std::cout <<"dphi(M,P) "<< deltaPhi(etrk_[st].phi_gemsh_even, gv_sh_even[st_csc].phi())<<" momentum phi from simhits "<< gv_sh_even[st_csc].phi() <<" phi from gem-csc "<< phi_m << std::endl;
	}
	if (st==2 or st==3){
		if (odd) gp_gemsh_odd[1] = keygp;
		else gp_gemsh_even[1] = keygp;
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
      		if (odd and etrk_[1].phi_cscsh_odd>-9) etrk_[1].dphi_sh_odd = deltaPhi(etrk_[1].phi_cscsh_odd,float(keygp.phi()));
      		else if (etrk_[1].phi_cscsh_even>-9)     etrk_[1].dphi_sh_even = deltaPhi(etrk_[1].phi_cscsh_even,float(keygp.phi()));
		if (odd) etrk_[1].phiM_gemcsc_odd = phi_m;
		else etrk_[1].phiM_gemcsc_even = phi_m;
		if (odd) etrk_[1].dphipositionpt_gemsh_odd = deltaPhi(float(keygp.phi()),etrk_[1].ptphi_sh_odd);
		else etrk_[1].dphipositionpt_gemsh_even = deltaPhi(float(keygp.phi()),etrk_[1].ptphi_sh_even);
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

  if (verbose_) std::cout <<"GEMCSCAnalyzer step5 "<< std::endl;
  // placeholders for best mtching pads
  GlobalPoint best_pad_odd[12];
  GlobalPoint best_pad_even[12];

  // GEM digis and pads in superchambers
  for(auto d: match_gd.superChamberIdsDigi())
  {
    GEMDetId id(d);
    int MEStation = id.station();

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

	GlobalPoint keygp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id_tmp.rawId()));
        auto bestgem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, keygp);
      if (odd)
      {
        best_pad_odd[st] = bestgem_dg_and_gp.second;
        etrk_[st].has_gem_pad |= 1;
        etrk_[st].chamber_lct_odd = id.chamber();
        etrk_[st].pad_odd = digi_channel(bestgem_dg_and_gp.first);
        etrk_[st].hsfromgem_odd = match_gd.extrapolateHsfromGEMPad( d, digi_channel(bestgem_dg_and_gp.first));
        etrk_[st].phi_pad_odd = best_pad_odd[st].phi();
        etrk_[st].eta_pad_odd = best_pad_odd[st].eta();
	//question: ME1a is not included here
        if (is_valid(lct_odd[st]))
        {
	    auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_odd[st]);
	    best_pad_odd[st] = gem_dg_and_gp.second;
	    etrk_[st].bx_pad_odd = digi_bx(gem_dg_and_gp.first);
	    etrk_[st].dphi_pad_odd = deltaPhi(etrk_[st].phi_lct_odd, etrk_[st].phi_pad_odd);
	    etrk_[st].dphi_pad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, etrk_[st].phi_pad_odd);
	    //std::cout <<"GEMid "<< id << std::endl;
	    //std::cout <<"GEM sim phi "<< keygp.phi() <<" pad phi "<< etrk_[st].phi_pad_odd <<" phi_lct "<< etrk_[st].phi_lct_odd <<" phi_fit_lct "<< etrk_[st].phi_layer3_fit_odd<< std::endl;
	    etrk_[st].deta_pad_odd = etrk_[st].eta_lct_odd - etrk_[st].eta_pad_odd;
      	}
    	}
      else
      {
        best_pad_even[st] = bestgem_dg_and_gp.second;
      	etrk_[st].has_gem_pad |= 2;
        etrk_[st].pad_even = digi_channel(bestgem_dg_and_gp.first);
        etrk_[st].hsfromgem_even = match_gd.extrapolateHsfromGEMPad( d, digi_channel(bestgem_dg_and_gp.first));
        etrk_[st].phi_pad_even = best_pad_even[st].phi();
        etrk_[st].eta_pad_even = best_pad_even[st].eta();
        if (is_valid(lct_even[st]))
        {
	    auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_even[st]);
	    best_pad_even[st] = gem_dg_and_gp.second;
	    etrk_[st].bx_pad_even = digi_bx(gem_dg_and_gp.first);
	    etrk_[st].dphi_pad_even = deltaPhi(etrk_[st].phi_lct_even, etrk_[st].phi_pad_even);
	    etrk_[st].dphi_pad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, etrk_[st].phi_pad_even);
	    //std::cout <<"GEMid "<< id << std::endl;
	    //std::cout <<"GEM sim phi "<< keygp.phi() <<" pad phi "<< etrk_[st].phi_pad_even <<" phi_lct "<< etrk_[st].phi_lct_even <<" phi_fit_lct "<< etrk_[st].phi_layer3_fit_even<< std::endl;
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
    int MEStation = id.station();

    const int stations(detIdToMEStation(MEStation,id.ring()));
    int st;
    if (stations==2 or stations==3) st=1;//ME1a,ME1b->ME11
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
	GlobalPoint keygp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(id_tmp.rawId()));
        auto bestgem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, keygp);

      if (odd)
      {
        best_pad_odd[st] = bestgem_dg_and_gp.second;
      	etrk_[st].has_gem_pad |= 1;
        etrk_[st].pad_odd = digi_channel(bestgem_dg_and_gp.first);
        etrk_[st].hsfromgem_odd = match_gd.extrapolateHsfromGEMPad( d, digi_channel(bestgem_dg_and_gp.first));
        etrk_[st].phi_pad_odd = best_pad_odd[st].phi();
        etrk_[st].eta_pad_odd = best_pad_odd[st].eta();
      	if (is_valid(lct_odd[st]))
      	{
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_odd[st]);
        	best_pad_odd[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_odd = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].dphi_pad_odd = deltaPhi(etrk_[st].phi_lct_odd, etrk_[st].phi_pad_odd);
        	etrk_[st].dphi_pad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, etrk_[st].phi_pad_odd);
        	etrk_[st].deta_pad_odd = etrk_[st].eta_lct_odd - etrk_[st].eta_pad_odd;
      	}
    	}
      else
      {
        best_pad_even[st] = bestgem_dg_and_gp.second;
      	etrk_[st].has_gem_pad |= 2;
        etrk_[st].pad_even = digi_channel(bestgem_dg_and_gp.first);
        etrk_[st].hsfromgem_even = match_gd.extrapolateHsfromGEMPad( d, digi_channel(bestgem_dg_and_gp.first));
        etrk_[st].phi_pad_even = best_pad_even[st].phi();
        etrk_[st].eta_pad_even = best_pad_even[st].eta();
      	if (is_valid(lct_even[st]))
      	{
        	auto gem_dg_and_gp = match_gd.digiInGEMClosestToCSC(pads, gp_lct_even[st]);
        	best_pad_even[st] = gem_dg_and_gp.second;
        	etrk_[st].bx_pad_even = digi_bx(gem_dg_and_gp.first);
        	etrk_[st].dphi_pad_even = deltaPhi(etrk_[st].phi_lct_even, etrk_[st].phi_pad_even);
        	etrk_[st].dphi_pad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, etrk_[st].phi_pad_even);
        	etrk_[st].deta_pad_even = etrk_[st].eta_lct_even - etrk_[st].eta_pad_even;
      	}
    	}
      if (id_tmp.layer()==1) break;
    }
   }

  if (verbose_) std::cout <<"GEMCSCAnalyzer step6 "<< std::endl;
  for(auto d: match_gd.superChamberIdsCoPad())
  {
    GEMDetId id(d);
    int MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_gem_copad |= 1;
    else     etrk_[st].has_gem_copad |= 2;

    auto copads = match_gd.coPadsInSuperChamber(d);
    if (copads.size() == 0) continue;
    if (odd) etrk_[st].Copad_odd = digi_channel(copads.at(0));
    else etrk_[st].Copad_even = digi_channel(copads.at(0));
    std::cout <<"Matching GEMCopad detid "<< id <<" size "<< copads.size() << std::endl;

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

  //for GE11 and GE21, different phi segmenta
  // GEM digis and pads in superchambers
  for(auto d: match_gd.detIdsDigi())
  {
    GEMDetId id(d);
    int MEStation = id.station();

    const int st(detIdToMEStation(MEStation,id.ring()));
    bool odd(id.chamber()%2==1);
    auto GPsPad1 = match_gd.positionPad1InDetId(d);
    auto GPsPad2 = match_gd.positionPad2InDetId(d);
    auto GPsPad4 = match_gd.positionPad4InDetId(d);
    auto GPsPad8 = match_gd.positionPad8InDetId(d);
    GEMDetId id_tmp(id.region(), id.ring(), id.station(), id.layer(), id.chamber(), 0);
    if (GPsPad1.size()==0) continue;
    float dphi1 = 99;
    float dphi2 = 99;
    float dphi4 = 99;
    float dphi8 = 99;

    float lct_phi = -9;
    float lct_fit_phi = -9;
    float phi_pad = -9;
    float phi_gemsh = -9;
    if (odd ){
	if (is_valid(lct_odd[st])){
	    lct_phi = etrk_[st].phi_lct_odd;
	    lct_fit_phi = etrk_[st].phi_layer3_fit_odd;
	}
	phi_pad = etrk_[st].phi_pad_odd;
	phi_gemsh = etrk_[st].phi_gemsh_odd;
    }else {
	if (is_valid(lct_even[st])){
	    lct_phi = etrk_[st].phi_lct_even;
	    lct_fit_phi = etrk_[st].phi_layer3_fit_even;
	}
	phi_pad = etrk_[st].phi_pad_odd;
	phi_gemsh = etrk_[st].phi_gemsh_even;
    }
    GlobalPoint gp1, gp2, gp4, gp8;
    //std::cout <<"GEM sim phi "<< keygp.phi() <<" 1strippad phi "<< GPsPad1.front().phi()<<" 2strippad phi "<< GPsPad2.front().phi()<<" 4strippad phi "<< GPsPad4.front().phi() <<" 8strippad phi "<< GPsPad8.front().phi() << std::endl;
    for (auto gp : GPsPad1){
    	if (fabs(gp.phi()-phi_gemsh) < dphi1){
	    dphi1 = fabs(gp.phi() - phi_gemsh);
	    gp1 = GlobalPoint(gp);
	}
    }
    for (auto gp : GPsPad2){
    	if (fabs(gp.phi()-phi_gemsh) < dphi2){
	    dphi2 = fabs(gp.phi() - phi_gemsh);
	    gp2 = GlobalPoint(gp);
	}
    }
    for (auto gp : GPsPad4){
    	if (fabs(gp.phi()-phi_gemsh) < dphi4){
	    dphi4 = fabs(gp.phi() - phi_gemsh);
	    gp4 = GlobalPoint(gp);
	}
    }
    for (auto gp : GPsPad8){
    	if (fabs(gp.phi()-phi_gemsh) < dphi8){
	    dphi8 = fabs(gp.phi() - phi_gemsh);
	    gp8 = GlobalPoint(gp);
	}
    }
    if (odd){
 	    etrk_[st].phi_1strippad_odd = gp1.phi();
 	    etrk_[st].phi_2strippad_odd = gp2.phi();
 	    etrk_[st].phi_4strippad_odd = gp4.phi();
 	    etrk_[st].phi_8strippad_odd = gp8.phi();
      	    if (is_valid(lct_odd[st])){
	    	etrk_[st].dphi_1strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp1.phi()));
	    	etrk_[st].dphi_1strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp1.phi()));
	    	etrk_[st].dphi_2strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp2.phi()));
	    	etrk_[st].dphi_2strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp2.phi()));
	    	etrk_[st].dphi_4strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp4.phi()));
	    	etrk_[st].dphi_4strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp4.phi()));
	    	etrk_[st].dphi_8strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp8.phi()));
	    	etrk_[st].dphi_8strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp8.phi()));
	    }
    }else{
 	    etrk_[st].phi_1strippad_even = gp1.phi();
 	    etrk_[st].phi_2strippad_even = gp2.phi();
 	    etrk_[st].phi_4strippad_even = gp4.phi();
 	    etrk_[st].phi_8strippad_even = gp8.phi();
      	    if (is_valid(lct_even[st])){
	    	etrk_[st].dphi_1strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp1.phi()));
	    	etrk_[st].dphi_1strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp1.phi()));
	    	etrk_[st].dphi_2strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp2.phi()));
	    	etrk_[st].dphi_2strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp2.phi()));
	    	etrk_[st].dphi_4strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp4.phi()));
	    	etrk_[st].dphi_4strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp4.phi()));
	    	etrk_[st].dphi_8strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp8.phi()));
	    	etrk_[st].dphi_8strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp8.phi()));
	    }
    }

    if (odd and MEStation==1 ){
 	    etrk_[1].phi_1strippad_odd = gp1.phi();
 	    etrk_[1].phi_2strippad_odd = gp2.phi();
 	    etrk_[1].phi_4strippad_odd = gp4.phi();
 	    etrk_[1].phi_8strippad_odd = gp8.phi();
      	    if (is_valid(lct_odd[st])){
	    	etrk_[1].dphi_1strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp1.phi()));
	    	etrk_[1].dphi_1strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp1.phi()));
	    	etrk_[1].dphi_2strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp2.phi()));
	    	etrk_[1].dphi_2strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp2.phi()));
	    	etrk_[1].dphi_4strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp4.phi()));
	    	etrk_[1].dphi_4strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp4.phi()));
	    	etrk_[1].dphi_8strippad_odd = deltaPhi(etrk_[st].phi_lct_odd, float(gp8.phi()));
	    	etrk_[1].dphi_8strippad_fit_odd = deltaPhi(etrk_[st].phi_layer3_fit_odd, float(gp8.phi()));
	    }
    }else if(MEStation==1){
 	    etrk_[1].phi_1strippad_even = gp1.phi();
 	    etrk_[1].phi_2strippad_even = gp2.phi();
 	    etrk_[1].phi_4strippad_even = gp4.phi();
 	    etrk_[1].phi_8strippad_even = gp8.phi();
      	    if (is_valid(lct_even[st])){
	    	etrk_[1].dphi_1strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp1.phi()));
	    	etrk_[1].dphi_1strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp1.phi()));
	    	etrk_[1].dphi_2strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp2.phi()));
	    	etrk_[1].dphi_2strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp2.phi()));
	    	etrk_[1].dphi_4strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp4.phi()));
	    	etrk_[1].dphi_4strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp4.phi()));
	    	etrk_[1].dphi_8strippad_even = deltaPhi(etrk_[st].phi_lct_even, float(gp8.phi()));
	    	etrk_[1].dphi_8strippad_fit_even = deltaPhi(etrk_[st].phi_layer3_fit_even, float(gp8.phi()));
	    }
    }

    if (verbose_){
	std::cout <<"GEM id "<< id << std::endl;
	std::cout <<" final GEM sim phi "<< phi_gemsh <<" 1strippad "<< gp1.phi()<<" 2strippad "<< gp2.phi()<<" 4strippad "<< gp4.phi()<<" 8strippad "<< gp8.phi() <<" pad phi "<< phi_pad <<" lct phi "<< lct_phi <<" lct fit phi "<< lct_fit_phi << std::endl;
    }

  }

  //ME0 digis
  auto me0digis(match_me0digi.superChamberIds());
  std::cout <<"me0 digis , chamber id size "<< me0digis.size() << std::endl;
  for (auto d: me0digis){
    const ME0DetId id(d);
    int nlayers = match_me0digi.nLayersWithDigisInSuperChamber(d);
    std::cout <<"ME0 Detid "<< id <<" nlayer hits "<< nlayers << std::endl;
    if (nlayers < minNHitsChamberME0Digi_) continue;
    bool odd(id.chamber()%2 == 1);
    if (odd) etrk_[ME0].has_csc_strips |= 1;
    else etrk_[ME0].has_csc_strips |= 2;
    if (odd) etrk_[ME0].nlayers_st_dg_odd = nlayers;
    else etrk_[ME0].nlayers_st_dg_even = nlayers;
  }


  //ME0 Rechits
  /*auto me0RecHitsChamber(match_me0rh.chamberIdsME0RecHit());
  //std::cout <<"me0 rechits , chamber id size "<< me0RecHitsChamber.size() << std::endl;
  for (auto d: me0RecHitsChamber){
    const ME0DetId id(d);
    auto me0RecHits(match_me0rh.me0RecHitsInChamber(d));
    std::cout <<"ME0 Detid "<< id<<" matched Rechits num "<< me0RecHits.size() << std::endl;
    //std::cout << "Available rechits in chamber" << std::endl;
    //for (auto rh: me0RecHits){
      //std::cout << rh.me0Id() << " " << rh << std::endl;
    //}
  }*/
  auto me0RecHitsSuperChamber(match_me0rh.superChamberIdsME0RecHit());
  //std::cout <<"me0 rechits , superchamber id size "<< me0RecHitsSuperChamber.size() << std::endl;
  for (auto d: me0RecHitsSuperChamber){
    const ME0DetId id(d);
    int nlayers = match_me0rh.nLayersWithRecHitsInSuperChamber(d);
    std::cout <<"ME0 Detid "<< id<<" nlayers with rechits " << nlayers << std::endl; 
    if (nlayers < minNHitsChamberME0RecHit_) continue;
    bool odd(id.chamber()%2 == 1);
    if (odd) etrk_[ME0].has_rechits |= 1;
    else etrk_[ME0].has_rechits |= 2;
    if (odd) etrk_[ME0].nlayers_rechits_odd = nlayers;
    else etrk_[ME0].nlayers_rechits_even = nlayers;
  }

  //ME0 Segments
  std::vector<ME0Segment> allmatchedSegments;
  auto me0Segs(match_me0rh.superChamberIdsME0Segment_bydR());
  std::cout <<"me0 Segs , chamber id size "<< me0Segs.size() << std::endl;
  for (auto d: me0Segs){
    const ME0DetId id(d);
    std::cout <<"ME0 Detid "<< id << std::endl;
    auto allSegs(match_me0rh.me0SegmentsInSuperChamber_bydR(d));
    allmatchedSegments.insert(allmatchedSegments.begin(), allSegs.begin(), allSegs.end());
    auto me0Segment(match_me0rh.bestME0Segment_bydR(d));
    if (fabs(me0Segment.chi2()) < 0.0000001 and fabs(me0Segment.time()) < 0.00000001 and fabs(me0Segment.timeErr()) < 0.000001 and fabs(me0Segment.deltaPhi())<.0001 and me0Segment.specificRecHits().size()<3){
      std::cout <<"has valid ME0 id "<< id <<" but empty ME0 segment, allSegs size  "<< allSegs.size() << std::endl;
      for (auto  seg : allSegs)
	  std::cout <<" matched seg , detid "<< seg.me0DetId() <<" "<< seg << std::endl;
      continue;
    }
    
    float dPhi = match_me0rh.me0DeltaPhi(me0Segment);
    double chi2 = me0Segment.chi2(); float time = me0Segment.time();
    float timeErr = me0Segment.timeErr(); int nRecHits = me0Segment.nRecHits();
    bool odd(id.chamber()%2 == 1);
    GlobalPoint gp(match_me0rh.globalPoint(me0Segment));
    GlobalPoint gp_propagated(match_me0rh.propagateToZ(AVERAGE_ME0_Z*id.region()));
    GlobalPoint gp_ME0_st2(match_me0rh.propagateFromME0ToCSC(me0Segment, t.momentum().pt(), 2, odd)); 
    float dR = deltaR(float(gp.eta()), float(gp.phi()), float(gp_propagated.eta()), float(gp_propagated.phi()));
    cout << "me0Segment " << me0Segment <<" gp eta "<< gp.eta() <<" phi "<< gp.phi() << " dphi  "<<  me0Segment.deltaPhi() <<" dphi in matcher "<< dPhi << endl;
    std::cout <<"Sim Pt "<< t.momentum().pt() <<" gp_ME0_st2 "<<  gp_ME0_st2 <<" eta "<< gp_ME0_st2.eta() <<" phi  "<< gp_ME0_st2.phi() << std::endl;

    if (odd) etrk_[ME0].chamber_lct_odd = id.chamber();
    else etrk_[ME0].chamber_lct_even = id.chamber();
    if (odd) etrk_[ME0].has_lct |= 1;
    else etrk_[ME0].has_lct |= 2;

    if (me0_simhits.find(d) == me0_simhits.end())
	std::cout <<"failed to find simhits in ME0 "<< id <<" but find Me0segment " << me0Segment << std::endl;
    if (odd)
    {
      etrk_[ME0].phi_lct_odd = gp.phi();
      etrk_[ME0].eta_lct_odd = gp.eta();
      etrk_[ME0].perp_lct_odd = gp.perp();
      etrk_[ME0].phi_me0_st2_odd = gp_ME0_st2.phi();
      etrk_[ME0].eta_me0_st2_odd = gp_ME0_st2.eta();
      //if (fabs(etrk_[ME0].perp_lct_odd-etrk_[ME0].perp_cscsh_odd)>5.0 and etrk_[ME0].perp_cscsh_odd>10.0)
	  //std::cout <<"CSCid "<< id <<" perp_cscsh_odd "<< etrk_[ME0].perp_cscsh_odd<<" perp_lct_odd "<< etrk_[ME0].perp_lct_odd<<" csc_phi "<< etrk_[ME0].phi_cscsh_odd<<" phi_lct_odd "<< etrk_[ME0].phi_lct_odd << std::endl;
      etrk_[ME0].dphi_lct_odd = dPhi;
      etrk_[ME0].chi2_lct_odd = chi2;
      etrk_[ME0].bx_lct_odd = time;
      etrk_[ME0].timeErr_lct_odd = timeErr;
      etrk_[ME0].nHits_lct_odd = nRecHits;
      etrk_[ME0].chamber_lct_odd = id.chamber() ;
      etrk_[ME0].dR_sim_lct_odd = dR;
      //deltaPhi(lct, me0)
      if (fabs(etrk_[1].phi_lct_even) < 4 and abs(etrk_[1].chamber_lct_even/2- id.chamber()) <= 1){
	  etrk_[ME0].me0_st1_dphi_odd = deltaPhi(etrk_[1].phi_lct_even, gp.phi());
	  etrk_[ME0].me0_st1_isEven_odd = true;
      }else if (fabs(etrk_[1].phi_lct_odd) < 4 and abs(etrk_[1].chamber_lct_odd/2- id.chamber()) <= 1){
	  etrk_[ME0].me0_st1_dphi_odd = deltaPhi(etrk_[1].phi_lct_odd, gp.phi());
	  etrk_[ME0].me0_st1_isEven_odd = false;
      }else etrk_[ME0].me0_st1_dphi_odd = -9;
    }
    else
    {
      etrk_[ME0].phi_lct_even = gp.phi();
      etrk_[ME0].eta_lct_even = gp.eta();
      etrk_[ME0].perp_lct_even = gp.perp();
      etrk_[ME0].phi_me0_st2_even = gp_ME0_st2.phi();
      etrk_[ME0].eta_me0_st2_even = gp_ME0_st2.eta();
      //if (fabs(etrk_[ME0].perp_lct_even-etrk_[ME0].perp_cscsh_even)>5.0 and etrk_[ME0].perp_cscsh_even>10.0)
	  //std::cout <<"CSCid "<< id <<" perp_cscsh_even "<< etrk_[ME0].perp_cscsh_even <<" perp_lct_even "<< etrk_[ME0].perp_lct_even <<" csc_phi "<< etrk_[ME0].phi_cscsh_even <<" phi_lct_even "<< etrk_[ME0].phi_lct_even << std::endl;
      etrk_[ME0].dphi_lct_even = dPhi;
      etrk_[ME0].chi2_lct_even = chi2;
      etrk_[ME0].timeErr_lct_even = timeErr;
      etrk_[ME0].nHits_lct_even = nRecHits;
      etrk_[ME0].bx_lct_even = time;
      etrk_[ME0].chamber_lct_even = id.chamber() ;
      etrk_[ME0].dR_sim_lct_even = dR;
      if (fabs(etrk_[1].phi_lct_even) < 4 and abs(etrk_[1].chamber_lct_even/2- id.chamber()) <= 1){
	  etrk_[ME0].me0_st1_dphi_even = deltaPhi(etrk_[1].phi_lct_even, gp.phi());
	  etrk_[ME0].me0_st1_isEven_even = true;
      }else if (fabs(etrk_[1].phi_lct_odd) < 4 and abs(etrk_[1].chamber_lct_odd/2- id.chamber()) <= 1){
	  etrk_[ME0].me0_st1_dphi_even = deltaPhi(etrk_[1].phi_lct_odd, gp.phi());
	  etrk_[ME0].me0_st1_isEven_even = false;
      }else etrk_[ME0].me0_st1_dphi_even = -9;
    }
  }

  if (verbose_) std::cout <<"GEMCSCAnalyzer step7 "<< std::endl;
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

  if (verbose_) std::cout <<"GEMCSCAnalyzer step8 "<< std::endl;
  GlobalPoint gp1,gp2, gp3, gp4, gp_ge11, gp_ge21;
  GlobalVector gv1,gv2;
  if (etrk_[1].has_csc_sh>0 and etrk_[6].has_csc_sh>0){
     etrk_[0].meRing =1;
     if ((etrk_[1].has_csc_sh&1)>0 and (etrk_[6].has_csc_sh&2)>0){
        gp1=gp_sh_odd[1];
        gp2=gp_sh_even[6];
        gv1=gv_sh_odd[1];
        gv2=gv_sh_even[6];
	if ((etrk_[1].has_gem_sh&1)>0)
	    gp_ge11 = gp_gemsh_odd[1];
	if ((etrk_[6].has_gem_sh&2)>0)
	    gp_ge21 = gp_gemsh_even[6];
	if ((etrk_[8].has_csc_sh&2)>0)
	    gp3=gp_sh_even[8];
     }else if ((etrk_[1].has_csc_sh&1)>0 and (etrk_[6].has_csc_sh&1)>0 ){
        gp1=gp_sh_odd[1];
        gp2=gp_sh_odd[6];
        gv1=gv_sh_odd[1];
        gv2=gv_sh_odd[6];
	if ((etrk_[1].has_gem_sh&1)>0)
	    gp_ge11 = gp_gemsh_odd[1];
	if ((etrk_[6].has_gem_sh&1)>0)
	    gp_ge21 = gp_gemsh_odd[6];
	if ((etrk_[8].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[8];
     }else if ((etrk_[1].has_csc_sh&2)>0 and (etrk_[6].has_csc_sh&2)>0 ){
        gp1=gp_sh_even[1];
        gp2=gp_sh_even[6];
        gv1=gv_sh_even[1];
        gv2=gv_sh_even[6];
	if ((etrk_[1].has_gem_sh&2)>0)
	    gp_ge11 = gp_gemsh_even[1];
	if ((etrk_[6].has_gem_sh&2)>0)
	    gp_ge21 = gp_gemsh_even[6];
	if ((etrk_[8].has_csc_sh&2)>0)
	    gp3=gp_sh_even[8];
     }else if ((etrk_[1].has_csc_sh&2)>0 and (etrk_[6].has_csc_sh&1)>0 ){
        gp1=gp_sh_even[1];
        gp2=gp_sh_odd[6];
        gv1=gv_sh_even[1];
        gv2=gv_sh_odd[6];
	if ((etrk_[1].has_gem_sh&2)>0)
	    gp_ge11 = gp_gemsh_even[1];
	if ((etrk_[6].has_gem_sh&1)>0)
	    gp_ge21 = gp_gemsh_odd[6];
	if ((etrk_[8].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[8];
     }
  }
     //low eta region , ME12, ME22
  if (etrk_[4].has_csc_sh>0 and etrk_[7].has_csc_sh>0){
     etrk_[0].meRing =2;
     if ((etrk_[4].has_csc_sh&1)>0 and (etrk_[7].has_csc_sh&2)>0){
        gp1=gp_sh_odd[4];
        gp2=gp_sh_even[7];
        gv1=gv_sh_odd[4];
        gv2=gv_sh_even[7];
	if ((etrk_[9].has_csc_sh&2)>0)
 	    gp3=gp_sh_even[9];
     }else if ((etrk_[4].has_csc_sh&1)>0 and (etrk_[7].has_csc_sh&1)>0){
        gp1=gp_sh_odd[4];
        gp2=gp_sh_odd[7];
        gv1=gv_sh_odd[4];
        gv2=gv_sh_odd[7];
	if ((etrk_[9].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[9];
     }else if ((etrk_[4].has_csc_sh&2)>0 and (etrk_[7].has_csc_sh&2)>0){
        gp1=gp_sh_even[4];
        gp2=gp_sh_even[7];
        gv1=gv_sh_even[4];
        gv2=gv_sh_even[7];
	if ((etrk_[9].has_csc_sh&2)>0)
	    gp3=gp_sh_even[9];
     }else if ((etrk_[4].has_csc_sh&2)>0 and (etrk_[7].has_csc_sh&1)>0){
        gp1=gp_sh_even[4];
        gp2=gp_sh_odd[7];
        gv1=gv_sh_even[4];
        gv2=gv_sh_odd[7];
	if ((etrk_[9].has_csc_sh&1)>0)
	    gp3=gp_sh_odd[9];
     }
  }


  DisplacedMuonTriggerPtassignment displacedMuonL1Pt(match_lct.allLctsMatched2SimMuon(), 
	  					     //match_gd.allGempadsMatch2SimMuon_2strip(), 
						     gemPadDigiInput_,
						     me0SegmentInput_,
						     displacedMuPt_cfg_, 
						     match_lct.eventSetup(), 
						     match_lct.event());
  displacedMuonL1Pt.setCharge(etrk_[0].charge);
  if (displacedMuonL1Pt.getNParity() >= 0)  std::cout <<"DisplacedMuon get npar "<< displacedMuonL1Pt.getNParity()<<" ring "<<displacedMuonL1Pt.getMeRing() << std::endl;


  int npar =displacedMuonL1Pt.getNParity();
  DisplacedMuonTriggerPtassignment displacedMuonL1Pt_sim(gp1,
	  						 gp2, 
							 gp3, 
							 gp4, 
							 gp_ge11, 
							 gp_ge21, 
							 npar, 
							 displacedMuPt_cfg_,
							 match_lct.eventSetup(), 
							 match_lct.event());
  displacedMuonL1Pt_sim.setCharge(etrk_[0].charge);

  if ((etrk_[1].has_csc_sh || etrk_[4].has_csc_sh) and (etrk_[6].has_csc_sh || etrk_[7].has_csc_sh)) {
     if ((etrk_[8].has_csc_sh || etrk_[9].has_csc_sh) || (etrk_[10].has_csc_sh || etrk_[11].has_csc_sh))
     	etrk_[0].hasSt3orSt4_sh=true;
     if (displacedMuonL1Pt.getNParity()>=0 and (fabs(displacedMuonL1Pt.getRadiusSt(1)-gp1.perp())>.02*gp1.perp()
		 or fabs(displacedMuonL1Pt.getRadiusSt(2)-gp2.perp())>.02*gp2.perp()
		 or fabs(displacedMuonL1Pt.getRadiusSt(3)-gp3.perp())>.02*gp3.perp()))
	    std::cout <<" warning, difference between fitting and sim is large, module, npar  "<< displacedMuonL1Pt.getNParity()<<" ring "<< displacedMuonL1Pt.getMeRing() << std::endl;
     if (etrk_[0].meRing == 1 and displacedMuonL1Pt.getNParity()>=0 and displacedMuonL1Pt.runDirectionbased(true)){
     	etrk_[0].phiM_st1_L1 = displacedMuonL1Pt.getlocalPhiDirection(1);
     	etrk_[0].phiM_st2_L1 = displacedMuonL1Pt.getlocalPhiDirection(2);
     	etrk_[0].phiM_st12_L1 = displacedMuonL1Pt.getlocalPhiDirection(12);
     	etrk_[0].phiM_st23_L1 = displacedMuonL1Pt.getlocalPhiDirection(23);
	etrk_[0].dphi_dir_st1_st2_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 2);
	etrk_[0].dphi_dir_st1_st12_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 12);
	etrk_[0].dphi_dir_st2_st23_L1 = displacedMuonL1Pt.getdeltaPhiDirection(2, 23);
	etrk_[0].dphi_dir_st12_st23_L1 = displacedMuonL1Pt.getdeltaPhiDirection(12, 23);
	etrk_[0].direction_ge21_pt = displacedMuonL1Pt.getDirectionPt();
	if (displacedMuonL1Pt.checkME0Region() and displacedMuonL1Pt.runDirectionbased(true, false)){
	   etrk_[0].phiM_ME11only_L1 = displacedMuonL1Pt.getlocalPhiDirection(1);
	   etrk_[0].dphi_dir_ME11only_st2_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 2);
	}
	//std::cout <<"Direction based ge21 pt "<< etrk_[0].direction_ge21_pt  << std::endl;
     }

     if (displacedMuonL1Pt.getNParity()>=0 and displacedMuonL1Pt.runDirectionbased(false)){
     	etrk_[0].phiM_ME21only_L1 = displacedMuonL1Pt.getlocalPhiDirection(2);
     	etrk_[0].phiM_st12_L1 = displacedMuonL1Pt.getlocalPhiDirection(12);
     	etrk_[0].phiM_st23_L1 = displacedMuonL1Pt.getlocalPhiDirection(23);
	if (etrk_[0].meRing == 2){
	    etrk_[0].phiM_st1_L1 = displacedMuonL1Pt.getlocalPhiDirection(1);
	    etrk_[0].dphi_dir_st1_st2_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 2);
	    etrk_[0].dphi_dir_st1_st12_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 12);
	    etrk_[0].dphi_dir_st2_st23_L1 = displacedMuonL1Pt.getdeltaPhiDirection(2, 23);
	    etrk_[0].dphi_dir_st12_st23_L1 = displacedMuonL1Pt.getdeltaPhiDirection(12, 23);
	}

	//assign values to these variables both in ring1 and ring2
	etrk_[0].dphi_dir_st1_st2_L1_csc = displacedMuonL1Pt.getdeltaPhiDirection(1, 2);
	etrk_[0].dphi_dir_st1_st12_L1_csc = displacedMuonL1Pt.getdeltaPhiDirection(1, 12);
	etrk_[0].dphi_dir_st2_st23_L1_csc = displacedMuonL1Pt.getdeltaPhiDirection(2, 23);
	etrk_[0].dphi_dir_st12_st23_L1 = displacedMuonL1Pt.getdeltaPhiDirection(12, 23);

	etrk_[0].direction_noge21_pt = displacedMuonL1Pt.getDirectionPt();
	//std::cout <<"Direction based no ge21 pt "<< etrk_[0].direction_noge21_pt  << std::endl;
	
	if (etrk_[0].meRing == 1 and displacedMuonL1Pt.checkME0Region() and displacedMuonL1Pt.runDirectionbased(false, false)){
	   etrk_[0].dphi_dir_ME11only_ME21only_L1 = displacedMuonL1Pt.getdeltaPhiDirection(1, 2);
	}

     }

     //etrk_[0].pt_direction_sh=PtassignmentHelper::Ptassign_Direction(csc_bending_angle_12, gp2.eta(), npar);
     etrk_[0].eta_st2_L1 = displacedMuonL1Pt.getTrackEta();
     etrk_[0].eta_st1_L1 = displacedMuonL1Pt.getCSCEta(1);//get station 1 eta
     etrk_[0].isMe0Region = displacedMuonL1Pt.checkME0Region();
     etrk_[0].npar_L1 = displacedMuonL1Pt.getNParity();
     if (etrk_[8].has_csc_sh || etrk_[9].has_csc_sh){
     	etrk_[0].deltay12_sim = PtassignmentHelper::deltaYcalculation(gp1, gp2);
     	etrk_[0].deltay23_sim = -PtassignmentHelper::deltaYcalculation(gp3, gp2);
	etrk_[0].pt_position_sh=PtassignmentHelper::Ptassign_Position_gp(gp1, gp2, gp3, gp2.eta(), npar); //t.momentum().eta()
	if (verbose_){
		std::cout <<"GEMCSCAnalyer sim gp1 x"<< gp1.x()<<" y "<< gp1.y() <<" gp2 x "<< gp2.x()<<" y "<< gp2.y()<<" phi "<< gp2.phi()<< std::endl;
		std::cout <<"at Sim deltay12 "<<  etrk_[0].deltay12_sim<<" deltay23 "<<  etrk_[0].deltay23_sim << std::endl;
	}
     }

  }


  if (etrk_[1].has_lct>0 and etrk_[6].has_lct>0 and etrk_[8].has_lct>0){
     int npar=-1;
     GlobalPoint gp1,gp2, gp3;
     GlobalPoint gp1_fit, gp2_fit, gp3_fit;
     if ((etrk_[1].has_lct&1)>0 and (etrk_[6].has_lct&2)>0 and (etrk_[8].has_lct&2)>0){
        gp1=gp_lct_odd[1];
        gp2=gp_lct_even[6];
        gp3=gp_lct_even[8];
        gp1_fit=gp_fit_odd[1];
        gp2_fit=gp_fit_even[6];
        gp3_fit=gp_fit_even[8];
	npar=0;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==1){
		etrk_[1].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[6].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[8].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 1 "<< std::endl;

     }else if ((etrk_[1].has_lct&1)>0 and (etrk_[6].has_lct&1)>0 and (etrk_[8].has_lct&1)>0){
        gp1=gp_lct_odd[1];
        gp2=gp_lct_odd[6];
        gp3=gp_lct_odd[8];
        gp1_fit=gp_fit_odd[1];
        gp2_fit=gp_fit_odd[6];
        gp3_fit=gp_fit_odd[8];
	npar=1;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==1){
		etrk_[1].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[6].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[8].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 1 "<< std::endl;
     }else if ((etrk_[1].has_lct&2)>0 and (etrk_[6].has_lct&2)>0 and (etrk_[8].has_lct&2)>0){
        gp1=gp_lct_even[1];
        gp2=gp_lct_even[6];
        gp3=gp_lct_even[8];
        gp1_fit=gp_fit_even[1];
        gp2_fit=gp_fit_even[6];
        gp3_fit=gp_fit_even[8];
	npar=2;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==1){
		etrk_[1].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[6].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[8].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 1 "<< std::endl;
     }else if ((etrk_[1].has_lct&2)>0 and (etrk_[6].has_lct&1)>0 and (etrk_[8].has_lct&1)>0){
        gp1=gp_lct_even[1];
        gp2=gp_lct_odd[6];
        gp3=gp_lct_odd[8];
        gp1_fit=gp_fit_even[1];
        gp2_fit=gp_fit_odd[6];
        gp3_fit=gp_fit_odd[8];
	npar=3;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==1){
		etrk_[1].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[6].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[8].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 1 "<< std::endl;
     }
     etrk_[0].hasSt1St2St3=true;
     etrk_[0].npar_lct= npar;
     etrk_[0].deltay12_lct = PtassignmentHelper::deltaYcalculation(gp1, gp2);
     etrk_[0].deltay23_lct = -PtassignmentHelper::deltaYcalculation(gp3, gp2);
     etrk_[0].deltay12_fit = PtassignmentHelper::deltaYcalculation(gp1_fit, gp2_fit);
     etrk_[0].deltay23_fit = -PtassignmentHelper::deltaYcalculation(gp3_fit, gp2_fit);
     etrk_[0].pt_position=PtassignmentHelper::Ptassign_Position_gp(gp1, gp2, gp3, gp2.eta(), npar);
     etrk_[0].pt_position_fit=PtassignmentHelper::Ptassign_Position_gp(gp1_fit, gp2_fit, gp3_fit, gp2_fit.eta(), npar);
     if (verbose_){
	std::cout <<"LCT gp1 x "<< gp1.x()<<" y "<< gp1.y() <<" phi "<< gp1.phi()<<" perp "<< gp1.perp() <<" gp2 x "<< gp2.x()<<" y "<< gp2.y()<<" phi "<< gp2.phi()<<" perp "<< gp2.perp() << std::endl;
     	std::cout <<"at L1 deltay12 "<<  etrk_[0].deltay12_lct <<" deltay23 "<<  etrk_[0].deltay23_lct << std::endl;
	std::cout <<"fit gp1 x "<< gp1_fit.x()<<" y "<< gp1_fit.y() <<" phi "<< gp1_fit.phi() <<" perp "<< gp1_fit.perp() <<" gp2 x "<< gp2_fit.x()<<" y "<< gp2_fit.y()<<" phi "<< gp2_fit.phi()<<" perp "<< gp2_fit.perp() << std::endl;
     	std::cout <<"at digi with fit deltay12 "<<  etrk_[0].deltay12_fit <<" deltay23 "<<  etrk_[0].deltay23_fit << std::endl;
     }

  }

  if (etrk_[4].has_lct>0 and etrk_[7].has_lct>0 and etrk_[9].has_lct>0){
     int npar=-1;
     GlobalPoint gp1,gp2, gp3;
     GlobalPoint gp1_fit, gp2_fit, gp3_fit;
     if ((etrk_[4].has_lct&1)>0 and (etrk_[7].has_lct&2)>0 and (etrk_[9].has_lct&2)>0){
        gp1=gp_lct_odd[4];
        gp2=gp_lct_even[7];
        gp3=gp_lct_even[9];
        gp1_fit=gp_fit_odd[4];
        gp2_fit=gp_fit_even[7];
        gp3_fit=gp_fit_even[9];
	npar=0;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==2){
		etrk_[4].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[7].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[9].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 2 "<< std::endl;
     }else if ((etrk_[4].has_lct&1)>0 and (etrk_[7].has_lct&1)>0 and (etrk_[9].has_lct&1)>0){
        gp1=gp_lct_odd[4];
        gp2=gp_lct_odd[7];
        gp3=gp_lct_odd[9];
        gp1_fit=gp_fit_odd[4];
        gp2_fit=gp_fit_odd[7];
        gp3_fit=gp_fit_odd[9];
	npar=1;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==2){
		etrk_[4].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[7].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[9].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 2 "<< std::endl;
     }else if ((etrk_[4].has_lct&2)>0 and (etrk_[7].has_lct&2)>0 and (etrk_[9].has_lct&2)>0){
        gp1=gp_lct_even[4];
        gp2=gp_lct_even[7];
        gp3=gp_lct_even[9];
        gp1_fit=gp_fit_even[4];
        gp2_fit=gp_fit_even[7];
        gp3_fit=gp_fit_even[9];
	npar=2;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==2){
		etrk_[4].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[7].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[9].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 2 "<< std::endl;
     }else if ((etrk_[4].has_lct&2)>0 and (etrk_[7].has_lct&1)>0 and (etrk_[9].has_lct&1)>0){
        gp1=gp_lct_even[4];
        gp2=gp_lct_odd[7];
        gp3=gp_lct_odd[9];
        gp1_fit=gp_fit_even[4];
        gp2_fit=gp_fit_odd[7];
        gp3_fit=gp_fit_odd[9];
	npar=3;
	if (displacedMuonL1Pt.getNParity()==npar and displacedMuonL1Pt.getMeRing()==2){
		etrk_[4].fitperp_lct_even = displacedMuonL1Pt.getRadiusSt(1);
		etrk_[7].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(2);
		etrk_[9].fitperp_lct_odd = displacedMuonL1Pt.getRadiusSt(3);
	}else
	    std::cout <<"displaced muon pt assignment, npar  "<< displacedMuonL1Pt.getNParity() <<" meRing "<< displacedMuonL1Pt.getMeRing()<<" in ana npar "<< npar <<" ring 2 "<< std::endl;
     }
     if (verbose_){
	std::cout <<"LCT gp1 x "<< gp1.x()<<" y "<< gp1.y() <<" phi "<< gp1.phi()<<" perp "<< gp1.perp() <<" gp2 x "<< gp2.x()<<" y "<< gp2.y()<<" phi "<< gp2.phi()<<" perp "<< gp2.perp() << std::endl;
     	std::cout <<"at L1 deltay12 "<<  etrk_[0].deltay12_lct <<" deltay23 "<<  etrk_[0].deltay23_lct << std::endl;
	std::cout <<"fit gp1 x "<< gp1_fit.x()<<" y "<< gp1_fit.y() <<" phi "<< gp1_fit.phi()<<" perp "<< gp1_fit.perp() <<" gp2 x "<< gp2_fit.x()<<" y "<< gp2_fit.y()<<" phi "<< gp2_fit.phi()<<" perp "<< gp2_fit.perp() << std::endl;
     	std::cout <<"at digi with fit deltay12 "<<  etrk_[0].deltay12_fit <<" deltay23 "<<  etrk_[0].deltay23_fit << std::endl;
     }
     etrk_[0].hasSt1St2St3=true;
     etrk_[0].npar_lct= npar;
     etrk_[0].deltay12_lct = PtassignmentHelper::deltaYcalculation(gp1, gp2);
     etrk_[0].deltay23_lct = -PtassignmentHelper::deltaYcalculation(gp3, gp2);
     etrk_[0].deltay12_fit = PtassignmentHelper::deltaYcalculation(gp1_fit, gp2_fit);
     etrk_[0].deltay23_fit = -PtassignmentHelper::deltaYcalculation(gp3_fit, gp2_fit);
     etrk_[0].pt_position=PtassignmentHelper::Ptassign_Position_gp(gp1, gp2, gp3, gp2.eta(), npar);
     etrk_[0].pt_position_fit=PtassignmentHelper::Ptassign_Position_gp(gp1_fit, gp2_fit, gp3_fit, gp2_fit.eta(), npar);
  }
     if (displacedMuonL1Pt.getNParity()>=0 and displacedMuonL1Pt.runPositionbased()){
    	etrk_[0].deltay12_test = displacedMuonL1Pt.getdeltaY12();
    	etrk_[0].deltay23_test = displacedMuonL1Pt.getdeltaY23();
    	etrk_[0].deltay123_test = displacedMuonL1Pt.getdeltaY123();
	etrk_[0].position_pt = displacedMuonL1Pt.getPositionPt();
	//std::cout <<"DisplacedMuon L1, Position based deltay12 "<< etrk_[0].deltay12_test <<" deltay23 "<< etrk_[0].deltay23_test <<" deltay123 "<< etrk_[0].deltay123_test<<" position pt "<< etrk_[0].position_pt << std::endl;
	displacedMuonL1Pt.runHybrid(true);
	etrk_[0].hybrid_pt = float(displacedMuonL1Pt.getHybridPt());
	displacedMuonL1Pt.runHybrid(false);
	etrk_[0].hybrid_noge21_pt = float(displacedMuonL1Pt.getHybridPt());
	//std::cout <<"simpt "<< etrk_[0].pt <<" eta "<< etrk_[0].eta <<" npar "<< etrk_[0].npar_L1 <<" hybrid pt "<< etrk_[0].hybrid_pt <<" hybrid pt(NOGE21) "<< etrk_[0].hybrid_noge21_pt << std::endl;
     }


     //std::cout <<"DisplacedMuon SIM get npar "<< displacedMuonL1Pt_sim.getNParity()<<" ring "<<displacedMuonL1Pt_sim.getMeRing() << std::endl;
     if (verbose_ and displacedMuonL1Pt_sim.getNParity()>=0 and displacedMuonL1Pt_sim.runPositionbased())
     {
     	std::cout <<"DisplacedMuon SIM deltay12 "<< displacedMuonL1Pt_sim.getdeltaY12() <<" deltay23 "<< displacedMuonL1Pt_sim.getdeltaY23() <<" ddY123 "<< displacedMuonL1Pt_sim.getdeltaY123() <<" etrk deltay12 "<< etrk_[0].deltay12_sim <<" deltay23 "<< etrk_[0].deltay23_sim << std::endl;
     }
     if (verbose_ and displacedMuonL1Pt.getNParity()>=0 and displacedMuonL1Pt.runPositionbased() and etrk_[0].hasSt1St2St3)
     {
	if (fabs(displacedMuonL1Pt.getdeltaY23()-etrk_[0].deltay23_fit)>.1)
	    std::cout <<"sim pt "<< etrk_[0].pt<<" difference between module and ana "<< displacedMuonL1Pt.getdeltaY23()-etrk_[0].deltay23_fit << std::endl;
     	std::cout <<"DisplacedMuon deltay12 "<< displacedMuonL1Pt.getdeltaY12() <<" deltay23 "<< displacedMuonL1Pt.getdeltaY23() <<" ddY123 "<< displacedMuonL1Pt.getdeltaY123() <<" etrk deltay12, sim "<< etrk_[0].deltay12_sim <<" fit "<< etrk_[0].deltay12_fit <<" deltay23,sim "<< etrk_[0].deltay23_sim <<" fit "<< etrk_[0].deltay23_fit << std::endl;

     }else if (verbose_ and etrk_[0].hasSt1St2St3){
    	std::cout <<"Failed to run DisplacedMuonPtassignment, sim pt "<<  etrk_[0].pt <<" GEMCSCAnalyer sim deltay12 "<< etrk_[0].deltay12_sim<<" deltay23 "<< etrk_[0].deltay23_sim << " L1 deltay12 "<< etrk_[0].deltay12_fit <<" deltay23 "<< etrk_[0].deltay23_fit << std::endl;
     }

  if (verbose_) std::cout <<"GEMCSCAnalyzer step9 "<< std::endl;
 //general propagation
  /*auto propagate_odd_gp(match_track.simTrackPropagateGPs_odd());
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
  }*/
  if (match_track.tfTracks().size()) {
    etrk_[0].has_tfTrack = 1;
    TFTrack* besttrack = match_track.bestTFTrack();
    //std::cout <<"tfTracks size "<< match_track.tfTracks().size() <<" bestTrack "; besttrack->print();
    etrk_[0].trackpt = besttrack->pt();
    etrk_[0].tracketa = besttrack->eta();
    etrk_[0].trackphi = besttrack->phi();
  //  quality_packed;
    if (besttrack->tracktype() == 0){
       etrk_[0].pt_packed = besttrack->ptPacked();
       etrk_[0].eta_packed = besttrack->etaPacked();
       etrk_[0].phi_packed = besttrack->phiPacked();
       etrk_[0].quality_packed = besttrack->qPacked();
    }
     // rank = 0;
   etrk_[0].deltaphi12 = besttrack->dPhi12();
   etrk_[0].deltaphi23 = besttrack->dPhi23();
   etrk_[0].hasME1 = besttrack->hasStubEndcap(1);
   etrk_[0].hasME2 = besttrack->hasStubEndcap(2);
   etrk_[0].nstubs = besttrack->nStubs();
   etrk_[0].deltaR = besttrack->dr();
   etrk_[0].chargesign = besttrack->chargesign();
   /*unsigned int lct1 = 999;
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
       //etrk_[0].dphiGE11 = ((besttrack->getTriggerDigis()).at(lct1))->getGEMDPhi();
       etrk_[0].ME1_hs = ((besttrack->getTriggerDigis()).at(lct1))->getStrip();
       etrk_[0].ME1_wg = ((besttrack->getTriggerDigis()).at(lct1))->getKeyWG();
       //etrk_[0].passGE11_simpt = match_lct.passDPhicut(id_me1, etrk_[0].chargesign, etrk_[0].dphiGE11, pt);
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
       //etrk_[0].dphiGE21 = ((besttrack->getTriggerDigis()).at(lct2))->getGEMDPhi();
       etrk_[0].ME2_hs = ((besttrack->getTriggerDigis()).at(lct2))->getStrip();
       etrk_[0].ME2_wg = ((besttrack->getTriggerDigis()).at(lct2))->getKeyWG();
       //etrk_[0].passGE21_simpt = match_lct.passDPhicut(id_me2, etrk_[0].chargesign, etrk_[0].dphiGE21, pt);
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

    }*/
  }

  auto l1GmtCands(match_track.gmts());
  if (l1GmtCands.size() and match_track.bestGMTCand()) {

    etrk_[0].has_gmtCand = 1;
    TFCand* bestGmtCand = match_track.bestGMTCand();
    //std::cout <<"size of l1GMTs "<< l1GmtCands.size() <<" bestGMTCand in GEMCSC "; bestGmtCand->print();
    etrk_[0].bestdRGmtCand = bestGmtCand->dr();
    etrk_[0].L1Mu_pt = bestGmtCand->pt();
    etrk_[0].L1Mu_eta = bestGmtCand->eta();
    etrk_[0].L1Mu_phi = bestGmtCand->phi();
    etrk_[0].L1Mu_charge = bestGmtCand->charge();
    etrk_[0].L1Mu_bx = bestGmtCand->bx();
    etrk_[0].L1Mu_quality = bestGmtCand->quality();

    float tfeta = bestGmtCand->eta();
    float tfphi = bestGmtCand->phi();
    float mindPhi = .40;

    for (auto me0Segment : allmatchedSegments){
      GlobalPoint gp_ME0_st2(match_me0rh.propagateFromME0ToCSC(me0Segment, etrk_[0].L1Mu_pt, 2, me0Segment.me0DetId().chamber()%2==1)); 
      GlobalPoint gp(match_me0rh.globalPoint(me0Segment));
      float dR = deltaR(tfeta, tfphi, float(gp.eta()), float(gp_ME0_st2.phi()));
      float dPhi_L1Mu_ME0_propagate = deltaPhi(tfphi, float(gp_ME0_st2.phi()));
      float dPhi_L1Mu_ME0 = deltaPhi(tfphi, float(gp.phi()));
      //std::cout <<"L1Mu eta "<< tfeta <<" phi "<< tfphi <<" propagated ME0 eta "<< gp_ME0_st2.eta() <<" phi "<< gp_ME0_st2.phi()<<" ME0 eta "<< gp.eta() <<" phi "<< gp.phi() <<" dPhi "<< dPhi_L1Mu_ME0 << std::endl;


      if (fabs(dPhi_L1Mu_ME0_propagate) < mindPhi and fabs(gp.eta() - tfeta)<0.25){
	  mindPhi = std::fabs(dPhi_L1Mu_ME0);
	  int chamber = me0Segment.me0DetId().chamber();
	  if (fabs(etrk_[1].phi_lct_even) < 4 and abs(etrk_[1].chamber_lct_even/2-chamber) <= 1){
	      etrk_[0].L1Mu_me0_st1_dphi = deltaPhi(etrk_[1].phi_lct_even, gp.phi());
	      etrk_[0].L1Mu_me0_st1_isEven = true;
	  }else if (fabs(etrk_[1].phi_lct_odd) < 4 and abs(etrk_[1].chamber_lct_odd/2-chamber) <= 1){
	      etrk_[0].L1Mu_me0_st1_dphi = deltaPhi(etrk_[1].phi_lct_odd, gp.phi());
	      etrk_[0].L1Mu_me0_st1_isEven = false;
	  }else etrk_[0].L1Mu_me0_st1_dphi = -9;

	  float dPhi = match_me0rh.me0DeltaPhi(me0Segment);
	  etrk_[0].L1Mu_me0_st2_eta = gp_ME0_st2.eta();
	  etrk_[0].L1Mu_me0_st2_phi = gp_ME0_st2.phi();
	  etrk_[0].L1Mu_me0_eta = gp.eta();
	  etrk_[0].L1Mu_me0_phi = gp.phi();
	  etrk_[0].L1Mu_me0_dR = dR;
	  etrk_[0].L1Mu_me0_mindPhi1 = dPhi_L1Mu_ME0;
	  etrk_[0].L1Mu_me0_mindPhi2 = dPhi_L1Mu_ME0_propagate;
	  etrk_[0].L1Mu_me0_dPhi = dPhi;
      }
    }
  }

  /*if (match_track.tfCands().size()) {
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
  }*/

  if (verbose_) std::cout <<"GEMCSCAnalyzer step10 "<< std::endl;
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
        dtrk_.dphi_sh = deltaPhi(float(csc_sh_gp.phi()), float(gem_sh_gp.phi()));
        dtrk_.dphi_dg = deltaPhi(float(csc_dg_gp.phi()), float(gem_dg_gp.phi()));
        dtrk_.dphi_pad = deltaPhi(float(csc_dg_gp.phi()), float(gem_pad_gp.phi()));
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
          dtrk_.dphi_lct_pad = deltaPhi(float(csc_lct_gp.phi()), float(gem_pad_gp.phi()));
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
  //const L1TrackMatcher& match_track = match.l1Tracks();
  const UpgradeL1TrackMatcher& match_track = match.l1Tracks();
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
    int MEStation = id.station();

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
    int MEStation = id.station();

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
    int MEStation = id.station();

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
     //std::cout << " propagated information " << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME1 << " phi " << etrk_[0].phi_propagated_ME1 << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME2 << " phi " << etrk_[0].phi_propagated_ME2 << std::endl;
    // std::cout << " eta " << etrk_[0].eta_propagated_ME3 << " phi " << etrk_[0].phi_propagated_ME3 << std::endl;
     //std::cout << " propagated phi in  ME1 " << etrk_[0].phi_propagated_ME1 <<" stub phi in ME1 " <<etrk_[0].phi_ME1_TF << std::endl;
     //std::cout << " propagated phi in  ME2 " << etrk_[0].phi_interStat12 <<" stub phi in ME2 " << etrk_[0].phi_ME2_TF << std::endl;
     ///std::cout << " propagated phi in  ME3 " << etrk_[0].phi_interStat23 <<" stub phi in ME3 " << etrk_[0].phi_ME3_TF << std::endl;


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

