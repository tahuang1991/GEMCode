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

#include "GEMCode/GEMValidation/src/SimTrackMatchManager.h"

#include "TTree.h"

#include <iomanip>
#include <sstream>
#include <memory>

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
};


struct MyTrackEff
{
  void init(); // initialize to default values
  TTree* book(TTree *t, const std::string & name = "trk_eff");

  Float_t pt, eta, phi;
  Char_t charge;
  Char_t endcap;
  Char_t chamber_odd; // bit1: has GEM pad   bit2: has CSC LCT
  Char_t chamber_even; // bit1: has GEM pad   bit2: has CSC LCT

  Char_t has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
  Char_t has_csc_strips; // #layers with comparator digis > minHitsChamber    bit1: in odd, bit2: even
  Char_t has_csc_wires; // #layers with wire digis > minHitsChamber    bit1: in odd, bit2: even

  Char_t has_clct; // bit1: in odd, bit2: even
  Char_t has_alct; // bit1: in odd, bit2: even
  Char_t has_lct; // bit1: in odd, bit2: even

  Char_t bend_lct_odd;
  Char_t bend_lct_even;
  Char_t bx_lct_odd;
  Char_t bx_lct_even;
  UChar_t hs_lct_odd;
  UChar_t hs_lct_even;

  Float_t phi_lct_odd;
  Float_t phi_lct_even;
  Float_t eta_lct_odd;
  Float_t eta_lct_even;
  Float_t dphi_lct_odd; // dphi stored as data member in LCT
  Float_t dphi_lct_even;

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

};

struct AllLCT
{
    
  void init(); // initialize to default values
  TTree* bookLCTTree(TTree *t);

  Float_t phi;
  Float_t eta;
  Float_t pt;
  Float_t charge;

  Bool_t match;
  Bool_t passdphicut;

  Int_t station;
  Int_t chamber;
  Int_t ring;
  Int_t endcap;

  Float_t dphi;
  Int_t hs;
  Int_t wg;
  Int_t bx;
  Int_t pattern;
  Int_t quality;
  
};


void MyTrackEff::init()
{
  pt = 0.;
  phi = 0.;
  eta = -9.;
  charge = -9;
  endcap = -9;
  chamber_odd = 0;
  chamber_even = 0;
  quality_odd = 0;
  quality_even = 0;

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
  phi_lct_odd = -9.;
  phi_lct_even = -9.;
  eta_lct_odd = -9.;
  eta_lct_even = -9.;
  dphi_lct_odd = -9.;
  dphi_lct_even = -9.;

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
}

    
void AllLCT::init(){ // initialize to default values

  phi=0.0;
  eta=0.0;
  pt=0.0;
  charge=0.0;
 
  match=false;
  passdphicut=false;

  station=-1;
  chamber=-1;
  ring=-1;
  endcap=-2;

  dphi=-99.0;
  hs=-1;
  wg=-1;
  bx=-1;
  pattern=-1;
  quality=-1;

  }

TTree* MyTrackEff::book(TTree *t, const std::string & name)
{
  edm::Service< TFileService > fs;
  t = fs->make<TTree>(name.c_str(), name.c_str());

  t->Branch("pt", &pt);
  t->Branch("eta", &eta);
  t->Branch("phi", &phi);
  t->Branch("charge", &charge);
  t->Branch("endcap", &endcap);
  t->Branch("chamber_odd", &chamber_odd);
  t->Branch("chamber_even", &chamber_even);
  t->Branch("quality_odd", &quality_odd);
  t->Branch("quality_even", &quality_even);
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
  t->Branch("phi_lct_odd", &phi_lct_odd);
  t->Branch("phi_lct_even", &phi_lct_even);
  t->Branch("eta_lct_odd", &eta_lct_odd);
  t->Branch("eta_lct_even", &eta_lct_even);
  t->Branch("dphi_lct_odd", &dphi_lct_odd);
  t->Branch("dphi_lct_even", &dphi_lct_even);

  
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
  
  return t;
}

TTree* AllLCT::bookLCTTree(TTree *t)
{

  edm::Service< TFileService > fs;
  t = fs->make<TTree>("LCTTree", "LCTTree");

  t->Branch("pt", &pt);
  t->Branch("eta", &eta);
  t->Branch("phi", &phi);
  t->Branch("match", &match);
  t->Branch("passdphicut", &passdphicut);
  t->Branch("charge", &charge);
  t->Branch("endcap", &endcap);
  t->Branch("station", &station);
  t->Branch("chamber", &chamber);
  t->Branch("ring", &ring);
  t->Branch("endcap", &endcap);
  t->Branch("dphi", &dphi);
  t->Branch("hs", &hs);
  t->Branch("wg", &wg);
  t->Branch("bx", &bx);
  t->Branch("pattern", &pattern);
  t->Branch("quality", &quality);

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
  

  const double ME11GEMdPhi[9][3] = {
    {-2 , 1.0, 1.0 },
    {3 , 0.03971647, 0.01710244 },
    {5 , 0.02123785, 0.00928431 },
    {7 , 0.01475524, 0.00650928 },
    {10, 0.01023299, 0.00458796 },
    {15, 0.00689220, 0.00331313 },
    {20, 0.00535176, 0.00276152 },
    {30, 0.00389050, 0.00224959 },
    {40, 0.00329539, 0.00204670 }
  };
  const double ME21GEMdPhi[9][3] = {
    {-2 , 1.0, 1.0 },
    {3 , 0.01832829, 0.01003643 },
    {5 , 0.01095490, 0.00631625 },
    {7 , 0.00786026, 0.00501017 },
    {10, 0.00596349, 0.00414560 },
    {15, 0.00462411, 0.00365550 },
    {20, 0.00435298, 0.00361550 },
    {30, 0.00465160, 0.00335700 },
    {40, 0.00372145, 0.00366262 }
  };

  void bookSimTracksDeltaTree();

  void analyzeTrackChamberDeltas(SimTrackMatchManager& match, int trk_no);
  void analyzeTrackEff(SimTrackMatchManager& match, int trk_no);
  void printout(SimTrackMatchManager& match, int trk_no);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);
  
  edm::ParameterSet cfg_;
  edm::InputTag simInputLabel_;
  double simTrackMinPt_;
  double simTrackMinEta_;
  double simTrackMaxEta_;
  double simTrackOnlyMuon_;
  int verbose_;
  bool ntupleTrackChamberDelta_;
  bool ntupleTrackEff_;
  bool matchprint_;
  std::vector<string> cscStations_;
  std::vector<std::pair<int,int> > cscStationsCo_;
  std::set<int> stations_to_use_;

  TTree *tree_eff_[12]; // for up to 9 stations
  TTree *tree_delta_;
  TTree *tree_lct;
    
  MyTrackEff  etrk_[12];
  MyTrackChamberDelta dtrk_;

  AllLCT LCTStub;

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
  matchprint_ = cfg_.getParameter<bool>("matchprint");

  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  simInputLabel_ = simTrack.getParameter<edm::InputTag>("input");
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

  /*
  auto tfTrack = cfg_.getParameter<edm::ParameterSet>("tfTrack");
  auto tfCand = cfg_.getParameter<edm::ParameterSet>("tfCand");
  auto gmtCand = cfg_.getParameter<edm::ParameterSet>("gmtCand");
  auto l1Extra = cfg_.getParameter<edm::ParameterSet>("l1Extra");
  */
  if (ntupleTrackChamberDelta_) bookSimTracksDeltaTree();
  tree_lct = LCTStub.bookLCTTree(tree_lct);
  if (ntupleTrackEff_)
  {
    vector<int> stations = ps.getParameter<vector<int> >("stationsToUse");
    copy(stations.begin(), stations.end(), inserter(stations_to_use_, stations_to_use_.end()) );

    for(auto s: stations_to_use_)
    {
      stringstream ss;
      ss << "trk_eff_"<< cscStations_[s];
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
  //
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
  edm::Handle<edm::SimVertexContainer> sim_vertices;

  ev.getByLabel(simInputLabel_, sim_tracks);
  ev.getByLabel(simInputLabel_, sim_vertices);
  const edm::SimVertexContainer & sim_vert = *sim_vertices.product();

  /*
  {  // print out 1st strip coordinates for rolls in GE1/1 chamber
    edm::ESHandle<GEMGeometry> gem_g;
    es.get<MuonGeometryRecord>().get(gem_g);
    const GEMGeometry * gem_geo = &*gem_g;
    for (int r=1; r<7; ++r)
    {
      GEMDetId p(1, 1, 1, 1, 1, r);
      auto roll = gem_geo->etaPartition(p);
      auto lp = roll->centreOfStrip(1);
      GlobalPoint gp = gem_geo->idToDet(p())->surface().toGlobal(lp);
      cout<<setprecision(9)<<"rollp "<<r<<" "<<gp.phi()<<" "<<gp.perp()<<" "<<roll->localPitch(lp)<<" "<<roll->localPitch(lp)/gp.perp()<<endl;
    }
  }
  */

  /*
  // print out 1st strip coordinates for ME1/b chamber
  edm::ESHandle<CSCGeometry> csc_g;
  es.get<MuonGeometryRecord>().get(csc_g);
  const CSCGeometry * csc_geo = &*csc_g;
  for (int nmb=0;nmb<4;++nmb) for (int la=1; la<7; ++la){
    CSCDetId id(1, 1, 1, 1, la);
    if (nmb==1) id = CSCDetId(1,1,4,1,la);
    if (nmb==2) id = CSCDetId(1,1,2,1,la);
    if (nmb==3) id = CSCDetId(1,2,1,1,la);
    auto strip_topo = csc_geo->layer(id)->geometry()->topology();
    MeasurementPoint mp_top(0.25, 0.5);
    MeasurementPoint mp_bot(0.25, -0.5);
    LocalPoint lp = strip_topo->localPosition(0.25);
    LocalPoint lp_top = strip_topo->localPosition(mp_top);
    LocalPoint lp_bot = strip_topo->localPosition(mp_bot);
    GlobalPoint gp = csc_geo->idToDet(id)->surface().toGlobal(lp);
    GlobalPoint gp_top = csc_geo->idToDet(id)->surface().toGlobal(lp_top);
    GlobalPoint gp_bot = csc_geo->idToDet(id)->surface().toGlobal(lp_bot);
    cout<<id<<endl;
    cout<<setprecision(6)<<"glayer "<<la<<" "<<gp.phi()<<" "<<gp.perp()<<" "<<strip_topo->localPitch(lp)<<" "<<strip_topo->localPitch(lp)/gp.perp()
        <<"  "<<gp_top.phi()<<" "<<gp_top.perp()<<" "<<strip_topo->localPitch(lp_top)<<" "<<strip_topo->localPitch(lp_top)/gp_top.perp()
        <<"  "<<gp_bot.phi()<<" "<<gp_bot.perp()<<" "<<strip_topo->localPitch(lp_bot)<<" "<<strip_topo->localPitch(lp_bot)/gp_bot.perp()
        <<endl;
  }
  */
  int trk_no=0;
  for (auto& t: *sim_tracks.product())
  {
    if (!isSimTrackGood(t)) continue;

    // match hits and digis to this SimTrack
    SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);

    if (ntupleTrackChamberDelta_) analyzeTrackChamberDeltas(match, trk_no);
    if (ntupleTrackEff_) analyzeTrackEff(match, trk_no);
    // if (matchprint_) printout(match, trk_no);
    
        
          int st = detIdToMEStation(2, 1);
          //bool has_alct(etrk_[st].has_alct>0); bool has_clct(etrk_[st].has_clct>0) ;
	  bool has_stubs(etrk_[st].has_alct>0 && etrk_[st].has_gem_copad>0);
	  bool no_lct(etrk_[st].has_lct==0);
          // if (has_csc_sh_odd || has_csc_sh_even)  std::cout <<"st1 has_csc_sh " << std::endl;
          // if (has_alct_odd || has_alct_even)   std::cout <<"  st1 has_alct " << std::endl;
          bool Debug(has_stubs && no_lct);
          if (matchprint_ and Debug ) printout(match, trk_no);
          trk_no++;
   
  }
}



void GEMCSCAnalyzer::analyzeTrackEff(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const RPCDigiMatcher& match_rd = match.rpcDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  //const TrackMatcher& match_track = match.tracks();
  const SimTrack &t = match_sh.trk();
   
  for (auto s: stations_to_use_)
  {
    etrk_[s].init();
    etrk_[s].pt = t.momentum().pt();
    etrk_[s].phi = t.momentum().phi();
    etrk_[s].eta = t.momentum().eta();
    etrk_[s].charge = t.charge();
    etrk_[s].endcap = (etrk_[s].eta > 0.) ? 1 : -1;
  }

  // SimHits
  for(auto d: match_sh.chamberIdsCSC(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    const int nlayers(match_sh.nLayersWithHitsInSuperChamber(d));
    if (nlayers < minNHitsChamberCSCSimHit_) continue;

    const bool odd(id.chamber()%2==1);
    if (odd) etrk_[st].has_csc_sh |= 1;
    else etrk_[st].has_csc_sh |= 2;

    if (odd) etrk_[st].nlayers_csc_sh_odd = nlayers;
    else etrk_[st].nlayers_csc_sh_even = nlayers;
    
    // std::cout<<"nlayer "<<nlayers <<" odd: "<< etrk_[st].nlayers_csc_sh_odd<<" even: "<<etrk_[st].nlayers_csc_sh_even
//	 <<" "<< (odd ? "odd":"even")<<" csc det " <<id <<std::endl;
  //  std::cout<<" "<<((etrk_[st].has_csc_sh&1)>0 ? "odd true":"odd false" ) <<" "<<((etrk_[st].has_csc_sh&2)>0 ? "even true":"even false")<<std::endl;
    // case ME11
    if (st==2 or st==3){
      if (odd) etrk_[1].has_csc_sh |= 1;
      else etrk_[1].has_csc_sh |= 2;

      if (odd) etrk_[1].nlayers_csc_sh_odd = nlayers;
      else etrk_[1].nlayers_csc_sh_even = nlayers;
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
      etrk_[st].dphi_lct_odd = digi_dphi(lct);
      etrk_[st].bx_lct_odd = digi_bx(lct);
      etrk_[st].hs_lct_odd = digi_channel(lct);
      etrk_[st].chamber_odd |= 2;
      etrk_[st].quality_odd = digi_quality(lct);
    }
    else
    {
      lct_even[st] = lct;
      gp_lct_even[st] = gp;
      etrk_[st].bend_lct_even = bend;
      etrk_[st].phi_lct_even = gp.phi();
      etrk_[st].eta_lct_even = gp.eta();
      etrk_[st].dphi_lct_even = digi_dphi(lct);
      etrk_[st].bx_lct_even = digi_bx(lct);
      etrk_[st].hs_lct_even = digi_channel(lct);
      etrk_[st].chamber_even |= 2;
      etrk_[st].quality_even = digi_quality(lct);
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
        etrk_[1].dphi_lct_odd = digi_dphi(lct);
        etrk_[1].bx_lct_odd = digi_bx(lct);
        etrk_[1].hs_lct_odd = digi_channel(lct);
        etrk_[1].chamber_odd |= 2;
        etrk_[1].quality_odd = digi_quality(lct);
      }
      else
      {
        lct_even[1] = lct;
        gp_lct_even[1] = gp;
        etrk_[1].bend_lct_even = bend;
        etrk_[1].phi_lct_even = gp.phi();
        etrk_[1].eta_lct_even = gp.eta();
        etrk_[1].dphi_lct_even = digi_dphi(lct);
        etrk_[1].bx_lct_even = digi_bx(lct);
        etrk_[1].hs_lct_even = digi_channel(lct);
        etrk_[1].chamber_even |= 2;
        etrk_[1].quality_even = digi_quality(lct);
      }
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
    
    int ring = id.ring();
    if (std::fabs(t.momentum().eta())>2.1 && MEStation==1) ring = 4;
    const int st(detIdToMEStation(MEStation,ring));
    if (stations_to_use_.count(st) == 0) continue;

    const bool odd(id.chamber()%2==1);
    if (match_sh.hitsInSuperChamber(d).size() > 0)
    {
      if (odd) etrk_[st].has_gem_sh |= 1;
      else     etrk_[st].has_gem_sh |= 2;

      auto sh_gp = match_sh.simHitsMeanPosition(match_sh.hitsInSuperChamber(d));
      if (odd) etrk_[st].eta_gemsh_odd = sh_gp.eta();
      else     etrk_[st].eta_gemsh_even = sh_gp.eta();

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

      auto sh_gp = match_sh.simHitsMeanPosition(match_sh.hitsInSuperChamber(d));
      if (odd) etrk_[1].eta_gemsh_odd = sh_gp.eta();
      else     etrk_[1].eta_gemsh_even = sh_gp.eta();

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
  for(auto d: match_gd.superChamberIds())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
   
    int ring = id.ring(); 
    if (std::fabs(t.momentum().eta())>2.1 && MEStation==1) ring = 4;
    const int st(detIdToMEStation(MEStation,ring));
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

    auto pads = match_gd.padsInSuperChamber(d);
    if(pads.size() == 0) continue;
    if (odd)
    {
      etrk_[st].has_gem_pad |= 1;
      etrk_[st].chamber_odd |= 1;
      etrk_[st].pad_odd = digi_channel(pads.at(0));
      etrk_[st].hsfromgem_odd = match_gd.extrapolateHsfromGEMPad( d, ring, digi_channel(pads.at(0)));
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
      etrk_[st].hsfromgem_even = match_gd.extrapolateHsfromGEMPad( d, ring, digi_channel(pads.at(0)));
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
  }

//ME11Case
  for(auto d: match_gd.superChamberIds())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
    
    int ring = id.ring();
    if (std::fabs(t.momentum().eta())>2.1 && MEStation==1) ring = 4;
    const int stations(detIdToMEStation(MEStation,ring));
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

    auto pads = match_gd.padsInSuperChamber(d);
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
   }

  for(auto d: match_gd.superChamberIdsWithCoPads())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
    
    int ring = id.ring();
    if (std::fabs(t.momentum().eta())>2.1 && MEStation==1) ring = 4;
    const int st(detIdToMEStation(MEStation,ring));
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

  auto rpc_ch_ids = match_sh.chamberIdsRPC();
  for (auto d:rpc_ch_ids)
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


  rpc_ch_ids = match_rd.detIds(); 
  //rpc_ch_ids = match_rd.chamberIds(); 
  for (auto d:rpc_ch_ids)
  {
    RPCDetId id(d);
    const int st(detIdToMEStation(id.station(), id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    //meanstrip in rpc 
    auto rpcdigis = match_rd.digisInDetId(id); 
    int rpc_medianstrip(match_rd.median(rpcdigis));
    int cscchamber = CSCTriggerNumbering::chamberFromTriggerLabels(id.sector(), 0, id.station(), id.subsector());
    cscchamber = (cscchamber+16)%18+1;
    //std::cout <<"rpc detid " << id << " csc chamebr:"<< cscchamber << std::endl;
    bool odd(cscchamber%2 == 1);
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

//fill tree_lct
  for(auto d: match_lct.chamberIdsAllLCT(0))
  {
    CSCDetId id(d);
    
    //if (t.momentum().eta()*id.endcap()<0) continue;
    //std::cout<<"eta "<<t.momentum().eta()<<" endcap "<<id.endcap()<<std::endl;
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;

    auto lcts = match_lct.allLCTsInChamber(d);
    auto lct = match_lct.lctInChamber(d);
    for (auto p : lcts)   
    {
      LCTStub.init();
      if (p==lct)  LCTStub.match = true;
      else LCTStub.match = false;
      LCTStub.eta = t.momentum().eta();
      LCTStub.phi = t.momentum().phi();
      LCTStub.pt = t.momentum().pt();
      LCTStub.charge = t.charge();

      LCTStub.endcap = id.endcap();
      LCTStub.station = id.station();
      LCTStub.chamber = id.chamber();
      LCTStub.ring = id.ring();

      LCTStub.dphi = digi_dphi(p);
      LCTStub.bx = digi_bx(p);
      LCTStub.hs = digi_channel(p);
      LCTStub.wg = digi_wg(p);
      LCTStub.quality = digi_quality(p);
      bool is_odd = (LCTStub.chamber%2==1);
     if (LCTStub.station == 1)
     {
      for (int b = 0; b < 9; b++)
      { // cutting on gem csc dPhi
	 if (double(LCTStub.pt) >= ME11GEMdPhi[b][0])
	 {
             if ((is_odd && ME11GEMdPhi[b][1] > fabs(LCTStub.dphi)) || 
		    (!is_odd && ME11GEMdPhi[b][2] > fabs(LCTStub.dphi)))
		  LCTStub.passdphicut = true;
	     else LCTStub.passdphicut = false;
	 }
      }
      if (LCTStub.dphi<-50)  LCTStub.passdphicut = true;
      if (!LCTStub.passdphicut && std::abs(LCTStub.eta)>2.05 && LCTStub.match>0) printout(match, trk_no); 
     }

     if (LCTStub.station == 2)
     {
      for (int b = 0; b < 9; b++)
      { // cutting on gem csc dPhi
	 if (double(LCTStub.pt) >= ME21GEMdPhi[b][0])
	 {
             if ((is_odd && ME21GEMdPhi[b][1] > fabs(LCTStub.dphi)) || 
		    (!is_odd && ME21GEMdPhi[b][2] > fabs(LCTStub.dphi)))
		  LCTStub.passdphicut = true;
	     else LCTStub.passdphicut = false;
	 }
      }
      if (LCTStub.dphi<-50)  LCTStub.passdphicut = true;
     }
     
      tree_lct->Fill(); 
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
    cout<<"n_gd_ids "<<match_gd.detIds().size()<<endl;
    cout<<"n_gd_ids_copad "<<match_gd.detIdsWithCoPads().size()<<endl;
    auto gem_gd_sch_ids = match_gd.superChamberIds();
    cout<<"n_gd_ids_sch "<<gem_gd_sch_ids.size()<<endl;
    cout<<"n_gd_ids_cosch "<<match_gd.superChamberIdsWithCoPads().size()<<endl;
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
    auto gdids = match_gd.detIds();
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
    auto gem_d_sch_ids = match_gd.superChamberIds();
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
}


 void GEMCSCAnalyzer::printout(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const RPCDigiMatcher& match_rd = match.rpcDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  //  const TrackMatcher& match_track = match.tracks();
  const SimTrack &t = match_sh.trk();

  
  std::cout << "======================== matching information ========================= " << std::endl;
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
  
  for (auto d: match_sh.chamberIdsRPC())
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
  for(auto d: match_gd.superChamberIds())
  {
    GEMDetId id(d);
    int MEStation;
    if (id.station() == 3) MEStation = 2;
    else if (id.station() == 2) continue;
    else MEStation = id.station();
    
    int ring = id.ring();
    if (std::fabs(t.momentum().eta())>2.1 && MEStation==1) ring = 4;
    const int st(detIdToMEStation(MEStation,ring));
    if (stations_to_use_.count(st) == 0) continue;
    
    int nlayers = match_gd.nLayersWithDigisInSuperChamber(d);
    auto digis = match_gd.digisInSuperChamber(d);
    int median_strip = match_gd.median(digis);
    int hs = match_gd.extrapolateHsfromGEMStrip( d, ring, median_strip);
    std::cout <<"GEM Chamber: "<<d<<" "<<id<<" layerswithhits:"<<nlayers
              <<" Medianstrip in Digi:" <<median_strip<<" hs:" << hs<<std::endl;
    // std::cout <<"GEM Pads:"  ;
    auto pads = match_gd.padsInSuperChamber(d);
    for ( auto p=pads.begin(); p != pads.end(); p++)
      std::cout << "  "<< *p <<std::endl; 
  }
  
  std::cout << "######matching Copad to simtrack " << std::endl;
  for (auto d: match_gd.superChamberIdsWithCoPads())
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
      {  std::cout <<"  "<< *p ; }
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

  std::cout << "######matching CLCT to Simtrack " << std::endl;
  for(auto d: match_lct.chamberIdsAllCLCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    auto clcts = match_lct.allCLCTsInChamber(d);
//    auto clct = match_lct.clctInChamber(d);
//    if (std::find(clcts.begin(),clcts.end(),clct) != clcts.end())  std::cout<<"the matching clct ";
//    else std::cout <<" another clct "; 
    for (auto p : clcts)    
       std::cout<<id<<" "<< p <<std::endl;
    
  }

  std::cout << "######matching ALCT to Simtrack " << std::endl;
  for(auto d: match_lct.chamberIdsAllALCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    auto alcts = match_lct.allALCTsInChamber(d);
    for (auto p : alcts)    
       std::cout<<id<<" "<< p <<std::endl;
    
  }

  std::cout << "######matching LCT to Simtrack " << std::endl;
  for(auto d: match_lct.chamberIdsAllLCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    auto lcts = match_lct.allLCTsInChamber(d);
    for (auto p : lcts)    
       std::cout<<id<<" " <<p <<std::endl;
    
  }


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

