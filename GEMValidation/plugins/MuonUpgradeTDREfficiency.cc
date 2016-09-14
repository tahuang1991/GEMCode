#ifndef GEMCode_GEMValidation_MuonUpgradeTDREfficiency
#define GEMCode_GEMValidation_MuonUpgradeTDREfficiency

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

#include "TTree.h"

#include <iomanip>
#include <sstream>
#include <memory>

using namespace std;
using namespace matching;


// "signed" LCT bend pattern
const int LCT_BEND_PATTERN[11] = { -99,  -5,  4, -4,  3, -3,  2, -2,  1, -1,  0};


struct TrackEff
{
  void init(); // initialize to default values
  TTree* book(TTree *t, const std::string & name = "trk_eff");

  Int_t lumi;
  Int_t run;
  Int_t event;

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
  UChar_t wg_lct_odd;
  UChar_t hs_lct_even;
  UChar_t wg_lct_even;

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

};

void TrackEff::init()
{
  lumi = -99;
  run = -99;
  event = -99;

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
  wg_lct_odd = 0;
  wg_lct_even = 0;
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
}


TTree* TrackEff::book(TTree *t, const std::string & name)
{
  edm::Service< TFileService > fs;
  t = fs->make<TTree>(name.c_str(), name.c_str());

  t->Branch("lumi", &lumi);
  t->Branch("run", &run);
  t->Branch("event", &event);

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
  t->Branch("wg_lct_even", &wg_lct_even);
  t->Branch("wg_lct_odd", &wg_lct_odd);
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

  return t;
}

// --------------------------- MuonUpgradeTDREfficiency ---------------------------

class MuonUpgradeTDREfficiency : public edm::EDAnalyzer
{
public:

  explicit MuonUpgradeTDREfficiency(const edm::ParameterSet&);

  ~MuonUpgradeTDREfficiency() {}
  
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  void analyzeTrackEff(SimTrackMatchManager& match, int trk_no);
  void printout(SimTrackMatchManager& match, int trk_no);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);
  
  edm::ParameterSet cfg_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertexInput_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackInput_;

  int verboseSimTrack_;
  double simTrackMinPt_;
  double simTrackMinEta_;
  double simTrackMaxEta_;
  double simTrackOnlyMuon_;
  int verbose_;
  bool ntupleTrackEff_;
  bool matchprint_;
  std::vector<string> cscStations_;
  std::vector<std::pair<int,int> > cscStationsCo_;
  std::set<int> stations_to_use_;

  TTree *tree_eff_[12]; // for up to 9 stations
  
  TrackEff  etrk_[12];

  int minNHitsChamberCSCSimHit_;
  int minNHitsChamberCSCWireDigi_;
  int minNHitsChamberCSCStripDigi_;
  int minNHitsChamberCLCT_;
  int minNHitsChamberALCT_;
  int minNHitsChamberLCT_;
  int minNHitsChamberMPLCT_;
};


MuonUpgradeTDREfficiency::MuonUpgradeTDREfficiency(const edm::ParameterSet& ps)
: cfg_(ps.getParameterSet("simTrackMatching"))
, verbose_(ps.getUntrackedParameter<int>("verbose", 0))
{
  cscStations_ = cfg_.getParameter<std::vector<string> >("cscStations");
  ntupleTrackEff_ = cfg_.getParameter<bool>("ntupleTrackEff");
  matchprint_ = false; //cfg_.getParameter<bool>("matchprint");

  auto simVertex = cfg_.getParameter<edm::ParameterSet>("simVertex");
  simVertexInput_ = consumes<edm::SimVertexContainer>(simVertex.getParameter<edm::InputTag>("validInputTags"));

  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrack.getParameter<int>("verbose");
  simTrackInput_ = consumes<edm::SimTrackContainer>(simTrack.getParameter<edm::InputTag>("validInputTags"));
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


int MuonUpgradeTDREfficiency::detIdToMEStation(int st, int ri)
{
  auto p(std::make_pair(st, ri));
  return std::find(cscStationsCo_.begin(), cscStationsCo_.end(), p) - cscStationsCo_.begin();
}


void MuonUpgradeTDREfficiency::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
}


bool MuonUpgradeTDREfficiency::isSimTrackGood(const SimTrack &t)
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


void MuonUpgradeTDREfficiency::analyze(const edm::Event& ev, const edm::EventSetup& es)
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

    SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es, consumesCollector());

    if (ntupleTrackEff_) analyzeTrackEff(match, trk_no);
    ++trk_no;
  }
}



void MuonUpgradeTDREfficiency::analyzeTrackEff(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  const SimTrack &t = match_sh.trk();
   
  for (auto s: stations_to_use_)
  {

    etrk_[s].init();
    etrk_[s].run = match.simhits().event().id().run();
    etrk_[s].lumi = match.simhits().event().id().luminosityBlock();
    etrk_[s].event = match.simhits().event().id().event();
    etrk_[s].pt = t.momentum().pt();
    etrk_[s].phi = t.momentum().phi();
    etrk_[s].eta = t.momentum().eta();
    etrk_[s].charge = t.charge();
    etrk_[s].endcap = (etrk_[s].eta > 0.) ? 1 : -1;
  }

  // SimHits
  auto csc_simhits(match_sh.chamberIdsCSC(0));
  for(auto d: csc_simhits)
  {

    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    int nlayers(match_sh.nLayersWithHitsInSuperChamber(d));

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
      etrk_[st].wg_lct_odd = digi_wg(lct);
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
      etrk_[st].wg_lct_even = digi_wg(lct);
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
        etrk_[1].wg_lct_odd = digi_wg(lct);
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
        etrk_[1].wg_lct_even = digi_wg(lct);
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

    const int st(detIdToMEStation(MEStation,id.ring()));
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

    auto pads = match_gd.padsInSuperChamber(d);
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
        etrk_[st].dphi_pad_odd = reco::deltaPhi((float)etrk_[st].phi_lct_odd, (float)etrk_[st].phi_pad_odd);
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
        etrk_[st].dphi_pad_even = reco::deltaPhi((float)etrk_[st].phi_lct_even, (float)etrk_[st].phi_pad_even);
        etrk_[st].deta_pad_even = etrk_[st].eta_lct_even - etrk_[st].eta_pad_even;
      }
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
        etrk_[st].dphi_pad_odd = reco::deltaPhi((float)etrk_[st].phi_lct_odd, (float)etrk_[st].phi_pad_odd);
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
        etrk_[st].dphi_pad_even = reco::deltaPhi((float)etrk_[st].phi_lct_even, (float)etrk_[st].phi_pad_even);
        etrk_[st].deta_pad_even = etrk_[st].eta_lct_even - etrk_[st].eta_pad_even;
      }
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

  for (auto s: stations_to_use_)
  {
    tree_eff_[s]->Fill();
  }
}


 void MuonUpgradeTDREfficiency::printout(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  const GEMDigiMatcher& match_gd = match.gemDigis();
  const CSCDigiMatcher& match_cd = match.cscDigis();
  const CSCStubMatcher& match_lct = match.cscStubs();
  //  const L1TrackMatcher& match_track = match.l1Tracks();
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
       std::cout<<id<< p <<std::endl;
    
  }

  std::cout << "######matching ALCT to Simtrack " << std::endl;
  for(auto d: match_lct.chamberIdsAllALCT(0))
  {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stations_to_use_.count(st) == 0) continue;
    auto alcts = match_lct.allALCTsInChamber(d);
    for (auto p : alcts)    
       std::cout<<id<< p <<std::endl;
    
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
    
  }


  // std::cout << "######  matching Tracks to Simtrack " << std::endl;
  // if (match_track.tfTracks().size()) {
  //   TFTrack* besttrack = match_track.bestTFTrack();
  //   std::cout << "       Best TFTrack                  " << std::endl;
  //   besttrack->print();


  // }
  // else std::cout << "NO matched TFtracks"  << std::endl;


  std::cout << "==========================  end of printing ========================\n\n" << std::endl;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonUpgradeTDREfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MuonUpgradeTDREfficiency);

#endif
