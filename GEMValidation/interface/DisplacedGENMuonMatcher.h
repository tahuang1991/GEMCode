#ifndef GEMCode_GEMValidation_DisplacedGENMuonMatcher_h
#define GEMCode_GEMValidation_DisplacedGENMuonMatcher_h

/**\class DisplacedGENMuonMatcher

 Description: Matching of displaced GEN muon to SimTrack 

 Original Author:  "Sven Dildick"
*/

#include "GEMCode/GEMValidation/interface/BaseMatcher.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <vector>
#include <map>
#include <set>

class DisplacedGENMuonMatcher : public BaseMatcher
{
public:
  
  DisplacedGENMuonMatcher(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es);
  
  ~DisplacedGENMuonMatcher();

  const reco::GenParticle* getMatchedGENMuon() const {return matchedGENMuon_;} 
  std::vector<const reco::GenParticle*> getMatchedGENMuons() const {return matchedGENMuons_;}
  const reco::GenParticle* getMatchedDarkBoson() const {return matchedDarkBoson_;}

  int darkBosonIndex() const {return darkBosonIndex_;}
  int genMuonIndex() const {return genMuonIndex_;}
  bool checkRunOK() const {return runOK_;}

  double darkBosonDeltaR() const {return genGd0Gd1_m;}
  double darkBosonInvM() const {return genGd0Gd1_dR;}
  double matchedGenMudxy() const {return matchedGenMu_dxy;}
  double matchedGenMudR() const {return matchedGenMu_dR;}

  void testprint() const ;

private:

  void init();

  std::vector<edm::InputTag> input_;
  int verbose_;
  bool run_;
  enum {MuonGun, DarkSUSY};
  int sampleType_;
 
  void matchDisplacedGENMuonMatcherToSimTrack(const reco::GenParticleCollection& genParticles);
  void matchDisplacedGENMuonFromMuonGunMatcherToSimTrack(const reco::GenParticleCollection& genParticles);

  /* bool PtOrder (const reco::GenParticle* p1, const reco::GenParticle* p2); */
  double dxy(double px, double py, double vx, double vy, double pt);
  double phiHeavyCorr(double pt, double eta, double phi, double q);
  double invariantMass(const reco::Candidate* p1, const reco::Candidate* p2);

  Float_t beamSpot_x;
  Float_t beamSpot_y;
  Float_t beamSpot_z;

  const int nGlu = 2;
  const int nGd = 2;
  const int nGdMu = 2;

  Float_t genGlu_p[2];
  Float_t genGlu_pt[2];
  Float_t genGlu_px[2];
  Float_t genGlu_py[2];
  Float_t genGlu_pz[2];
  Float_t genGlu_eta[2];
  Float_t genGlu_phi[2];

  Float_t genH_m;
  Float_t genH_p;
  Float_t genH_pt;
  Float_t genH_px;
  Float_t genH_py;
  Float_t genH_pz;
  Float_t genH_eta;
  Float_t genH_phi;
  Float_t genH_vx;
  Float_t genH_vy;
  Float_t genH_vz;

  // Dark photon  
  Float_t genGd_m[2];
  Float_t genGd_E[2];
  Float_t genGd_p[2];
  Float_t genGd_pt[2];
  Float_t genGd_px[2];
  Float_t genGd_py[2];
  Float_t genGd_pz[2];
  Float_t genGd_eta[2];
  Float_t genGd_phi[2];
  Float_t genGd_vx[2];
  Float_t genGd_vy[2];
  Float_t genGd_vz[2];
  Float_t genGd_vLx[2];
  Float_t genGd_vLy[2];
  Float_t genGd_vLz[2];
  Float_t genGd_lxy[2];
  Float_t genGd_l[2];
  Float_t genGd_dxy[2];
  Float_t genGdMu_dxy_max[2];

  // Gen level muon
  Float_t genGdMu_q[2][2];
  Float_t genGdMu_p[2][2];
  Float_t genGdMu_pt[2][2];
  Float_t genGdMu_px[2][2];
  Float_t genGdMu_py[2][2];
  Float_t genGdMu_pz[2][2];
  Float_t genGdMu_eta[2][2];
  Float_t genGdMu_phi[2][2];
  Float_t genGdMu_phi_corr[2][2];
  Float_t genGdMu_vx[2][2];
  Float_t genGdMu_vy[2][2];
  Float_t genGdMu_vz[2][2];
  Float_t genGdMu_dxy[2][2];

  Float_t genGd0Gd1_m;
  Float_t genGd0Gd1_dR;

  Float_t matchedGenMu_dR;
  Float_t matchedGenMu_dxy;

  const reco::GenParticle* matchedGENMuon_;
  std::vector<const reco::GenParticle*> matchedGENMuons_;
  const reco::GenParticle* matchedDarkBoson_;

  bool runOK_;
  int darkBosonIndex_;
  int genMuonIndex_;
};

#endif
