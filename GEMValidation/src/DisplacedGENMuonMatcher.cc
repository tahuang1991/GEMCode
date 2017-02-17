#include "GEMCode/GEMValidation/interface/DisplacedGENMuonMatcher.h"

DisplacedGENMuonMatcher::DisplacedGENMuonMatcher(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es)
: BaseMatcher(t, v, ps, ev, es)
{
  auto displacedGenMu_= conf().getParameter<edm::ParameterSet>("displacedGenMu");
  input_ = displacedGenMu_.getParameter<std::vector<edm::InputTag>>("validInputTags");
  verbose_ = displacedGenMu_.getParameter<int>("verbose");
  run_ = displacedGenMu_.getParameter<bool>("run");
  sampleType_ = displacedGenMu_.getParameter<int>("sampleType");

  edm::Handle<reco::GenParticleCollection> genParticles;
  runOK_ = false;
  bool getProduct = gemvalidation::getByLabel(input_, genParticles, event());
  if(getProduct and run_ and sampleType_ == DarkSUSY) matchDisplacedGENMuonMatcherToSimTrack(*genParticles.product());
  else if (getProduct and run_ and sampleType_ == MuonGun) matchDisplacedGENMuonFromMuonGunMatcherToSimTrack(*genParticles.product());
  
  matchedGENMuon_ = NULL;
  matchedDarkBoson_ = NULL;

}

DisplacedGENMuonMatcher::~DisplacedGENMuonMatcher()
{

}

void
DisplacedGENMuonMatcher::matchDisplacedGENMuonMatcherToSimTrack(const reco::GenParticleCollection& genParticles)
{
  double eq = 0.000001;

  // edm::Handle<reco::BeamSpot> beamSpot;
  // match_sh.event().getByLabel("offlineBeamSpot",beamSpot);
  beamSpot_x = 0;//beamSpot->position().x();
  beamSpot_y = 0;//beamSpot->position().y();
  beamSpot_z = 0;//beamSpot->position().z();

  // Loop over all genParticles and save prompt muons from particles with codes 36 (a1) or 3000022 (gammaD) in vector genMuons
  std::vector<const reco::GenParticle*> genGlu_unsorted;
  std::vector<const reco::GenParticle*> genGlu;
  std::vector<const reco::GenParticle*> genH;
  std::vector<const reco::GenParticle*> genGd_unsorted;
  std::vector<const reco::GenParticle*> genGd;
  std::vector<const reco::GenParticle*> genMuons;
  std::vector<const reco::Candidate*>   genMuonMothers;

  int counterGenParticle = 0;
  for(reco::GenParticleCollection::const_iterator iGenParticle = genParticles.begin();  iGenParticle != genParticles.end();  ++iGenParticle) {
    counterGenParticle++;
    // Check if gen particle is muon (pdgId = +/-13) and stable (status = 1)
    //    if (verbose_) std::cout << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;
    if ( fabs( iGenParticle->pdgId() ) == 13 and iGenParticle->status() == 1 ) {
      if (verbose_) std::cout << "Muon " << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;
      // Mother of the muon can be muon. Find the last muon in this chain: genMuonCand
      // Example: a1 -> mu+ (status = 3) mu- (status = 3)
      //          mu- (status = 3) -> mu- (status = 2) -> mu- (status = 1)
      const reco::Candidate *genMuonCand = &(*iGenParticle);
      bool isMuonMother = true;
      while(isMuonMother) {
        isMuonMother = false;
        for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
          if ( fabs( genMuonCand->mother(iMother)->pdgId() ) == 13 ) {
            isMuonMother = true;
            genMuonCand = genMuonCand->mother(iMother);
          }
        }
      }
      // Loop over all real (non-muon) mothers of the muon (here we use genMuonCand)
      for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
        // Check if mother is CP-odd Higgs (PdgId = 36) or gamma_Dark (PdgId = 3000022)
	//        if ( genMuonCand->mother(iMother)->pdgId() == 36 || genMuonCand->mother(iMother)->pdgId() == 3000022 || genMuonCand->mother(iMother)->pdgId() == 443 ) {
        if ( genMuonCand->mother(iMother)->pdgId() == 36 || genMuonCand->mother(iMother)->pdgId() == 3000022 ) {
          // Store the muon (stable, first in chain) into vector
          genMuons.push_back(&(*iGenParticle));
          // Store mother of the muon into vector. We need this to group muons into dimuons later
          genMuonMothers.push_back(genMuonCand->mother(iMother));
        }
      }
    }
    // Check if gen particle is
    if (    ( iGenParticle->status() == 3 and iGenParticle->pdgId() == 25 ) // decaying (status = 3) SM Higgs (pdgId = 25)
	    || ( iGenParticle->status() == 3 and iGenParticle->pdgId() == 35 ) // decaying (status = 3) CP-even Higgs (pdgId = 35)
	    ) {
      if (iGenParticle->pdgId() == 25) if (verbose_) std::cout << "SM Higgs ";
      if (iGenParticle->pdgId() == 35) if (verbose_) std::cout << "CP-even Higgs ";
      if (verbose_) std::cout << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;    
      genH.push_back(&(*iGenParticle)); // Store the Higgs into vector
    }
    // Check if gen particle is
    if (    ( iGenParticle->status() == 3 and iGenParticle->pdgId() == 36      ) // decaying (status = 3) CP-odd Higgs (pdgId = 36)
	    || ( iGenParticle->status() == 3 and iGenParticle->pdgId() == 3000022 ) // decaying (status = 3) gamma_Dark (pdgId = 3000022)
	    //         || ( iGenParticle->status() == 2 and iGenParticle->pdgId() == 443   ) // decaying (status = 2) J/psi (pdgId = 443)
	    ) {
      if (iGenParticle->pdgId() == 36) if (verbose_) std::cout << "SM Higgs ";
      if (iGenParticle->pdgId() == 3000022) if (verbose_) std::cout << "CP-even Higgs ";
      if (verbose_) std::cout << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;   
      genGd_unsorted.push_back(&(*iGenParticle));
    }
    // Check if gen particle is gluon 
    if (iGenParticle->status() == 3 and iGenParticle->pdgId() == 21){
      if (verbose_) std::cout << "Gluon " << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;    
      genGlu_unsorted.push_back(&(*iGenParticle));
    }
  }
  if ( genH.size() == 1 ) {
    genH_m   = genH[0]->mass();
    genH_p   = genH[0]->p();
    genH_pt  = genH[0]->pt();
    genH_px  = genH[0]->px();
    genH_py  = genH[0]->py();
    genH_pz  = genH[0]->pz();
    genH_eta = genH[0]->eta();
    genH_phi = genH[0]->phi();
    genH_vx  = genH[0]->vx() - beamSpot_x;
    genH_vy  = genH[0]->vy() - beamSpot_y;
    genH_vz  = genH[0]->vz() - beamSpot_z;
  } else {
    if (verbose_) std::cout << "WARNING! genH.size() != 1" << std::endl;
  }
  bool sortGammaDark(false);
  if (sortGammaDark) {
    if ( genGd_unsorted.size() >= 2 ) {
      // Sort genGd by pt (leading pt first)
      std::sort (genGd_unsorted.begin(), genGd_unsorted.end(), gemvalidation::PtOrder);
      // Remove duplicates from genGd
      //    Float_t A_pt = genGd_unsorted[0]->pt();
      //    for ( unsigned int i = 1; i < genGd_unsorted.size(); i++ ) {
      //      if ( fabs( genGd_unsorted[i]->pt() - A_pt) > eq ) {
      //        A_pt = genGd_unsorted[i]->pt();
      //        genGd.push_back( genGd_unsorted[i] );
      //      }
      //    }
    }
  }
  genGd = genGd_unsorted;
  // check again that the gluons are in fact the mother particles
  if ( genGlu_unsorted.size() >= 2 ) {
    // Sort genGlu by pt (leading pt first)
    std::sort (genGlu_unsorted.begin(), genGlu_unsorted.end(), gemvalidation::PtOrder);
  }

  genGlu = genGlu_unsorted;

  if ( genGd.size() >= 2 ) {
    for (int i=0; i<2; ++i){
      genGd_m[i]   = genGd[i]->mass();
      genGd_E[i]   = genGd[i]->energy();
      genGd_p[i]   = genGd[i]->p();
      genGd_pt[i]  = genGd[i]->pt();
      genGd_px[i]  = genGd[i]->px();
      genGd_py[i]  = genGd[i]->py();
      genGd_pz[i]  = genGd[i]->pz();
      genGd_eta[i] = genGd[i]->eta();
      genGd_phi[i] = genGd[i]->phi();
      genGd_vx[i]  = genGd[i]->vx() - beamSpot_x;
      genGd_vy[i]  = genGd[i]->vy() - beamSpot_y;
      genGd_vz[i]  = genGd[i]->vz() - beamSpot_z;
    }
  } else {
    if (verbose_) std::cout << "WARNING! genGd.size() < 2" << std::endl;
  }
  if ( genGlu.size() >= 2 ) {    
    for (int i=0; i<2; ++i){
      genGlu_p[i] = genGlu[i]->p();
      genGlu_pt[i]  = genGlu[i]->pt();
      genGlu_px[i]  = genGlu[i]->px();
      genGlu_py[i]  = genGlu[i]->py();
      genGlu_pz[i]  = genGlu[i]->pz();
      genGlu_eta[i] = genGlu[i]->eta();
      genGlu_phi[i] = genGlu[i]->phi();
    }
  } else {
    if (verbose_) std::cout << "WARNING! genGlu.size() < 2" << std::endl;
  }
  // Group muons with the same mother
  std::vector< std::vector<const reco::GenParticle*> > genMuonGroupsUnsorted;
  std::vector<const reco::Candidate*> genMuonGroupsUnsortedMothers;
  std::vector<const reco::GenParticle*> genMuonsTMP1       = genMuons;
  std::vector<const reco::Candidate*>   genMuonMothersTMP1 = genMuonMothers;
  unsigned int nMuonGroup = 0;
  while ( genMuonsTMP1.size() > 0 ) {
    std::vector<const reco::GenParticle*> genMuonsTMP2;
    std::vector<const reco::Candidate*>   genMuonMothersTMP2;
    std::vector<const reco::GenParticle*> genMuonsSameMother;
    std::vector<const reco::Candidate*>   genMuonMothersSame;
    for ( unsigned int j = 0; j < genMuonsTMP1.size(); j++ ) {
      // Check if mothers are the same particle
      if ( fabs( genMuonMothersTMP1[0]->pt() - genMuonMothersTMP1[j]->pt() ) < eq ) {
        genMuonsSameMother.push_back( genMuonsTMP1[j] );
      } else {
        genMuonsTMP2.push_back( genMuonsTMP1[j] );
        genMuonMothersTMP2.push_back( genMuonMothersTMP1[j] );
      }
    }
    genMuonGroupsUnsorted.push_back(genMuonsSameMother);
    genMuonGroupsUnsortedMothers.push_back(genMuonMothersTMP1[0]);
    genMuonsTMP1       = genMuonsTMP2;
    genMuonMothersTMP1 = genMuonMothersTMP2;
    nMuonGroup++;
  }

  // Sort muon groups to match order of genGd vector
  std::vector< std::vector<const reco::GenParticle*> > genMuonGroups;
  std::vector<const reco::Candidate*> genMuonGroupsMothers;
  for (unsigned int iA = 0; iA < genGd.size(); iA++ ) {
    bool isMuGroupMatchedToA = false;
    int  nMuGroup = -1;
    for ( unsigned int iMuGroup = 0; iMuGroup < genMuonGroupsUnsortedMothers.size(); iMuGroup++ ) {
      if ( fabs ( genGd[iA]->pt() - genMuonGroupsUnsortedMothers[iMuGroup]->pt() ) < 0.000001 ) {
        isMuGroupMatchedToA = true;
        nMuGroup = iMuGroup;
        break;
      }
    }
    if ( isMuGroupMatchedToA and nMuGroup >= 0 ) {
      genMuonGroups.push_back( genMuonGroupsUnsorted[nMuGroup] );
      genMuonGroupsMothers.push_back( genMuonGroupsUnsortedMothers[nMuGroup] );
    } else {
      if (verbose_) std::cout << "Error! Muon group has no matched boson A" << std::endl;
    }
  }
  
  float minDeltaR = 99.0;
  int index1=-99;
  int index2=-99;
  if ( genMuonGroups.size() == 2 and genMuonGroups[0].size() == 2 and genMuonGroups[1].size() == 2 ) {
    std::sort( genMuonGroups[0].begin(), genMuonGroups[0].end(), gemvalidation::PtOrder );
    std::sort( genMuonGroups[1].begin(), genMuonGroups[1].end(), gemvalidation::PtOrder );
    
    for (int i=0; i<2; ++i){ 
      for (int j=0; j<2; ++j){
        genGdMu_q[i][j] = genMuonGroups[i][j]->charge();
        genGdMu_p[i][j] = genMuonGroups[i][j]->p();
        genGdMu_pt[i][j] = genMuonGroups[i][j]->pt();
        genGdMu_px[i][j] = genMuonGroups[i][j]->px();
        genGdMu_py[i][j] = genMuonGroups[i][j]->py();
        genGdMu_pz[i][j] = genMuonGroups[i][j]->pz();
        genGdMu_eta[i][j] = genMuonGroups[i][j]->eta();
        genGdMu_phi[i][j] = genMuonGroups[i][j]->phi();
        genGdMu_phi_corr[i][j] = phiHeavyCorr(genMuonGroups[i][j]->pt(), 
					      genMuonGroups[i][j]->eta(), 
					      genMuonGroups[i][j]->phi(), 
					      genMuonGroups[i][j]->charge());
        genGdMu_vx[i][j] = genMuonGroups[i][j]->vx();
        genGdMu_vy[i][j] = genMuonGroups[i][j]->vy();
        genGdMu_vz[i][j] = genMuonGroups[i][j]->vz();
	const double deltaR(reco::deltaR(genGdMu_eta[i][j], genGdMu_phi[i][j], trk().momentum().eta(), trk().momentum().phi()));
	if (deltaR < minDeltaR) {
	  minDeltaR = deltaR;
	  index1 = i;
	  index2 = j;
	}	
      }
    }
  }

  if (index1 == -99 or index2 == -99 ){
        if (verbose_) std::cout <<" failed to find matched dark boson or matched gen muon " << std::endl;
      	runOK_ = false;
  	return; 
  }
  if (verbose_) std::cout <<" Find matched dark boson and matched gen muon index1 "<< index1 << " index2 "<< index2 << std::endl;
  runOK_ = true;

  if (verbose_){
    for (int i=0; i<2; ++i){ 
      for (int j=0; j<2; ++j){
	if (i==index1 and j==index2) std::cout <<"=======  this is gen Mu matched to sim Mu ======= "<< std::endl;
	std::cout << "genGdMu_pt_" << i << j << ": " << genGdMu_pt[i][j] << std::endl;
	std::cout << "genGdMu_eta_" << i << j << ": " << genGdMu_eta[i][j] << std::endl;
	std::cout << "genGdMu_phi_" << i << j << ": " << genGdMu_phi[i][j] << std::endl;
	std::cout << "genGdMu_phi_corr_" << i << j << ": " << genGdMu_phi_corr[i][j] << std::endl;
	std::cout << "deltaR (" << i << j << ", simTrack): "<< reco::deltaR(genGdMu_eta[i][j], genGdMu_phi[i][j], trk().momentum().eta(), trk().momentum().phi()) << std::endl;
      }
    }
    std::cout << "genGd0Gd1_m " << genGd0Gd1_m << std::endl;
    std::cout << "genGd0Gd1_dR " << genGd0Gd1_dR << std::endl;
  }

  genGd0Gd1_m = invariantMass(genGd[0],genGd[1]);
  genGd0Gd1_dR = reco::deltaR(genGd_eta[0], genGd_phi[0], genGd_eta[1], genGd_phi[1]);
  darkBosonIndex_ = index1;
  genMuonIndex_ = index2;
  matchedGENMuon_ = genMuonGroups[index1][index2];
  matchedGenMu_dR = reco::deltaR(genGdMu_eta[index1][index2], genGdMu_phi[index1][index2], trk().momentum().eta(), trk().momentum().phi());

  matchedGenMu_dxy = dxy(genGdMu_px[index1][index2], genGdMu_py[index1][index2], genGdMu_vx[index1][index2], genGdMu_vy[index1][index2], genGdMu_pt[index1][index2]);
  //matchedGenMu_dxy = dxy(matchedGENMuon_->px(), matchedGENMuon_->py(), matchedGENMuon_->vx(), matchedGENMuon_->vy(), matchedGENMuon_->pt());

  matchedGENMuons_.clear();
  matchedGENMuons_.push_back(genMuonGroups[index1][0]);
  matchedGENMuons_.push_back(genMuonGroups[index1][1]);
  matchedDarkBoson_ = genGd[index1];
}


void
DisplacedGENMuonMatcher::matchDisplacedGENMuonFromMuonGunMatcherToSimTrack(const reco::GenParticleCollection& genParticles)
{

  beamSpot_x = 0;//beamSpot->position().x();
  beamSpot_y = 0;//beamSpot->position().y();
  beamSpot_z = 0;//beamSpot->position().z();

  float minDeltaR = 99;
  int counterGenParticle = 0;
  for(reco::GenParticleCollection::const_iterator iGenParticle = genParticles.begin();  iGenParticle != genParticles.end();  ++iGenParticle) {
    counterGenParticle++;
    // Check if gen particle is muon (pdgId = +/-13) and stable (status = 1)
    //    if (verbose_) std::cout << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;
    if ( fabs( iGenParticle->pdgId() ) == 13 and iGenParticle->status() == 1 ) {
      if (verbose_) std::cout << "Muon " << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;
	const double deltaR(reco::deltaR(iGenParticle->eta(), iGenParticle->phi(), trk().momentum().eta(), trk().momentum().phi()));
	if (deltaR < minDeltaR) {
	  minDeltaR = deltaR;
	  matchedGENMuon_ = iGenParticle->clone();
	  matchedGenMu_dR = deltaR;
	}
    }


  }

}
double 
DisplacedGENMuonMatcher::dxy(double px, double py, double vx, double vy, double pt)
{
  //Source: https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrackReco/interface/TrackBase.h#119
  return (- vx * py + vy * px ) / pt;
}

double 
DisplacedGENMuonMatcher::phiHeavyCorr(double pt, double eta, double phi, double q)
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
DisplacedGENMuonMatcher::invariantMass(const reco::Candidate* p1, const reco::Candidate* p2) 
{
  return  sqrt((p1->energy() + p2->energy())*(p1->energy() + p2->energy()) -
               (p1->px() + p2->px())*(p1->px() + p2->px()) -
               (p1->py() + p2->py())*(p1->py() + p2->py()) -
               (p1->pz() + p2->pz())*(p1->pz() + p2->pz()) );
}


void DisplacedGENMuonMatcher::testprint() const {

	std::cout <<"here is a test for displacedGenMuon matcher "<< std::endl;

}
