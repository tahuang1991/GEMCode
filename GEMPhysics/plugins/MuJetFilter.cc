/*
 * Muon-jet filter:
 * - 2 muons in barrel
 * - 2 muons in endcap
 * - all muon pT > 10GeV
 */

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;

class MuJetFilter : public edm::EDFilter {
   public:
      explicit MuJetFilter(const edm::ParameterSet&);
      ~MuJetFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
};

MuJetFilter::MuJetFilter(const edm::ParameterSet& iConfig)
{
}


MuJetFilter::~MuJetFilter()
{
}

bool MuJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  
  std::vector<const reco::GenParticle*> genMuons;
  std::vector<const reco::Candidate*>   genMuonMothers;
  // Loop over all gen particles
  int counterGenParticle = 0;
  for(auto iGenParticle = genParticles->begin();  iGenParticle != genParticles->end();  ++iGenParticle) {
    counterGenParticle++;
    //    std::cout << counterGenParticle << " " << iGenParticle->status() << " " << iGenParticle->pdgId() << " " << iGenParticle->vx() << " " << iGenParticle->vy() << " " << iGenParticle->vz() << std::endl;
    // Check if gen particle is muon (pdgId = +/-13) and stable (status = 1)
    if (fabs(iGenParticle->pdgId()) == 13 && iGenParticle->status() == 1) {
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
  }
  
  // 2 muons in barrel, 2 muons in endcap
  int nMuBarrel = 0;
  int nMuEndcap = 0;
  int nMuPTok = 0;
  for (auto muon: genMuons) {
    const float eta(fabs(muon->eta()));
    if (eta < 1.1) nMuBarrel++;
    if (eta > 1.1 and eta < 2.4) nMuEndcap++;
    if (muon->pt() > 10) nMuPTok++;
  }

  return nMuBarrel>=2 and nMuEndcap>=2 and nMuPTok>=1;
}

void MuJetFilter::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuJetFilter);
