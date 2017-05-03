import FWCore.ParameterSet.Config as cms

def randomizeMuonGun(process):
    ## filter for gen level
    process.genMuons = cms.EDFilter("PdgIdCandViewSelector",
        src = cms.InputTag("genParticles"),
        pdgId = cms.vint32(13, -13)
    )

    ## filter sequence
    process.gen_mu_select = cms.Sequence(process.genMuons)
    process.generation_step += process.gen_mu_select
    process.simulation_step += process.gen_mu_select
    return process
    
def randomizeMuonGunGEM(process):
    ## filter for gen level
    process.genMuons = cms.EDFilter("PdgIdCandViewSelector",
        src = cms.InputTag("genParticles"),
        pdgId = cms.vint32(13, -13)
    )

    ## filter for sim level
    process.genMuonsGEM = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("genMuons"),
        cut = cms.string("pt > 0.5 & abs(eta)<2.8 & abs(eta)>2.0"),   
        filter = cms.bool(True)
    )

    ## filter sequence
    process.gen_mu_select = cms.Sequence(process.genMuons * process.genMuonsGEM)
    process.generation_step += process.gen_mu_select
    process.simulation_step += process.gen_mu_select
    return process
