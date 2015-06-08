import FWCore.ParameterSet.Config as cms

process = cms.Process("FILTER")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:/eos/uscms/store/user/lpcgem/dildick/DisplacedMuGEMCSCILT/dildick/DarkSUSY_mH_125_mGammaD_0400_ctau_00_14TeV_madgraph452_bridge224_LHE_pythia6_GEN_SIM/DarkSUSY_mH_125_mGammaD_0400_ctau_00_14TeV_madgraph452_bridge224_LHE_pythia6_DIGI_L1/75d0e84afc642f0aede2fcf263b4fa2e/out_L1_1_1_O8R.root'
    )
)

process.USER = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring( 'p')
    ),
    fileName = cms.untracked.string('out_genfilter.root')
)

process.muJetFilter = cms.EDFilter('MuJetFilter')

process.p = cms.Path(process.muJetFilter)

process.outpath = cms.EndPath(process.USER)













