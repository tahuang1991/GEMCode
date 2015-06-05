import FWCore.ParameterSet.Config as cms

process = cms.Process("FILTER")
process.maxEvents = cms.untracked.PSet(
	    output = cms.untracked.int32(100)
	)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:out_reco.root'
    )
)

process.USER = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_genfilter.root')
)

process.MuJetFilter = cms.EDFilter('MuJetFilter')

process.p = cms.Path(process.MuJetFilter)

process.outpath = cms.EndPath(process.USER)













