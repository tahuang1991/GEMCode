import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMRECOANA")

## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

## global tag for 2019 upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

# the analyzer configuration
from GEMCode.GEMValidation.simTrackMatching_cfi import *
process.GEMRecHitAnalyzer = cms.EDAnalyzer("GEMRecHitAnalyzer",
    simTrackMatching = SimTrackMatching
)
## simhits and rechits
process.GEMRecHitAnalyzer = noDigiCollections(process.GEMRecHitAnalyzer)
process.GEMRecHitAnalyzer = noL1Collections(process.GEMRecHitAnalyzer)
process.GEMRecHitAnalyzer = noTrackCollections(process.GEMRecHitAnalyzer)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:out_rechit.root')
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("gem_localrec_ana.root")
)

process.p = cms.Path(process.GEMRecHitAnalyzer)

## messages                                                                                                                                                   
print
print 'Input files:'
print '----------------------------------------'
print process.source.fileNames
print
print 'Output file:'
print '----------------------------------------'
print process.TFileService.fileName
print
