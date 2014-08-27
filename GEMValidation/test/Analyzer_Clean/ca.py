import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMCSCANA")

## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2019_cff')

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:out_L1.root')
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

## input
from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import *
from GEMCode.GEMValidation.InputFileHelpers import *
suffix = 'PU140_TMB_Baseline_2'
process = useInputDir(process, eosfiles[suffix], True)
events=100000
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(events)
)
outputFileName = 'hp_CMSSW14_usingL1SLHC7_'   + suffix + '_eff.root'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFileName)
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

##Change to handle old partition geometry for SLHC7 and arround

#from Geometry.GEMGeometry.gemGeometryCustoms import custom_GE11_6partitions_v1
#process = custom_GE11_6partitions_v1(process)

## global tag for upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
# the analyzer configuration
def enum(*sequential, **named):
  enums = dict(zip(sequential, range(len(sequential))), **named)
  return type('Enum', (), enums)
Stations = enum('ALL','ME11','ME1a','ME1b','ME12','ME13','ME21','ME22','ME31','ME32','ME41','ME42')

from GEMCode.GEMValidation.simTrackMatching_cfi import SimTrackMatching
process.GEMCSCAnalyzer = cms.EDAnalyzer("GEMCSCAnalyzer",
    verbose = cms.untracked.int32(0),
    stationsToUse = cms.vint32(Stations.ME11,Stations.ME1a,Stations.ME1b,
                               Stations.ME21,Stations.ME31,Stations.ME41),
    simTrackMatching = SimTrackMatching
)
matching = process.GEMCSCAnalyzer.simTrackMatching
matching.simTrack.minPt = 1.5
matching.gemRecHit.input = ""
"""
matching.cscTfTrack.input = ""
matching.tfCand.input = ""
matching.gmtCand.input = ""
matching.l1Extra.input = ""
"""
doGem = False
if doGem:
  matching.cscSimHit.minNHitsChamber = 3
  matching.cscStripDigi.minNHitsChamber = 3
  matching.cscWireDigi.minNHitsChamber = 3
  matching.cscCLCT.minNHitsChamber = 3
  matching.cscALCT.minNHitsChamber = 3
  matching.cscLCT.minNHitsChamber = 3
  matching.cscLCT.matchAlctGem = True
  matching.cscMPLCT.minNHitsChamber = 3


process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.p = cms.Path(process.GEMCSCAnalyzer)
mat = process.GEMCSCAnalyzer.simTrackMatching
mat.simTrack.minEta = cms.double(1.4)
mat.me0Segment.run = cms.bool(False)
mat.me0Muon.run = cms.bool(False)
mat.me0SimHit.run = cms.bool(False)
mat.me0SimHit.simMuOnly = cms.bool(False)
mat.me0RecHit.run = cms.bool(False)


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
