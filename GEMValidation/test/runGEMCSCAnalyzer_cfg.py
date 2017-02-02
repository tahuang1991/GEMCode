import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMCSCANA")

## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D1_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

process.source = cms.Source("PoolSource",
	#fileNames = cms.untracked.vstring('file:out_L1.root'),
	fileNames = cms.untracked.vstring('/store/user/dildick/DarkSUSY_mH_125_mGammaD_20000_cT_0_14TeV_GEN_SIM_90X/DarkSUSY_mH_125_mGammaD_20000_cT_0_14TeV_PU0_DIGI_L1/170116_230113/0000/step2_1.root')
)

from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import *
from GEMCode.GEMValidation.InputFileHelpers import *
InputDir = ['/eos/uscms/store/user/rdking/SingleMu1to100_14TeV_GEN_SIM_90X/DarkSUSY_mH_125_mGammaD_20000_cT_0_14TeV_PU140_DIGI_L1/170125_170819/0000/']
process = useInputDir(process, InputDir, True)
#process = useInputDir(process, files['_gem98_pt2-50_PU0_pt0_new'], False)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_ana.root")
)

## global tag for upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# the analyzer configuration
def enum(*sequential, **named):
  enums = dict(zip(sequential, range(len(sequential))), **named)
  return type('Enum', (), enums)
Stations = enum('ALL','ME11','ME1a','ME1b','ME12','ME13','ME21','ME22','ME31','ME32','ME41','ME42')

from GEMCode.GEMValidation.simTrackMatching_cfi import SimTrackMatching
process.GEMCSCAnalyzer = cms.EDAnalyzer("GEMCSCAnalyzer",
    verbose = cms.untracked.int32(1),
    stationsToUse = cms.vint32(Stations.ALL,Stations.ME11,Stations.ME1a,Stations.ME1b,
                              Stations.ME21,Stations.ME31,Stations.ME41),
    simTrackMatching = SimTrackMatching
)
matching = process.GEMCSCAnalyzer.simTrackMatching
matching.simTrack.minPt = 1.5
matching.matchprint = cms.bool(False)
"""
matching.gemRecHit.input = ""
matching.cscTfTrack.input = ""
matching.tfCand.input = ""
matching.gmtCand.input = ""
matching.l1Extra.input = ""
"""
doGem = True
if doGem:
  matching.cscSimHit.minNHitsChamber = 3
  matching.cscStripDigi.minNHitsChamber = 3
  matching.cscWireDigi.minNHitsChamber = 3
  matching.cscCLCT.minNHitsChamber = 3
  matching.cscALCT.minNHitsChamber = 3
  matching.cscLCT.minNHitsChamber = 3
  matching.cscMPLCT.minNHitsChamber = 3

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.p = cms.Path(process.GEMCSCAnalyzer)

## messages
print
#print 'Input files:'
print '----------------------------------------'
#print process.source.fileNames
#print
print 'Output file:'
print '----------------------------------------'
print process.TFileService.fileName
print
