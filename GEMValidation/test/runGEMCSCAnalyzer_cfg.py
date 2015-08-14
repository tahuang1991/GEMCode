import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMCSCANA")

## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:out_hlt.test.root'),
)

from GEMCode.GEMValidation.InputFileHelpers import *
from L1Trigger.L1IntegratedMuonTrigger.hltSamples import *
"""
eosfiles["cT_100_14TeV_PU140_HLT_fullScope"] 
eosfiles["cT_100_14TeV_PU140_HLT_descope235MCHF"] 
eosfiles["cT_100_14TeV_PU140_HLT_descope235MCHFaging"] 
eosfiles["cT_100_14TeV_PU140_HLT_descope200MCHF"] 
eosfiles["cT_100_14TeV_PU140_HLT_descope200MCHFaging"] 

eosfiles["cT_1000_14TeV_PU140_HLT_fullScope"] 
eosfiles["cT_1000_14TeV_PU140_HLT_descope235MCHFaging"] 
eosfiles["cT_1000_14TeV_PU140_HLT_descope200MCHF"] 
eosfiles["cT_1000_14TeV_PU140_HLT_descope200MCHFaging"] 
"""


process = useInputDir(process, eosfiles["cT_100_14TeV_PU140_HLT_fullScope"])
print process.source.fileNames

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_ana_DarkSUSY_mH_125_mGammaD_20000_cT_100_14TeV_PU140_HLT_fullScope.root")
)

## global tag for upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

from GEMCode.GEMValidation.simTrackMatching_cfi import SimTrackMatching
process.GEMCSCAnalyzer = cms.EDAnalyzer("GEMCSCAnalyzer",
    verbose = cms.untracked.int32(0),
    simTrackMatching = SimTrackMatching
)
matching = process.GEMCSCAnalyzer.simTrackMatching
matching.simTrack.minPt = 1.5
matching.matchprint = cms.bool(False)

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
doRpc = False
if doRpc:
  matching.cscLCT.matchAlctRpc = True

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.p = cms.Path(process.GEMCSCAnalyzer)

## messages
print
#print 'Input files:'
print '----------------------------------------'
#print process.source.fileNames
print
print 'Output file:'
print '----------------------------------------'
print process.TFileService.fileName
print
