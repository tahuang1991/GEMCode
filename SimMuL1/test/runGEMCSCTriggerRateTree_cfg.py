import FWCore.ParameterSet.Config as cms
import os

## initialization
process = cms.Process('GEMCSCTRGANA')

## CMSSW RELEASE
cmssw = os.getenv( "CMSSW_VERSION" )

## steering
pileup = 140
events = 1000000
#sample = 'dimu'
sample = 'minbias'
globalTag = 'upgradePLS3'
suffix = 'rate'

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.Digi_cff')
process.load("Configuration.StandardSequences.SimL1Emulator_cff")
process.load("Configuration.StandardSequences.L1Extra_cff")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

## upgrade CSC geometry 
from SLHCUpgradeSimulations.Configuration.muonCustoms import unganged_me1a_geometry
process = unganged_me1a_geometry(process)

## upgrade CSC digitizer
from SLHCUpgradeSimulations.Configuration.muonCustoms import digitizer_timing_pre3_median
process = digitizer_timing_pre3_median(process)

## upgrade CSC L1 customizations
from SLHCUpgradeSimulations.Configuration.muonCustoms import customise_csc_L1Stubs 
process = customise_csc_L1Stubs(process)

## upgrade CSC TrackFinder
from SLHCUpgradeSimulations.Configuration.muonCustoms import customise_csc_L1TrackFinder
process = customise_csc_L1TrackFinder(process)

## upgrade L1Extra step
from SLHCUpgradeSimulations.Configuration.muonCustoms import customise_csc_L1Extra_allsim
process = customise_csc_L1Extra_allsim(process)
process.l1extraParticles.centralBxOnly = cms.bool(True)
process.l1extraParticles.produceMuonParticles = cms.bool(True)
process.l1extraParticles.produceCaloParticles = cms.bool(False)
process.l1extraParticles.ignoreHtMiss = cms.bool(False)

process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:out_L1_Nu_SLHC12_2023Muon_PU140_geonmo.root'
    )
)

## input
from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import *
from GEMCode.GEMValidation.InputFileHelpers import *
#process = useInputDir(process, files['_gem98_pt2-50_PU0_pt0_new'], False)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(events)
)

## analyzer - use unganged TF!!!
from GEMCode.GEMValidation.simTrackMatching_cfi import SimTrackMatching
from L1Trigger.CSCTrackFinder.csctfTrackDigisUngangedME1a_cfi import csctfTrackDigisUngangedME1a
process.GEMCSCTriggerRateTree = cms.EDAnalyzer("GEMCSCTriggerRateTree",
    simTrackMatching = SimTrackMatching,
    sectorProcessor = csctfTrackDigisUngangedME1a.SectorProcessor,
)

## output
outputFileName = 'hp_' + sample + "_" + cmssw + "_" + globalTag + "_pu%d"%(pileup) + '_w3_' + suffix + '.root'
print "outputFile:", outputFileName

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFileName)
)

## Sequence and schedule
process.ana_seq = cms.Sequence(process.GEMCSCTriggerRateTree)
process.ana_step = cms.Path(process.ana_seq)
process.schedule = cms.Schedule(process.ana_step)

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
