## pick your scenario:
## 1: 2019
## 2: 2019WithGem
## 3: 2023Muon

scenario = 3

## This configuration runs the DIGI+L1Emulator step
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("MUTRG")

## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
if scenario is 1 or scenario is 2:
    process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2019_cff')
elif scenario is 3:
    process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
else:
    print 'Something wrong with geometry'
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.SimL1Emulator_cff")
process.load("Configuration.StandardSequences.L1Extra_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
if scenario is 1 or scenario is 2:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
elif scenario is 3:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')    

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## calibration
from CalibMuon.CSCCalibration.CSCIndexer_cfi import CSCIndexerESProducer
process.CSCIndexerESProducer= CSCIndexerESProducer

from CalibMuon.CSCCalibration.CSCChannelMapper_cfi import CSCChannelMapperESProducer
process.CSCChannelMapperESProducer= CSCChannelMapperESProducer

## input commands
process.source = cms.Source("PoolSource",
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  inputCommands = cms.untracked.vstring('keep  *_*_*_*'),
  fileNames = cms.untracked.vstring('file:out_digi.root')
)

## input
from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import eosfiles
from GEMCode.GEMValidation.InputFileHelpers import useInputDir
dataset = '_Nu_SLHC12_2023Muon_PU140'
dataset = "_pt2-50_SLHC11_2023Muon_PU140"
dataset = '_Nu_SLHC12_2023Muon_PU140_geonmo'
process = useInputDir(process, eosfiles[dataset], True)


physics = False
if not physics:
    ## drop all unnecessary collections
    process.source.inputCommands = cms.untracked.vstring(
        'keep  *_*_*_*',
        'drop *_simCscTriggerPrimitiveDigis_*_*',
        'drop *_simDtTriggerPrimitiveDigis_*_*',
        'drop *_simRpcTriggerDigis_*_*',
        'drop *_simCsctfTrackDigis_*_*',
        'drop *_simDttfDigis_*_*',
        'drop *_simCsctfDigis_*_*',
        'drop *_simGmtDigis_*_*',
        'drop *_l1extraParticles_*_*'
        )
    
## output commands 
theOutDir = ''
theFileName = 'out_L1'+ dataset + '.root'
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(theOutDir + theFileName),
    outputCommands = cms.untracked.vstring('keep  *_*_*_*')
)

physics = False
if not physics:
    ## drop all unnecessary collections
    process.output.outputCommands = cms.untracked.vstring(
        'keep  *_*_*_*',
        # drop all CF stuff
        'drop *_mix_*_*',
        # drop tracker simhits
        'drop PSimHits_*_Tracker*_*',
        # drop calorimetry stuff
        'drop PCaloHits_*_*_*',
        'drop L1Calo*_*_*_*',
        'drop L1Gct*_*_*_*',
        # drop calorimetry l1extra
        'drop l1extraL1Em*_*_*_*',
        'drop l1extraL1Jet*_*_*_*',
        'drop l1extraL1EtMiss*_*_*_*',
        # clean up simhits from other detectors
        'drop PSimHits_*_Totem*_*',
        'drop PSimHits_*_FP420*_*',
        'drop PSimHits_*_BSC*_*',
        # drop some not useful muon digis and links
        'drop *_*_MuonCSCStripDigi_*',
        'drop *_*_MuonCSCStripDigiSimLinks_*',
        'drop *SimLink*_*_*_*',
        'drop *RandomEngineStates_*_*_*',
        'drop *_randomEngineStateProducer_*_*'
        )

## custom sequences
process.mul1 = cms.Sequence(
  process.SimL1MuTriggerPrimitives *
  process.SimL1MuTrackFinders *
  process.simRpcTriggerDigis *
  process.simGmtDigis *
  process.L1Extra
)

process.muL1Short = cms.Sequence(
  process.simCscTriggerPrimitiveDigis * 
  process.SimL1MuTrackFinders *
  process.simGmtDigis *
  process.L1Extra
)

## define path-steps
shortRun = False
if shortRun:
    process.L1simulation_step = cms.Path(process.muL1Short)
else: 
    process.L1simulation_step = cms.Path(process.mul1)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

## Schedule definition
process.schedule = cms.Schedule(
    process.L1simulation_step,
    process.endjob_step,
    process.out_step
)

## customization
if scenario is 1:
    from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019
    process = cust_2019(process)
elif scenario is 2:
    from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019WithGem
    process = cust_2019WithGem(process)
elif scenario is 3:
    from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muon
    process = cust_2023Muon(process)

## some extra L1 customs
process.l1extraParticles.centralBxOnly = cms.bool(True)
process.l1extraParticles.produceMuonParticles = cms.bool(True)
process.l1extraParticles.produceCaloParticles = cms.bool(False)
process.l1extraParticles.ignoreHtMiss = cms.bool(False)

## messages
print
print 'Input files:'
print '----------------------------------------'
print process.source.fileNames
print
print 'Output file:'
print '----------------------------------------'
print process.output.fileName
print 
