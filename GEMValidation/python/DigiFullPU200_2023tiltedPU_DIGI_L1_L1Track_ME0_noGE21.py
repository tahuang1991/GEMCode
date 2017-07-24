# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step2 -n 10 --era Phase2C2_timing -s DIGI:pdigi_valid,L1 --conditions 91X_upgrade2023_realistic_v3 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --pileup AVE_200_BX_25ns --nThreads 4 --geometry Extended2023D17 --customise=GEMCode/GEMValidation/customizeL1TrackME0ForMuonTDR.customizeL1TrackME0ForMuonTDR --python DigiFullPU200_2023tiltedPU.py --no_exec --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2C2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:step1.root'),
    inputCommands = cms.untracked.vstring('keep *',
        'drop *_genParticles_*_*',
        'drop *_genParticlesForJets_*_*',
        'drop *_kt4GenJets_*_*',
        'drop *_kt6GenJets_*_*',
        'drop *_iterativeCone5GenJets_*_*',
        'drop *_ak4GenJets_*_*',
        'drop *_ak7GenJets_*_*',
        'drop *_ak8GenJets_*_*',
        'drop *_ak4GenJetsNoNu_*_*',
        'drop *_ak8GenJetsNoNu_*_*',
        'drop *_genCandidatesForMET_*_*',
        'drop *_genParticlesForMETAllVisible_*_*',
        'drop *_genMetCalo_*_*',
        'drop *_genMetCaloAndNonPrompt_*_*',
        'drop *_genMetTrue_*_*',
        'drop *_genMetIC5GenJs_*_*'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(10485760),
    fileName = cms.untracked.string('file:step2.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TTTrack*_Level1TTTracks_*')
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_L1TkMuons_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*Muon*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_me0*_*_*')

# don't keep the track clusters and stubs
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTCluster*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTStub*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simEcal*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simHcal*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_mix_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiStripDigis_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiPixelDigis_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simCalo*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simCalo*_TrackerHits*_*')

process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')

## customization provided by Sebastien Viret

## This crash is related to a geometry point. The branch you merged was
## based on a tracker geometry with 11 rings in the TIB1, whereas I suppose
## that in D17 you're using the last one, with 12 rings. Therefore the stub
## producer joboption is bad.

# Overwrite stub builder config
#
# Baseline Threshold
# 2Gev thresh 99% in the barrel flat and tilted except in TIB1 where 3GeV 99% and 95% in two last tilted rings
# 2GeV thresh 99% in the disks except in ring1 1 and 2 of disk 1 and 3 where 95% 3GeV
#
process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(True)
process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatching2S = cms.bool(True)
process.TTStubAlgorithm_official_Phase2TrackerDigi_.NTiltedRings = cms.vdouble( 0., 12., 12., 12., 0., 0., 0.) #Number of tilted rings per side in barrel layers (for tilted geom only)
process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 2., 2., 3, 4, 5, 6.5)
process.TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
    cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 2, 2, 2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1., 1.) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 2.5, 2.5, 2.5, 2.5, 2, 2, 2.5, 2.5, 2, 2, 2, 2) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 3.5, 3.5, 3, 3, 3, 3, 2.5, 2.5, 2.5, 2, 2, 2) ),
    )
process.TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
    cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2.5, 3, 3, 3.5, 4, 2.5, 3, 3.5, 4.5, 5.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2, 2.5, 3, 3, 3, 2, 3, 4, 5, 5.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 2, 2, 2.5, 2.5, 2.5, 3.5, 2.5, 5, 5.5, 6) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 2, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 1.5, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
    )

process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)

# L1 tracking
process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

  ## ME0 segments
process.simMuonME0ReDigis192 = process.simMuonME0ReDigis.clone(
    numberOfStrips = cms.uint32(192)
    )

process.RandomNumberGeneratorService.simMuonME0ReDigis192 = cms.PSet(
    initialSeed = cms.untracked.uint32(2234567),
    engineName = cms.untracked.string('HepJamesRandom')
    )

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.me0RecHits192 = process.me0RecHits.clone(
    me0DigiLabel = cms.InputTag("simMuonME0ReDigis192")
    )

process.me0Segments192 = process.me0Segments.clone(
    me0RecHitLabel = cms.InputTag("me0RecHits192")
    )

process.me0Segments192.algo_psets[1].algo_pset.maxPhiAdditional = cms.double(1.2*0.35/192)
process.me0Segments192.algo_psets[1].algo_pset.maxPhiSeeds = cms.double(1.2*0.35/192)

process.me0DigiRecoSequence = cms.Sequence(
    process.simMuonME0ReDigis192 *
    process.me0RecHits192 *
    process.me0Segments192
    )

process.load("L1Trigger.L1TTrackMatch.L1TkMuonProducer_cfi")
process.pL1TkMuon = cms.Path( process.L1TkMuons )

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '')

# Path and EndPath definitions
process.pdigi_valid = cms.Sequence(process.generatorSmeared+
                                   process.randomEngineStateProducer+
                                   process.mix+
                                   process.simMuonCSCDigis+
                                   process.simMuonDTDigis+
                                   process.simMuonRPCDigis+
                                   process.simMuonGEMDigis+
                                   process.simMuonGEMPadDigis+
                                   process.simMuonGEMPadDigiClusters+
                                   process.simMuonME0Digis+
                                   process.simMuonME0ReDigis+
                                   process.addPileupInfo)

process.digitisation_step = cms.Path(process.pdigi_valid)
process.simCscTriggerPrimitiveDigis.commonParam.runME21ILT = cms.bool(False)
process.L1simulation_step = cms.Path(process.SimL1TMuon*process.me0DigiRecoSequence)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,
                                process.L1TrackTrigger_step,
                                process.TTClusterStub,
                                process.TTTracksWithTruth,
                                process.L1simulation_step,
                                process.pL1TkMuon,
                                process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from GEMCode.GEMValidation.customizeL1TrackME0ForMuonTDR
from GEMCode.GEMValidation.customizeL1TrackME0ForMuonTDR import customizeL1TrackME0ForMuonTDR
from GEMCode.GEMValidation.addPileup import addPileup

process = addPileup(process)
#call to customisation function customizeL1TrackME0ForMuonTDR imported from GEMCode.GEMValidation.customizeL1TrackME0ForMuonTDR
#process = customizeL1TrackME0ForMuonTDR(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
