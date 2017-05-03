# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step2 --conditions auto:phase2_realistic -n 10 --era Phase2C2_timing --eventcontent FEVTDEBUGHLT -s DIGI:pdigi_valid,L1 --datatier GEN-SIM-DIGI --geometry Extended2023D4 --python DigiFullPU_2023tiltedPU.py --no_exec --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2C2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.simMuonME0ReDigis384 = process.simMuonME0ReDigis.clone(
    numberOfStrips = cms.uint32(384)
)
process.simMuonME0ReDigis192 = process.simMuonME0ReDigis.clone(
    numberOfStrips = cms.uint32(192)
)
process.simMuonME0ReDigis96 = process.simMuonME0ReDigis.clone(
    numberOfStrips = cms.uint32(96)
)

process.RandomNumberGeneratorService.simMuonME0ReDigis384 = cms.PSet(
    initialSeed = cms.untracked.uint32(1234567),
    engineName = cms.untracked.string('HepJamesRandom')
)
process.RandomNumberGeneratorService.simMuonME0ReDigis192 = cms.PSet(
    initialSeed = cms.untracked.uint32(2234567),
    engineName = cms.untracked.string('HepJamesRandom')
)
process.RandomNumberGeneratorService.simMuonME0ReDigis96 = cms.PSet(
    initialSeed = cms.untracked.uint32(3234567),
    engineName = cms.untracked.string('HepJamesRandom')
)

process.me0RecHits384 = process.me0RecHits.clone(
    me0DigiLabel = cms.InputTag("simMuonME0ReDigis384")
)
process.me0RecHits192 = process.me0RecHits.clone(
    me0DigiLabel = cms.InputTag("simMuonME0ReDigis192")
)
process.me0RecHits96 = process.me0RecHits.clone(
    me0DigiLabel = cms.InputTag("simMuonME0ReDigis96")
)

process.me0Segments384 = process.me0Segments.clone(
    me0RecHitLabel = cms.InputTag("me0RecHits384")
)
process.me0Segments192 = process.me0Segments.clone(
    me0RecHitLabel = cms.InputTag("me0RecHits192")
)
process.me0Segments96 = process.me0Segments.clone(
    me0RecHitLabel = cms.InputTag("me0RecHits96")
)


process.me0DigiRecoSequence = cms.Sequence(

    process.simMuonME0ReDigis384 *
    process.simMuonME0ReDigis192 *
    process.simMuonME0ReDigis96 *

    process.me0RecHits384 *
    process.me0RecHits192 *
    process.me0RecHits96 *

    process.me0Segments384 *
    process.me0Segments192 *
    process.me0Segments96
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/lpcgem/SingleMu_91X_FlatPt05_30_eta20_28_phase2_realistic_Extended2023D4_GEN_SIM/SingleMu_91X_FlatPt05_30_eta20_28_phase2_realistic_Extended2023D4_GEN_SIM/170425_023359/0000/step1_142.root'),
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
        dataTier = cms.untracked.string('GEN-SIM-DIGI'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(10485760),
    fileName = cms.untracked.string('file:step2.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*Muon*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_me0*_*_*')

# Additional output definition

# Other statements
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator + process.me0DigiRecoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
