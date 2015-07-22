#runevents = 1000
runevents = -1;

import FWCore.ParameterSet.Config as cms
process = cms.Process('L1analysis')

process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuGMTScalesConfig_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring("file:out_L1.root")
)
#out_L1muon2023GE21.root
fileOutputName = "Neutrino_SLHC25_PU140_0706"
#input: L1 sample
from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import *
from GEMCode.GEMValidation.InputFileHelpers import *
InputDir = ['/eos/uscms/store/user/tahuang/SLHC25_patch1_2023Neutrino_100k_L1_PU140/tahuang/SLHC25_patch1_2023Neutrino_gen_sim_100k/SLHC25_patch1_2023Neutrino_100k_L1_PU140/1bf93df4dfbb43dc918bd6e47dedbf79/']
#InputDir = ['/eos/uscms/store/user/tahuang/SLHC25_patch1_2023Muon_1M_L1_PU140_Pt2_50/tahuang/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/SLHC25_patch1_2023Muon_1M_L1_PU140_Pt2_50/1bf93df4dfbb43dc918bd6e47dedbf79/']
process = useInputDir(process, InputDir, True)

histofileName= "Rate_"+fileOutputName+".root"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(runevents)
)

#print "fileNames: ", process.source.fileNames
process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleMuPt100_cfi nevts:100'),
    name = cms.untracked.string('Applications')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
    histofileName
))

process.L1TTriggerRate = cms.EDAnalyzer('L1MuonTTriggerRate',
        minPt = cms.double(2.0),
        maxPt = cms.double(100.0),
        minEta = cms.double(1.6),
        maxEta = cms.double(2.4),
    )
process.pL1TAnalyser = cms.Path(process.L1TTriggerRate)
    
process.schedule = cms.Schedule(process.pL1TAnalyser)
