import FWCore.ParameterSet.Config as cms

fileNamesPU = cms.untracked.vstring()

f = open('fileNamesPU_MinBias_14TeV_91X_Phase2.txt')
for line in f:
    fileNamesPU.extend((line))
  
