import FWCore.ParameterSet.Config as cms
import os

fileNamesPU = cms.untracked.vstring()

## absolute path to txt file
path = os.getenv( "CMSSW_BASE" ) + "/src/GEMCode/GEMValidation/python/"

## load the official MinBias files for Muon TDR production
f = open(path + 'fileNamesPU_MinBias_14TeV_91X_Phase2.txt')
pu_files = f.read().split('\n')
f.close()

## add them to a list
pu_files = filter(lambda x: x.endswith('.root'),  pu_files)
fileNamesPU.extend(pu_files)
