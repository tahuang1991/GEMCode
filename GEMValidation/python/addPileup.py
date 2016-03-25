import FWCore.ParameterSet.Config as cms
from GEMCode.GEMValidation.fileNamesPU_RelValMinBias_13_76X_mcRun2_asymptotic_v9 import fileNamesPU

def addPileup(process):
    process.mix.input.fileNames = fileNamesPU
    return process
