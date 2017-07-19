import FWCore.ParameterSet.Config as cms
#from GEMCode.GEMValidation.fileNamesPU_RelValMinBias_13_76X_mcRun2_asymptotic_v9 import fileNamesPU
#from GEMCode.GEMValidation.fileNamesPU_MinBias_14_82X_Phase2 import fileNamesPU ## TDR samples with old Geometry D4
from GEMCode.GEMValidation.fileNamesPU_MinBias_14TeV_91X_Phase2 import fileNamesPU ## TDR samples with Geometry D17

def addPileup(process):
    process.mix.input.fileNames = fileNamesPU
    return process
