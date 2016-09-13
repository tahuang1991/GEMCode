# Instructions for GEM-CSC trigger development studies targeting Muon Upgrade TDR 2017 

Configurations for upgrade simulations with GEMs can be found here: 
https://twiki.cern.ch/twiki/bin/viewauth/CMS/GemSimulationsInstructionsCMSSW#Recipe_for_GEMs_in_latest_8XY_re

## Step 1: Setup ##
<PRE>
cmsrel CMSSW_8_1_0_pre11
cd CMSSW_8_1_0_pre11/src
cmsenv
git cms-init
git cms-merge-topic dildick:https://github.com/dildick/cmssw/tree/from-CMSSW_8_1_0_pre11-gem-trigger-cleanup
git cms-addpkg L1Trigger/CSCTriggerPrimitives
git clone git@github.com:gem-sw/GEMCode.git
cd GEMCode
git checkout -b for-CMSSW-81X origin/for-CMSSW-81X
scram b -j 9
</PRE>

## Step 2: GEN-SIM ##
<PRE>
cmsDriver.py SingleMuPt100_cfi \
--conditions auto:run2_mc -n 10 \
--era Phase2C1 \
--eventcontent FEVTDEBUG  \
-s GEN,SIM --datatier GEN-SIM \
--beamspot Realistic50ns13TeVCollision \
--customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted \
--geometry Extended2023D1 \
--python SingleMuPt100_2023tilted_GenSimFull.py \
--no_exec --fileout file:step1.root
</PRE>

## Step 3: DIGI-L1 ##
<PRE>
cmsDriver.py step2 \
--conditions auto:run2_mc \
--pileup_input das:/RelValMinBias_TuneZ2star_14TeV/1/GEN-SIM \
-n 1 --era Phase2C1 \
--eventcontent FEVTDEBUGHLT \
-s DIGI:pdigi_valid,L1 \
--datatier GEN-SIM-DIGI \
--pileup AVE_200_BX_25ns \
--customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted \
--geometry Extended2023D1 \
--python DigiFullPU_2023tiltedPU.py \
--no_exec \
--filein file:step1.root \
--fileout file:step2.root
</PRE>

When running without pileup, remove the pileup related flags. 

## Step 4: Analysis ##
Efficiency:
<PRE>
cmsRun runMuonUpgradeTDREfficiency_cfg.py
</PRE>
Rate:
<PRE>
cmsRun runMuonUpgradeTDRRate_cfg.py
</PRE>
