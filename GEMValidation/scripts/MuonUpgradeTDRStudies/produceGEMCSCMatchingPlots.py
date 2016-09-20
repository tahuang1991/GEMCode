import os

from ROOT import *

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT 
gROOT.SetBatch(1)

from GEMCSCValidation import *

def enum(*sequential, **named):
  enums = dict(zip(sequential, range(len(sequential))), **named)
  reverse = dict((value, key) for key, value in enums.iteritems())
  enums['reverse_mapping'] = reverse
  return type('Enum', (), enums)

class MuonUpgradeTDREfficiencyPlotter():
  def __init__(self):
    self.inputDir = os.getenv("CMSSW_BASE") + "/src/GEMCode/GEMValidation/scripts/MuonUpgradeTDRStudies/"
    self.inputFile = "../../test/out_ana_efficiency.root"
    self.targetDir = "./"
    self.ext = ".png"
    self.analyzer = "MuonUpgradeTDREfficiency"
    self.effSt = "trk_eff_"
    self.stations = enum('ALL','ME11','ME1a','ME1b','ME12','ME13','ME21','ME22','ME31','ME32','ME41','ME42')
    self.stationsToUse = [self.stations.ME11,
                          self.stations.ME1a,
                          self.stations.ME1b,
                          self.stations.ME12,
                          self.stations.ME13,
                          self.stations.ME21,
                          self.stations.ME22,
                          self.stations.ME31,
                          self.stations.ME32,
                          self.stations.ME41,
                          self.stations.ME42]
    self.file = TFile.Open(self.inputDir + self.inputFile)
    self.dirAna = (self.file).Get(self.analyzer)
    self.treeEffSt = []
    for x in self.stationsToUse:
      self.treeEffSt.append(self.dirAna.Get(self.effSt + self.stations.reverse_mapping[x]))
    self.yMin = 0.
    self.yMax = 1.1
    self.etaMin = 1.15
    self.etaMax = 2.45
    self.pu = 0
    self.matchAlctGem = True

print "Make new plotter"
plotter = MuonUpgradeTDREfficiencyPlotter()
for i in range(len(plotter.stationsToUse)):
  st = plotter.stationsToUse[i]
  print "Processing station ", plotter.stations.reverse_mapping[st]
  continue
"""
simTrackToCscSimHitMatching(plotter,st)
simTrackToCscStripsWiresMatching(plotter,st)
simTrackToCscStripsWiresMatching_2(plotter,st)
simTrackToCscAlctClctMatching(plotter,st)
simTrackToCscAlctClctMatching_2(plotter,st)
simTrackToCscLctMatching(plotter,st)
simTrackToCscMpLctMatching(plotter,st)
"""

  
