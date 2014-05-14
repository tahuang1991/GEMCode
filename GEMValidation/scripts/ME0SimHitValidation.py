from ROOT import *

from cuts import *
from drawPlots import *

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT 
ROOT.gROOT.SetBatch(1)

#_______________________________________________________________________________
def me0SimHitOccupancyXY(plotter,i):
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits,
           plotter.pre[i] + " SimHit occupancy: region-1, layer1;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(1),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l2" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1, layer2;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(2),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l3" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1, layer3;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(3),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l4" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1, layer4;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(4),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l5" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1, layer5;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(5),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rm1_l6" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1, layer6;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(-1),la(6),plotter.sel[i]), "COLZ")
  
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer1;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(1),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l2" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer2;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(2),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l3" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer3;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(3),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l4" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer4;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(4),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l5" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer5;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(5),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_xy_rp1_l6" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1, layer6;globalX [cm];globalY [cm]", 
           "h_", "(120,-280,280,120,-280,280)", "globalX:globalY", AND(re(1),la(6),plotter.sel[i]), "COLZ")

#_______________________________________________________________________________
def me0SimHitOccupancyRZ(plotter,i):
  draw_occ(plotter.targetDir, "sh_me0_zr_rm1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region-1;globalZ [cm];globalR [cm]", 
           "h_", "(80,-555,-515,120,20,280)", "sqrt(globalX*globalX+globalY*globalY):globalZ", AND(re(-1),plotter.sel[i]), "COLZ")
  draw_occ(plotter.targetDir, "sh_me0_zr_rp1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit occupancy: region1;globalZ [cm];globalR [cm]", 
           "h_", "(80,515,555,120,20,280)", "sqrt(globalX*globalX+globalY*globalY):globalZ", AND(re(1),plotter.sel[i]), "COLZ")

#_______________________________________________________________________________
def me0SimHitTOF(plotter,i):
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer1;Time of flight [ns];entries", 
          "h_", "(40,15,19)", "timeOfFlight", AND(re(-1),la(1),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l2" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer2;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(-1),la(2),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l3" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer3;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(-1),la(3),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l4" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer4;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(-1),la(4),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l5" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer5;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(-1),la(5),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rm1_l6" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region-1, layer6;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(-1),la(6),plotter.sel[i]))
  
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l1" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer1;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(1),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l2" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer2;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(2),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l3" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer3;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(3),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l4" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer4;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(4),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l5" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer5;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(5),plotter.sel[i]))
  draw_1D(plotter.targetDir, "sh_me0_tof_rp1_l6" + plotter.suff[i], plotter.ext, plotter.treeME0SimHits, plotter.pre[i] + " SimHit TOF: region1, layer6;Time of flight [ns];entries", 
          "h_", "(40,15,22)", "timeOfFlight", AND(re(1),la(6),plotter.sel[i]))
  
#_______________________________________________________________________________
def me0SimTrackToSimHitMatchingLX(plotter):
  pass

#_______________________________________________________________________________
def me0SimTrackToSimHitMatchingLY(plotter):
  pass

#_______________________________________________________________________________
def me0SimTrackToSimHitMatchingEta(plotter): 
  pass

#_______________________________________________________________________________
def me0SimTrackToSimHitMatchingPhi(plotter):
  pass
