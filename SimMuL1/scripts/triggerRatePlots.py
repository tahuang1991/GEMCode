## ROOT modules
from ROOT import *
from triggerPlotHelpers import getTree
import os

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT
ROOT.gROOT.SetBatch(1)


def draw_1D(target_dir, c_title, ext, t, title, h_name, h_bins, to_draw, cut = TCut(""), opt = ""):
  gStyle.SetStatStyle(0)
  gStyle.SetOptStat(1110)
  c = TCanvas("c","c",600,600)
  c.Clear()
  t.Draw(to_draw + ">>" + h_name + h_bins, cut)
  print t
  h = TH1F(gDirectory.Get(h_name).Clone(h_name))
  if not h:
    sys.exit('h does not exist')
  h.SetTitle(title)
  h.SetLineWidth(2)
  h.SetLineColor(kBlue)
  h.Draw(opt)
  h.SetMinimum(0.)
  c.SaveAs(target_dir + c_title + ext)

#_______________________________________________________________________________
if __name__ == "__main__":

  input_dir = "files/"
  output_dir = "plots_trigger/"
  ext = ".png"
  pwd = os.getenv("CMSSW_BASE")
  fileName = "hp_minbias_CMSSW_6_2_0_SLHC10_upgradePLS3_pu0_w3rate.root"
  path = os.path.join(pwd,"src",fileName)
  print path

  tree = getTree(fileName,"GEMCSCTriggerRateTree","GMTCand")
  draw_1D(output_dir, "rate_test", ext, tree, "pt;pt;pt", "h_", "(40,0,40)", "pt")
  
  

  """
  The simplest can be to take
  GMTCand tree
  and just plot the pt for the GE1/1 eta range for quality >=3 (double mu)
  and quality >=4 (single mu)
  for various SLHC and GEM-CSC ILT configurations.
  """
