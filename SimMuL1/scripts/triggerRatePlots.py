## ROOT modules
from ROOT import *
#from triggerPlotHelpers import getTree
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
  return h

#_______________________________________________________________________________
if __name__ == "__main__":

  input_dir = "files/"
  output_dir = "plots_trigger/"
  ext = ".png"
  pwd = os.getenv("CMSSW_BASE")
  fileName = "hp_minbias_CMSSW_6_2_0_SLHC11_upgradePLS3_pu0_w3_rate.root"
  path = os.path.join(pwd,"src",fileName)
  print path

  ##  tree = getTree(fileName,"GEMCSCTriggerRateTree","GMTCand")  
  file = TFile.Open(fileName)
  dirAna = file.Get("GEMCSCTriggerRateTree")
  tree = dirAna.Get("GMTCand")
  
##  draw_1D(output_dir, "rate_gmtreg_test", ext, dirAna.Get("GMTRegCand"), "pt;pt;pt", "h_", "(40,0,40)", "pt")

  ## base histogram
  h0 = draw_1D(output_dir, "rate_gmt_test", ext, dirAna.Get("GMTCand"), ";pt;events", "h_", "(100,0,100)", "pt")
  h1 = draw_1D(output_dir, "rate_gmt_test", ext, dirAna.Get("GMTCand"), ";pt;events", "h_", "(100,0,100)", "pt","pt>10")
  
  c = TCanvas("c","c",600,600)
  c.Clear()
  h0.Draw("")
  c.SaveAs("h0.png")

  c = TCanvas("c","c",600,600)
  c.Clear()
  h1.Draw("")
  c.SaveAs("h1.png")

  """
  ## produce integrated trigger plots
  ptbins = [0,5.01,6.01,10.01,15.01,20.01,30.01,40.01]
  h = []
  h.extend([h0])
  for i in range(len(ptbins)):
    hh = h0
    min = hh.FindBin(ptbins[i])
    ## get the min and max
    if i!=len(ptbins)-1:
      max = hh.FindBin(ptbins[i+1])
    else:
      max = hh.GetNbinsX()
    ## partial integration
    for b in xrange(min,max):
      h0.SetBinContent(b, hh.GetBinContent(b));
      h0.SetBinError(b, hh.GetBinError(b));
      
  c = TCanvas("c","c",600,600)
  c.Clear()
  h0.Draw("")
  c.SaveAs("temp.png")
  """

  """
  The simplest can be to take
  GMTCand tree
  and just plot the pt for the GE1/1 eta range for quality >=3 (double mu)
  and quality >=4 (single mu)
  for various SLHC and GEM-CSC ILT configurations.
  """
