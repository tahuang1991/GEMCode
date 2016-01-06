## ROOT modules
from ROOT import *
#from triggerPlotHelpers import getTree
import os

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT
ROOT.gROOT.SetBatch(1)

## return cumulative histogram
def getIntegral(hist):
  nBins = hist.GetXaxis().GetNbins()
  min = hist.GetXaxis().GetXmin()
  max = hist.GetXaxis().GetXmax()
  integral = TH1F("integral","", nBins, min, max)
  for i in range(0,nBins):
    total = 0
    for k in range(i,nBins):
      total = int(total + hist.GetBinContent(i))
    integral.SetBinContent(i,total);
  return integral

## scale to get rate
def applyScaling(hist):
  freq = 40000.
  bxWidth = 3
  ##need to check this
  bxFilling = 0.795
  hist.Scale(freq*bxFilling/nEvents/bxWidth);
  return hist

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

def draw_rate(target_dir, c_title, ext, t, title, h_name, h_bins, to_draw, cut = TCut(""), opt = ""):
  gStyle.SetStatStyle(0)
  gStyle.SetOptStat(1110)
  t.Draw(to_draw + ">>" + h_name + h_bins, cut)
  print t
  h = TH1F(gDirectory.Get(h_name).Clone(h_name))
  if not h:
    sys.exit('h does not exist')
  hist = getIntegral(hist)
  hist = applyScaling(hist)
  h.SetTitle(title)
  h.SetLineWidth(2)
  h.SetLineColor(kBlue)
  h.SetMinimum(0.)
  c = TCanvas("c","c",600,600)
  h.Draw(opt)
  c.Clear()
  c.SaveAs(target_dir + c_title + ext)
  return h


#_______________________________________________________________________________
if __name__ == "__main__":

  input_dir = "files/"
  output_dir = "plots_trigger/"
  ext = ".png"
  pwd = os.getenv("CMSSW_BASE")
  fileName = "hp_minbias_CMSSW_6_2_0_SLHC11_upgradePLS3_pu0_w3_rate.root"
  fileName = "hp_minbias_CMSSW_6_2_0_SLHC12_upgradePLS3_pu140_w3_rate.root"
  path = os.path.join(pwd,"src",fileName)
  print path

  ##  tree = getTree(fileName,"GEMCSCTriggerRateTree","GMTCand")  
  file = TFile.Open(fileName)
  dirAna = file.Get("GEMCSCTriggerRateTree")
  tree = dirAna.Get("GMTCand")

  nEvents =  dirAna.Get("h_events").GetEntries()
##  draw_1D(output_dir, "rate_gmtreg_test", ext, dirAna.Get("GMTRegCand"), "pt;pt;pt", "h_", "(40,0,40)", "pt")

  ## base histogram
  h0 = draw_1D(output_dir, "rate_gmt_test", ext, dirAna.Get("TFCand"), "TFTrack trigger rate;L1 candidate pt;Rate (unscaled)", "h_", "(100,0,100)", "pt")
  h1 = draw_1D(output_dir, "rate_gmt_test", ext, dirAna.Get("TFCand"), "TFTrack trigger rate;L1 candidate pt;Rate (unscaled)", "h_", "(100,0,100)", "pt","pt>10")
  
  c = TCanvas("c","c",600,600)
  c.Clear()
  gPad.SetLogy()
  h00 = getIntegral(h0)
#  h00 = applyScaling(h00)
  h00.GetXaxis().SetTitle("L1 Muon candidate p_{T} [GeV/c]")
  h00.GetYaxis().SetTitle("Rate [kHz]")
  h00.SetTitle("CSC Track Finder candidate")
  h00.Sumw2()
#  h00.Draw("l")
  h0.Draw("l")
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
