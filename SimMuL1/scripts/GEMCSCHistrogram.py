
## custom modules
from effFunctions import *
from cuts import *
from tdrStyle import *
from triggerPlotHelpers import *
from GEMCSCdPhiDict import *
from math import *
from mkdir import mkdir

## ROOT modules
from ROOT import *
import os.path
## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT
ROOT.gROOT.SetBatch(1)


gStyle.SetStatW(0.07)
gStyle.SetStatH(0.06)

gStyle.SetOptStat(0)

gStyle.SetTitleStyle(0)
gStyle.SetTitleAlign(13) ## coord in top left
gStyle.SetTitleX(0.)
gStyle.SetTitleY(1.)
gStyle.SetTitleW(1)
gStyle.SetTitleH(0.058)
gStyle.SetTitleBorderSize(0)

gStyle.SetPadLeftMargin(0.126)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadBottomMargin(0.13)

gStyle.SetMarkerStyle(1)

global input_file
global input_file1
global input_file2
global input_file3
global input_file4
global input_file5
global input_file6

#_______________________________________________________________________________
def getTree(fileName):
    """Get tree for given filename"""

    analyzer = "GEMCSCAnalyzer"
    trk_eff = "trk_eff_ME11"

    file = TFile.Open(fileName)
    if not file:
        sys.exit('Input ROOT file %s is missing.' %(fileName))

    dir = file.Get(analyzer)
    if not dir:
        sys.exit('Directory %s does not exist.' %(dir))
        
    tree = dir.Get(trk_eff)
    if not tree:
        sys.exit('Tree %s does not exist.' %(tree))

    return tree


#_______________________________________________________________________________
def getDphi(eff,pt,evenOdd):
    """Return the delta phi cut value given: (1) an efficiency, (2) a pt value and (3) choice for even/odd chambers"""

    return dphi_lct_pad["%s"%(eff)]["%s"%(pt)]["%s"%(evenOdd)]

#_______________________________________________________________________________
def draw_hist(t, title, h_bins, to_draw, den_cut, opt = "", color = kBlue, marker_st = 1, marker_sz = 1.):
    """Make an efficiency plot"""
    
    ## total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    debug = False
    if debug:
        print "Denominator cut", den_cut
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    nBins  = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])
    
    den = TH1F("den", "", nBins, minBin, maxBin)

    t.Draw(to_draw + ">>den", den_cut, "goff")

        
    den.SetLineWidth(2)
    den.SetLineColor(color)
    den.Draw(opt + " same")
    den.SetMarkerStyle(marker_st)
    den.SetMarkerColor(color)
    den.SetMarkerSize(marker_sz)
    
    return den
    

#____________________________________________________________________
def simTrackwithLCT(filesDir, plotDir, ext):

    gStyle.SetTitleStyle(0);
    gStyle.SetTitleAlign(13); ##coord in top left
    gStyle.SetTitleX(0.);
    gStyle.SetTitleY(1.);
    gStyle.SetTitleW(1);
    gStyle.SetTitleH(0.058);
    gStyle.SetTitleBorderSize(0);
    
    gStyle.SetPadLeftMargin(0.126);
    gStyle.SetPadRightMargin(0.04);
    gStyle.SetPadTopMargin(0.06);
    gStyle.SetPadBottomMargin(0.13);
    gStyle.SetOptStat(0);
    gStyle.SetMarkerStyle(1);
    

    etareb = 1
    yrange = [0.8,1.005]
    xrange = [1.4,2.5]    

    t1 = getTree("%s%s"%(filesDir, input_file))
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 11 + "LCT from GEMCopad and ALCT on any CSC chamber in ME11 " + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "Generated muon #eta"
    yTitle = "number"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot = "TMath::Abs(eta)"
    h_bins = "(40,1.5001,2.5001)"
#    h_bins = "(100,-2.5,2.5)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,600)
    c.Clear()
    base  = TH1F("base",title,nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
      
     
    Cut_den = AND(OR(ok_sh2,ok_sh1),OR(TCut("dphi_lct_odd==0"), TCut("dphi_lct_even==0")))
    h1 = draw_hist(t1, title, h_bins, toPlot, Cut_den,"same", kRed,21)
#    h2 = draw_geff(t2, title, h_bins, toPlot, Cut_den, Cut_num, "same", kBlue,2)
    h1.SetTitle(topTitle)
    h1.SetXTitle(xTitle)
    h1.SetYTitle(yTitle)
    print "the Entries: ",h1.GetEntries()
    h1.Draw()
     
    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
#    leg.AddEntry(h_eta_initial,"CSC simhits 3/6, by CSC Ana","pl")
#    leg.AddEntry(h_eta_initial,"CSC+GEM 4/8+ALCT-Copad, by CSC Ana","pl")
#    leg.AddEntry(h1,"CSC simhits 3/6, by GEM Ana","pl")
    leg.AddEntry(h1,"simtrack with LCT and GEMDPhi=0","pl")
    leg.Draw("same")
   
    
    tex = TLatex(.25,.5,"PU0")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")

    c.Print("%sPU0_99LCTs_GEMDPhi0%s"%(plotDir,ext))

#____________________________________________________________________
def simTrackwithLCTVsGEMDPhi(filesDir, plotDir, ext):

    gStyle.SetTitleStyle(0);
    gStyle.SetTitleAlign(13); ##coord in top left
    gStyle.SetTitleX(0.);
    gStyle.SetTitleY(1.);
    gStyle.SetTitleW(1);
    gStyle.SetTitleH(0.058);
    gStyle.SetTitleBorderSize(0);
    
    gStyle.SetPadLeftMargin(0.126);
    gStyle.SetPadRightMargin(0.04);
    gStyle.SetPadTopMargin(0.06);
    gStyle.SetPadBottomMargin(0.13);
    gStyle.SetOptStat(0);
    gStyle.SetMarkerStyle(1);
    

    etareb = 1
    yrange = [0.8,1.005]
    xrange = [1.4,2.5]    

    t1 = getTree("%s%s"%(filesDir, input_file))
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 11 + "Simtrack with LCT  in ME11" + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "Generated muon #GEMDPhi"
    yTitle = "number"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "dphi_lct_odd"
    toPlot2 = "dphi_lct_even"

    h_bins = "(50,-0.0250,0.0250)"
#    h_bins = "(100,-2.5,2.5)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,600)
    c.Clear()
    base  = TH1F("base",title,nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
      
     
    Cut_den1 = AND(ok_sh1,ok_lct1)
    Cut_den3 = AND(ok_sh1,ok_lct1,TCut("pt>10"))
    Cut_den2 = AND(ok_sh2,ok_lct2)
    Cut_den4 = AND(ok_sh2,ok_lct2,TCut("pt>10"))
    h1 = draw_hist(t1, title, h_bins, toPlot1, Cut_den1,"same", kRed,21)
    h3 = draw_hist(t1, title, h_bins, toPlot1, Cut_den3,"same", kMagenta,20)
    h2 = draw_hist(t1, title, h_bins, toPlot2, Cut_den2,"same", kBlue,2)
    h4 = draw_hist(t1, title, h_bins, toPlot2, Cut_den4,"same", kGreen,3)
    h1.SetTitle(topTitle)
    h1.SetXTitle(xTitle)
    h1.SetYTitle(yTitle)
    print "the Odd Entries: ",h1.GetEntries()
    print "the Even Entries: ",h2.GetEntries()
    h1.Draw()
    h2.Draw("same")
#    h3.Draw("same")
#    h4.Draw("same")
     
    leg = TLegend(0.40,0.24,.87,0.5, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.SetHeader("Flipped hs in LCT to get GEMDPhi")
#    leg.AddEntry(h_eta_initial,"CSC simhits 3/6, by CSC Ana","pl")
#    leg.AddEntry(h_eta_initial,"CSC+GEM 4/8+ALCT-Copad, by CSC Ana","pl")
    leg.AddEntry(h2,"simtrack with LCT for even chamber","pl")
#    leg.AddEntry(h4,"simtrack with LCT for even chamber with pt>10","pl")
    leg.AddEntry(h1,"simtrack with LCT for odd chamber","pl")
#    leg.AddEntry(h3,"simtrack with LCT for odd chamber with pt>10","pl")
    leg.Draw("same")
   
    
    tex = TLatex(.25,.5,"PU0")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")

    c.Print("%sPU0_LCTs_VS_GEMDPhi_flipped%s"%(plotDir,ext))



#____________________________________________________________________
def simTrackwithLCTVsQual(filesDir, plotDir, ext):

    gStyle.SetTitleStyle(0);
    gStyle.SetTitleAlign(13); ##coord in top left
    gStyle.SetTitleX(0.);
    gStyle.SetTitleY(1.);
    gStyle.SetTitleW(1);
    gStyle.SetTitleH(0.058);
    gStyle.SetTitleBorderSize(0);
    
    gStyle.SetPadLeftMargin(0.126);
    gStyle.SetPadRightMargin(0.04);
    gStyle.SetPadTopMargin(0.06);
    gStyle.SetPadBottomMargin(0.13);
    gStyle.SetOptStat(0);
    gStyle.SetMarkerStyle(1);
    

    etareb = 1
    yrange = [0.8,1.005]
    xrange = [1.4,2.5]    

    t1 = getTree("%s%s"%(filesDir, input_file))
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 11 + "Simtrack with LCT  in ME11 " + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "Generated muon Quality"
    yTitle = "number"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "TMath::Abs(quality_odd)"
    toPlot2 = "TMath::Abs(quality_even)"

    h_bins = "(18,0,18)"
#    h_bins = "(100,-2.5,2.5)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,600)
    c.Clear()
    base  = TH1F("base",title,nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
      
     
    Cut_den1 = AND(ok_sh1,ok_lct1)
    Cut_den3 = AND(ok_sh1,ok_lct1,TCut("pt>10"))
    Cut_den2 = AND(ok_sh2,ok_lct2)
    Cut_den4 = AND(ok_sh2,ok_lct2,TCut("pt>10"))
    h1 = draw_hist(t1, title, h_bins, toPlot1, Cut_den1,"same", kRed,21)
    h3 = draw_hist(t1, title, h_bins, toPlot1, Cut_den3,"same", kMagenta,20)
    h2 = draw_hist(t1, title, h_bins, toPlot2, Cut_den2,"same", kBlue,2)
    h4 = draw_hist(t1, title, h_bins, toPlot2, Cut_den4,"same", kGreen,3)
    h1.SetTitle(topTitle)
    h1.SetXTitle(xTitle)
    h1.SetYTitle(yTitle)
    print "the Odd Entries: ",h1.GetEntries()
    print "the Even Entries: ",h2.GetEntries()
    h1.Draw()
    h2.Draw("same")
    h3.Draw("same")
    h4.Draw("same")
     
    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
#    leg.AddEntry(h_eta_initial,"CSC simhits 3/6, by CSC Ana","pl")
#    leg.AddEntry(h_eta_initial,"CSC+GEM 4/8+ALCT-Copad, by CSC Ana","pl")
    leg.AddEntry(h2,"simtrack with LCT for even chamber","pl")
    leg.AddEntry(h4,"simtrack with LCT for even chamber with pt>10","pl")
    leg.AddEntry(h1,"simtrack with LCT for odd chamber","pl")
    leg.AddEntry(h3,"simtrack with LCT for odd chamber with pt>10","pl")
    leg.Draw("same")
   
    
    tex = TLatex(.25,.5,"PU0")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")

    c.Print("%sPU0_LCTs_VS_Qaul%s"%(plotDir,ext))



#________________________________________________________________________________
if __name__ == "__main__":
    

    input_dir = "files/"
    output_dir = "GEMCSC_Hist/"
    input_file = "gem-csc_stub_ana.root"
    if not os.path.exists(output_dir):
	os.makedirs(output_dir)
     
    ext = ".eps"
    simTrackwithLCT(input_dir, output_dir, ext)
    simTrackwithLCTVsGEMDPhi(input_dir, output_dir, ext)
    simTrackwithLCTVsQual(input_dir, output_dir, ext)
