
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
gStyle.SetPadRightMargin(0.10)
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
def getTree(fileName, trk_eff = "trk_eff_ME1a"):
    """Get tree for given filename"""

    analyzer = "GEMCSCAnalyzer"
#    trk_eff = "trk_eff_ME1a"

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
    topTitle = " " * 11 + "LCT from GEMCopad and ALCT on any CSC chamber in ME1a " + " " * 35 + "CMS Simulation Preliminary"
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

    c.Print("%sPU0_ME1a_LCTs_GEMDPhi0%s"%(plotDir,ext))

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
    topTitle = " " * 11 + "Simtrack with LCT  in ME1a" + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "Generated muon #GEMDPhi"
    yTitle = "number"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "dphi_lct_odd"
    toPlot2 = "dphi_lct_even"

    h_bins = "(50,-0.025,0.025)"
#    h_bins = "(100,-2.5,2.5)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,600)
    c.Clear()
    base  = TH1F("base",title,nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
    base.SetMaximum(10000)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
      
     
    Cut_den1 = AND(ok_sh1,ok_lct1)
    Cut_den3 = AND(ok_sh1,ok_lct1,TCut("pt>10 && abs(eta)<2.1"))
    Cut_den2 = AND(ok_sh2,ok_lct2)
    Cut_den4 = AND(ok_sh2,ok_lct2,TCut("pt>10 && abs(eta)<2.1"))
    h1 = draw_hist(t1, title, h_bins, toPlot1, Cut_den1,"same", kRed,21)
    h3 = draw_hist(t1, title, h_bins, toPlot1, Cut_den3,"same", kMagenta,20)
    h2 = draw_hist(t1, title, h_bins, toPlot2, Cut_den2,"same", kBlue,2)
    h4 = draw_hist(t1, title, h_bins, toPlot2, Cut_den4,"same", kGreen,3)
    h1.SetTitle(topTitle)
    h1.SetXTitle(xTitle)
    h1.SetYTitle(yTitle)
    h3.SetTitle(topTitle)
    h3.SetXTitle(xTitle)
    h3.SetYTitle(yTitle)
    print "the Odd Entries: ",h1.GetEntries()
    print "the Even Entries: ",h2.GetEntries()

    h1.SetBinContent(1, h1.GetBinContent(0) + h1.GetBinContent(1))
    h1.SetBinContent(50, h1.GetBinContent(50) + h1.GetBinContent(50+1))

    h2.SetBinContent(1, h2.GetBinContent(0) + h2.GetBinContent(1))
    h2.SetBinContent(50, h2.GetBinContent(50) + h2.GetBinContent(50+1))

    h3.SetBinContent(1, h3.GetBinContent(0) + h3.GetBinContent(1))
    h3.SetBinContent(50, h3.GetBinContent(50) + h3.GetBinContent(50+1))

    h4.SetBinContent(1, h4.GetBinContent(0) + h4.GetBinContent(1))
    h4.SetBinContent(50, h4.GetBinContent(50) + h4.GetBinContent(50+1))

    base.Draw()
#    h2.Draw("same")
    h3.Draw("same")
    h4.Draw("same")
     
    leg = TLegend(0.40,0.24,.87,0.5, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.SetHeader("PU0, Abs(eta)<2.1")
#    leg.AddEntry(h_eta_initial,"CSC simhits 3/6, by CSC Ana","pl")
#    leg.AddEntry(h_eta_initial,"CSC+GEM 4/8+ALCT-Copad, by CSC Ana","pl")
#    leg.AddEntry(h2,"simtrack with LCT for even chamber","pl")
    leg.AddEntry(h4,"simtrack with LCT for even chamber with pt>10","pl")
#    leg.AddEntry(h1,"simtrack with LCT for odd chamber","pl")
    leg.AddEntry(h3,"simtrack with LCT for odd chamber with pt>10","pl")
    leg.Draw("same")
   
    
    tex = TLatex(.25,.5,"PU0")
    tex.SetTextSize(0.05)
    tex.SetNDC()
#tex.Draw("same")

    c.Print("%sPU0_ME1a_LCTs_VS_GEMDPhi_jason_loweta%s"%(plotDir,ext))


#____________________________________________________________________
def simTrackwithLCTHsVsGEMDPhi(filesDir, plotDir, ext):

    gStyle.SetTitleStyle(0);
    gStyle.SetTitleAlign(13); ##coord in top left
    gStyle.SetTitleX(0.);
    gStyle.SetTitleY(1.);
    gStyle.SetTitleW(1);
    gStyle.SetTitleH(0.058);
    gStyle.SetTitleBorderSize(0);
    
    gStyle.SetPadLeftMargin(0.126);
    gStyle.SetPadRightMargin(0.14);
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
    topTitle = " " * 11 + "Simtrack with LCT  in ME1a" + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "halfstrip"
    yTitle = "GEMDPhi"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "dphi_lct_odd:hs_lct_odd >> base"
    toPlot2 = "dphi_lct_even:hs_lct_even >> base2"

    y_bins = "(100,-0.25,0.25)"
    x_bins = "(160,0,160)"
#    h_bins = "(100,-2.5,2.5)"
#    h_bins = "(100,-2.5,2.5)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])

    c1 = TCanvas("c1","c1",800,600)
    c1.Clear()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    base  = TH2F("base",title,xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    base2  = TH2F("base2",title,xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
#    base.SetMinimum(0.0)
    base.SetMaximum(150)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.GetXaxis().SetTitle(xTitle)
    base.GetYaxis().SetTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    base2.SetMaximum(150)
    base2.GetXaxis().SetTitle(xTitle)
    base2.GetYaxis().SetTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
      
    
    cut1 = TCut("pt>10 & (has_lct&1)>0 & (has_gem_pad&1)>0")
    cut2 = TCut("pt>10 & (has_lct&2)>0 & (has_gem_pad&2)>0")
    t1.Draw(toPlot1, cut1, "COLZ")
    base.Draw("colz")
    gPad.SetLogz()
    tex1 = TLatex(0.2,0.3,"PU140,Pt>10, odd chamber")
    tex1.SetTextSize(0.05)
#tex1.SetTextFont(62)
    tex1.SetNDC()
    tex1.Draw("same")
    c1.Print("%sPU140_ME1a_Hs_VS_GEMDPhi_odd_1%s"%(plotDir,ext))

    c2 = TCanvas("c2","c2",800,600)
    c2.Clear()
    c2.SetGridx()
    c2.SetGridy()
    c2.SetTickx()
    c2.SetTicky()
    t1.Draw(toPlot2, cut2, "COLZ")
    base2.Draw("colz")
    gPad.SetLogz()
    tex2 = TLatex(0.2,0.3,"PU140,Pt>10, even chamber")
    tex2.SetTextSize(0.05)
#   tex2.SetTextFont(62)
    tex2.SetNDC()
    tex2.Draw("same")
     

    c2.Print("%sPU140_ME1a_Hs_VS_GEMDPhi_even_1%s"%(plotDir,ext))

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
    topTitle = " " * 11 + "Simtrack with LCT  in ME1a " + " " * 35 + "CMS Simulation Preliminary"
    xTitle = "Generated muon Quality"
    yTitle = "number"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "quality_odd"
    toPlot2 = "quality_even"

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
    base.GetXaxis().SetTitle(xTitle)
    base.GetYaxis().SetTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    base2.GetXaxis().SetTitle(xTitle)
    base2.GetYaxis().SetTitle(yTitle)
      
     
    Cut_den1 = AND(ok_sh1,ok_lct1, TCut("pt>20"))
    Cut_den3 = AND(ok_sh1,ok_lct1,TCut("pt>10 && hs_lct_odd >5 && hs_lct_odd <123"))
    Cut_den2 = AND(ok_sh2,ok_lct2, TCut("pt>20"))
    Cut_den4 = AND(ok_sh2,ok_lct2,TCut("pt>10 && hs_lct_even >5 && hs_lct_even <123"))
    h1 = draw_hist(t1, title, h_bins, toPlot1, Cut_den1,"same", kRed,21)
    h3 = draw_hist(t1, title, h_bins, toPlot1, Cut_den3,"same", kMagenta,20)
    h2 = draw_hist(t1, title, h_bins, toPlot2, Cut_den2,"same", kBlue,2)
    h4 = draw_hist(t1, title, h_bins, toPlot2, Cut_den4,"same", kGreen,3)
    h1.SetTitle(topTitle)
    h1.SetXTitle(xTitle)
    h1.SetYTitle(yTitle)
    print "the Odd Entries: ",h1.GetEntries()
    print "the Even Entries: ",h2.GetEntries()

    h1.SetBinContent(1, h1.GetBinContent(0) + h1.GetBinContent(1))
    h1.SetBinContent(18, h1.GetBinContent(18) + h1.GetBinContent(18+1))

    h2.SetBinContent(1, h2.GetBinContent(0) + h2.GetBinContent(1))
    h2.SetBinContent(18, h2.GetBinContent(18) + h2.GetBinContent(18+1))

    h3.SetBinContent(1, h3.GetBinContent(0) + h3.GetBinContent(1))
    h3.SetBinContent(18, h3.GetBinContent(18) + h3.GetBinContent(18+1))

    h4.SetBinContent(1, h4.GetBinContent(0) + h4.GetBinContent(1))
    h4.SetBinContent(18, h4.GetBinContent(18) + h4.GetBinContent(18+1))

    h1.Draw()
    h2.Draw("same")
#    h3.Draw("same")
#    h4.Draw("same")
     
    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
#    leg.AddEntry(h_eta_initial,"CSC simhits 3/6, by CSC Ana","pl")
#    leg.AddEntry(h_eta_initial,"CSC+GEM 4/8+ALCT-Copad, by CSC Ana","pl")
    leg.AddEntry(h2,"simtrack with LCT for even chamber, with pt>20","pl")
#    leg.AddEntry(h4,"simtrack with LCT for even chamber with pt>10","pl")
    leg.AddEntry(h1,"simtrack with LCT for odd chamber, with pt>20","pl")
#    leg.AddEntry(h3,"simtrack with LCT for odd chamber with pt>10","pl")
    leg.Draw("same")
   
    
    tex = TLatex(.25,.5,"PU0")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")

    c.Print("%sPU0_ME1a_LCTs_VS_Qual%s"%(plotDir,ext))


#____________________________________________________________________
def simTrackwithLCTHsVsHsfromGEM(filesDir, plotDir, ext):

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
    topTitle = " " * 8 + "comparsion between hs in LCT and hs extrapolated from GEM in ME1a" + " " * 10 + "CMS Simulation Preliminary"
    xTitle = "halfstrip"
    yTitle = "extrapolated halfstrip"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)
    toPlot1 = "hsfromgem_odd:hs_lct_odd"
    toPlot2 = "hsfromgem_even:hs_lct_even"

    h_bins = "(160,0,160)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,600)
    c.Clear()
    gPad.SetTitle(topTitle)

    base  = TH2D("base"," ",nBins,minBin,maxBin, nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.SetXTitle(xTitle)
    base.SetYTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    toPlot1 = "hsfromgem_odd:hs_lct_odd >> base"
    toPlot2 = "hsfromgem_even:hs_lct_even >> base"
      
    
    Pt_cut = TCut("pt>10")
    c.Divide(2,1)
    c.cd(1)
    t1.Draw(toPlot1, Pt_cut, "COLZ")
    gPad.SetLogz()
    base.Draw("COLZ")
    leg1 = TLegend(0.40,0.24,.87,0.5, "", "NDC");
    leg1.SetHeader("Pt>10, Odd chamber")
    leg1.Draw("same")

    c.cd(2)
    t1.Draw(toPlot2, Pt_cut, "COLZ")
    gPad.SetLogz()
    base.Draw("COLZ")
    leg2 = TLegend(0.40,0.24,.87,0.5, "", "NDC");
    leg2.SetHeader("Pt>10, Even chamber")
    leg2.Draw("same")

     
    leg = TLegend(0.40,0.24,.87,0.5, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)

    c.Print("%sPU0_ME1a_Hs_VS_HsfromGEM%s"%(plotDir,ext))

#____________________________________________________________________
def simTrackwithLCTHsVsHsfromRPC(filesDir, plotDir, ext):

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
    topTitle = " " * 8 + "comparsion between hs in LCT and hs extraploted from RPC in ME1a" + " " * 5 + "CMS Simulation Preliminary"
    xTitle = "halfstrip"
    yTitle = "extraploted hs"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)

    h_bins = "(160,0,160)"
#    h_bins = "(100,-2.5,2.5)"
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])

    c = TCanvas("c","c",800,400)
    c.Clear()
    base  = TH2F("base"," ",nBins,minBin,maxBin,nBins,minBin,maxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.SetXTitle(xTitle)
    base.SetYTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    toPlot1 = "hsfromrpc_odd:hs_lct_odd >> base"
    toPlot2 = "hsfromrpc_even:hs_lct_even >> base"
      
    
    Pt_cut = TCut("pt>30 && eta<0")
    c.Divide(2,1)

    c.cd(1)
    t1.Draw(toPlot1, Pt_cut, "COLZ")
    gPad.SetLogz()
    base.Draw("COLZ")
    leg1 = TLegend(0.40,0.20,.80,0.30, "", "NDC");
    leg1.SetFillColor(ROOT.kWhite)
#leg1.SetHeader("Pt>10, Odd chamber")
    leg1.AddEntry(base,"Pt>30, Odd chamber","pl")
#leg1.Draw("same")
    tex1 = TLatex(.25,.5,"PU0, Pt>30, odd chamber")
    tex1.SetTextSize(0.05)
    tex1.SetNDC()
    tex1.Draw("same")

    c.cd(2)
    t1.Draw(toPlot2, Pt_cut, "COLZ")
    gPad.SetLogz()
    base.Draw("COLZ")
    leg2 = TLegend(0.40,0.20,.80,0.30, "", "NDC");
    leg2.SetFillColor(ROOT.kWhite)
#leg2.SetHeader("Pt>10, Even chamber")
    leg2.AddEntry(base,"Pt>30, Even chamber","pl")
#    leg2.Draw("same")
    tex2 = TLatex(.25,.5,"PU0, Pt>30, even chamber")
    tex2.SetTextSize(0.05)
    tex2.SetNDC()
    tex2.Draw("same")
     

    c.SaveAs("%sPU0_ME1a_Hs_VS_HsfromRPC_30Gev%s"%(plotDir,ext))
    c.SaveAs("%sPU0_ME1a_Hs_VS_HsfromRPC_30Gev%s"%(plotDir,".pdf"))

#____________________________________________________________________
def simTrackEtaVsPhi(filesDir, plotDir, ext):

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

    t1 = getTree("%s%s"%(filesDir, input_file),"trk_eff_ME1a")
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 8 + " " * 30 + "CMS Simulation Preliminary"
    xTitle = "#phi"
    yTitle = "#eta"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)

    x_bins = "(60,-3.14,3.14)"
    y_bins = "(100,1.5,2.5)"
    xnBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])

    ynBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    c = TCanvas("c","c",800,400)
   # c.Clear()
    base  = TH2F("base"," ",xnBins,xminBin,xmaxBin,ynBins,yminBin,ymaxBin)
    base2  = TH2F("base2"," ",xnBins,xminBin,xmaxBin,ynBins,yminBin,ymaxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.SetXTitle(xTitle)
    base.SetYTitle(yTitle)
    base2.SetXTitle(xTitle)
    base2.SetYTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    toPlot1 = "abs(eta):phi >> base"
    toPlot2 = "abs(eta):phi >> base2"
      
    
    cut1 = TCut("pt>10 && (dphi_lct_odd)==-99")
    cut2 = TCut("pt>10 && (dphi_lct_even)==-99")
#c.Divide(2,1)
    c.cd()
    pad1 = TPad("pad1"," ",0.0,0.0,0.5,0.9)
    pad2 = TPad("pad2"," ",0.5,0.0,1.0,0.9)
    pad1.Draw()
    pad2.Draw("same")
    tex3 = TLatex(.35,.90,"Eta Vs Phi in ME1a for (GEMDPhi)=-99")#title
    tex3.SetTextSize(0.05)
    tex3.SetNDC()
    tex3.Draw("same")


    pad1.cd()
    t1.Draw(toPlot1, cut1, "COLZ")
    pad1.SetLogz()
    base.Draw("COLZ")
    tex1 = TLatex(.25,.5,"PU140, Pt>10, odd chamber")
    tex1.SetTextSize(0.05)
    tex1.SetNDC()
    tex1.Draw("same")
    pad1.Update()

    pad2.cd()
    t1.Draw(toPlot2, cut2, "COLZ")
    pad2.SetLogz()
    base2.Draw("COLZ")
    tex2 = TLatex(.25,.5,"PU140, Pt>10, even chamber")
    tex2.SetTextSize(0.05)
    tex2.SetNDC()
    tex2.Draw("same")
    pad2.Update()
     
    c.Update()

    c.SaveAs("%sPU140_ME1a_Eta_Vs_Phi2%s"%(plotDir,ext))
    c.SaveAs("%sPU140_ME1a_Eta_Vs_Phi2%s"%(plotDir,".pdf"))
# c.SaveAs("%sPU140_ME1a_Eta_Vs_Phi2%s"%(plotDir,".C"))


#____________________________________________________________________
def simTrackPhiComparison(filesDir, plotDir, xaxis, yaxis, x_bins, y_bins, ext):

    gStyle.SetTitleStyle(0);
    gStyle.SetTitleAlign(13); ##coord in top left
    gStyle.SetTitleX(0.);
    gStyle.SetTitleY(1.);
    gStyle.SetTitleW(1);
    gStyle.SetTitleH(0.058);
    gStyle.SetTitleBorderSize(0);
    
    gStyle.SetPadLeftMargin(0.126);
    gStyle.SetPadRightMargin(0.10);
    gStyle.SetPadTopMargin(0.06);
    gStyle.SetPadBottomMargin(0.13);
    gStyle.SetOptStat(0);
    gStyle.SetMarkerStyle(1);
    

    etareb = 1
    yrange = [0.8,1.005]
    xrange = [1.4,2.5]    

    t1 = getTree("%s%s"%(filesDir, input_file),"trk_eff_ME1a")
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 8 + " " * 30 + "CMS Simulation Preliminary"
    xTitle = xaxis
    yTitle = yaxis
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)

#    x_bins = "(60,-3.14,3.14)"
#    y_bins = "(100,1.5,2.5)"
    xnBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])

    ynBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    c = TCanvas("c","c",600,400)
   # c.Clear()
    base  = TH2F("base"," ",xnBins,xminBin,xmaxBin,ynBins,yminBin,ymaxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base.SetXTitle(xTitle)
    base.SetYTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    toPlot1 = "%s:%s >> base"%(yaxis,xaxis)
      
    
#    cut1 = TCut("pt>10 && (has_csc_sh&1)>0 && (has_lct&1)>0")
    cut1 = TCut("pt>10 && (has_csc_sh&1)>0 && (has_lct&1)>0 && (has_gem_pad&1)>0")
#c.Divide(2,1)
    c.cd()
    t1.Draw(toPlot1, cut1, "COLZ")
    gPad.SetLogz()
    base.Draw("COLZ")
    tex1 = TLatex(.25,.3,"PU140, Pt>10,odd chamber")
    tex1.SetTextSize(0.05)
    tex1.SetNDC()
    tex1.Draw("same")
     
    c.Update()

    c.SaveAs("%sPU140_ME1a_%s_Vs_dphi%s"%(plotDir,xaxis,ext))
    c.SaveAs("%sPU140_ME1a_%s_Vs_dphi%s"%(plotDir,xaxis,".pdf"))
# c.SaveAs("%sPU140_ME1a_Eta_Vs_Phi2%s"%(plotDir,".C"))
    

#____________________________________________________________________
def simTrackGEMDPhiVsPt(filesDir, plotDir, ext):

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
#gStyle.SetOptStat(0111111);
    

    etareb = 1
    yrange = [0.8,1.005]
    xrange = [1.4,2.5]    

    t1 = getTree("%s%s"%(filesDir, input_file),"trk_eff_ME1a")
#    t2 = getTree("%s%s"%(filesDir, input_file4))

    ## variables for the plot
    topTitle = " " * 8 + " " * 30 + "CMS Simulation Preliminary"
    xTitle = "#mu P_{T} [GeV/c]"
    yTitle = "GEMDPhi"
    title = "%s;%s;%s"%(topTitle,xTitle,yTitle)

    x_bins = "(25,0,50)"
    y_bins = "(100,-0.05,0.05)"
    xnBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])

    ynBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    c = TCanvas("c","c",800,400)
   # c.Clear()
    base  = TH2F("base"," ",xnBins,xminBin,xmaxBin,ynBins,yminBin,ymaxBin)
    base2  = TH2F("base2"," ",xnBins,xminBin,xmaxBin,ynBins,yminBin,ymaxBin)
#    base.SetMinimum(0.0)
#    base.SetMaximum(1.02)
#    base.Draw("")
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetLabelSize(0.05)
    base2.GetXaxis().SetLabelSize(0.05)
    base2.GetYaxis().SetLabelSize(0.05)
    base.SetXTitle(xTitle)
    base.SetYTitle(yTitle)
    base2.SetXTitle(xTitle)
    base2.SetYTitle(yTitle)
#    base.GetYaxis().SetRangeUser(yrange[0],yrange[1])
    toPlot1 = "(phi_lct_odd-phi_pad_odd):pt >> base"#y:x
    toPlot2 = "(phi_lct_even-phi_pad_even):pt >> base2"
      
    
    cut1 = TCut("(has_lct&1)>0 && (has_gem_pad&1)>0 && pt>10")
    cut2 = TCut("(has_lct&2)>0 && (has_gem_pad&2)>0 && pt>10")
#c.Divide(2,1)
    c.cd()
    pad1 = TPad("pad1"," ",0.0,0.0,0.5,0.9)
    pad2 = TPad("pad2"," ",0.5,0.0,1.0,0.9)
    pad1.Draw()
    pad2.Draw("same")
    tex3 = TLatex(.40,.90,"GEMDPhi Vs Pt in ME1a")#title
    tex3.SetTextSize(0.05)
    tex3.SetNDC()
    tex3.Draw("same")

    dphi1  = TH2F("dphi1","dphi cut Vs pt under 98% ",xnBins,xminBin,xmaxBin,50,0,0.005)
    dphi2  = TH2F("dphi2","dphi cut Vs pt under 98% ",xnBins,xminBin,xmaxBin,50,0,0.005)
    
    dphi1.SetMarkerColor(kRed)
    dphi2.SetMarkerColor(kBlue)

    dphi1.SetXTitle(xTitle)
    dphi1.SetYTitle("dphi cut")
    dphi1.SetMarkerStyle(21)
    dphi2.SetMarkerStyle(19)

    pad1.cd()
    t1.Draw(toPlot1, cut1, "COLZ")
    num_odd = 0
    base.SetMinimum(2.5)
#    for i in range(2,xnBins):
#	temp = base.Integral(i,i,0,ynBins)
#	for j in range(0,ynBins):
#	    temp_j = base.Integral(i,i,j,ynBins-j)
#	    ratio = temp_j/temp
#	    print "total", temp, "  j",temp_j," ratio",ratio
#	    if ratio < 0.98:
#	          dphi1.Fill(i*2, ymaxBin-j*0.0004)
#		  break
#	    print "i ",i," j ",j,"  ", base.GetBinContent(i,j)
#    if base.GetBinContent(i,j) < 1.5 :
#		num_odd = num_odd + base.GetBinContent(i,j)
#   print "num_odd ",num_odd,"  total number odd", base.GetEntries()
#    base.SetMinimum(3.5)
    pad1.SetLogz()
    pad1.SetRightMargin(0.14);
    base.Draw("COLZ")
    tex1 = TLatex(.25,.2,"PU140, odd chamber")
    tex1.SetTextSize(0.05)
    tex1.SetNDC()
    tex1.Draw("same")
    pad1.Update()

    pad2.cd()
    t1.Draw(toPlot2, cut2, "COLZ")
    num_even = 0
    base2.SetMinimum(2.5)
#   for i in range(0,xnBins):
#	for j in range(0,ynBins):
#print "i ",i," j ",j,"  ", base2.GetBinContent(i,j)
#	    if base2.GetBinContent(i,j) < 1.5 :
#		num_even = num_even + base2.GetBinContent(i,j)
#print "num_even ",num_even,"  total number even", base2.GetEntries()
#for i in range(2,xnBins):
#	temp = base2.Integral(i,i,0,ynBins)
#	for j in range(0,ynBins):
#	    temp_j = base2.Integral(i,i,j,ynBins-j)
#	    ratio = temp_j/temp
#	    print "total", temp, "  j",temp_j," ratio",ratio
#	    if ratio < 0.98:
#	          dphi2.Fill(i*2, ymaxBin-j*0.0004)
#		  break
#	    print "i ",i," j ",j,"  ", base.GetBinContent(i,j)
#    base2.SetMinimum(2.5)
    pad2.SetLogz()
    pad2.SetRightMargin(0.14);
    base2.Draw("COLZ")
    tex2 = TLatex(.25,.2,"PU140, even chamber")
    tex2.SetTextSize(0.05)
    tex2.SetNDC()
    tex2.Draw("same")
    pad2.Update()
     
    c.Update()

    c.SaveAs("%sPU140_ME1a_GEMDPhi_Vs_Pt_1%s"%(plotDir,ext))
    c.SaveAs("%sPU140_ME1a_GEMDPhi_Vs_Pt_1%s"%(plotDir,".pdf"))
#c.SaveAs("%sPU140_ME1a_new_GEMDPhi_Vs_Pt%s"%(plotDir,".C"))
"""
    c2 = TCanvas("c2","c2",800,600)
    dphi1.Draw("l")
    dphi2.Draw("lsame")
    leg1 = TLegend(0.40,0.20,.80,0.30, "", "NDC");
    leg1.SetFillColor(ROOT.kWhite)
#  leg1.SetHeader("Pt>10, Odd chamber")
    leg1.AddEntry(dphi1,"Odd chamber","pl")
    leg1.AddEntry(dphi2,"even chamber","pl")
    leg1.Draw("same")
    tex1 = TLatex(.25,.5,"PU0, Pt>30, odd chamber")
    tex1.SetTextSize(0.05)
    tex1.SetNDC()
    c2.SaveAs("%sPU140_ME1a_98_dphicut_Vs_Pt%s"%(plotDir,ext))
    c2.SaveAs("%sPU140_ME1a_98_dphicut_Vs_Pt%s"%(plotDir,".pdf"))

"""
#________________________________________________________________________________
if __name__ == "__main__":
    

    input_dir = "files/"
    output_dir = "GEMCSC_Hist/"
#    input_file = "GEMCSC_Ana_PU0_100k_All.root"
# input_file = "TestRPC_PU0_100k_2023_fixeven_GEMCSCAna.root"
# input_file = "Test_PU0_100k_2023_GEMCSCAna.root"
#input_file = "PU140_100k_2023_TEST_GEMCSC.root"
    input_file = "PU140_100k_2023_FixBX_GEMCSC.root"
    if not os.path.exists(output_dir):
	os.makedirs(output_dir)
     
    ext = ".png"
#   simTrackwithLCT(input_dir, output_dir, ext)
#   simTrackwithLCTVsGEMDPhi(input_dir, output_dir, ext)
#    simTrackwithLCTHsVsGEMDPhi(input_dir, output_dir, ext)
#    simTrackwithLCTHsVsHsfromGEM(input_dir, output_dir, ext)
    #simTrackwithLCTHsVsHsfromRPC(input_dir, output_dir, ext)
#simTrackEtaVsPhi(input_dir, output_dir, ext)
#    simTrackGEMDPhiVsPt(input_dir, output_dir, ext)
#simTrackwithLCTVsQual(input_dir, output_dir, ext)   
    xaxis = "hs_lct_odd"
#yaxis = "(phi_lct_even-phi_pad_even)"
    yaxis = "dphi_lct_odd"
    x_bins = "(60,-3.14,3.14)"
    y_bins = "(60,-3.14,3.14)"
    hs_bins = "(96,0,96)"
    hs_bins1 = "(128,0,128)"
    dphi_bins = "(100,-0.25,0.25)"
    simTrackPhiComparison(input_dir, output_dir, xaxis, yaxis, hs_bins, dphi_bins, ext)
