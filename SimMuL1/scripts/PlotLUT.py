from ROOT import *
from cuts import *

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

global input_file1 
global input_dir 
global output_dir 
global ext 

#_______________________________________________________________________________
def getTree(fileName,trk_eff = "trk_eff_ME11"):
    """Get tree for given filename"""

    analyzer = "GEMCSCAnalyzer"

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
#______________________________________________________________________________________
def Plothseta():
    input_file1 = "PU0_AllImprovement_ByGEM_ANA.root"
    input_dir = "files/"
    output_dir = "csc_gem_matching_Combined/"
    ext = ".png"

# total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    t = getTree("%s%s"%(input_dir, input_file1))
    L_cut = ok_lct1
    to_draw = "hs_lct_odd:eta_lct_odd"
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    s_bins = "(127,0,127)"
    w_bins = "(50,1.5,2.5)"
    p_bins = "(96,0,96)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    wirenBins  = int(w_bins[1:-1].split(',')[0])
    wmin = float(w_bins[1:-1].split(',')[1])
    wmax = float(w_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    LUT = TH2F("LUT","", stripnBins, smin, smax, wirenBins, wmin, wmax)
    
    C = TCanvas("C","",800,600)
# t.Draw(to_draw + ">>LUT", L_cut, "goff")
    t.Draw("hs_lct_odd:TMath::Abs(eta_lct_odd) ", L_cut, "COLZ")
#    LUT.Draw("same")

    C.Print("%shalfstrip_eta%s"%(output_dir,ext))
     
#______________________________________________________________________________________
def PlothsWire():


    input_file1 = "PU0_AllImprovement_ByGEM_ANA.root"
    input_dir = "files/"
    output_dir = "csc_gem_matching_Combined/"
    ext = ".png"

# total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    t = getTree("%s%s"%(input_dir, input_file1))
    L_cut = TCut("TMath::Abs(eta)>2.1 & halfstrip_even>0 & wire_even>0")
    to_draw = "hs_lct_odd:eta_lct_odd"
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    s_bins = "(127,0,127)"
    w_bins = "(50,0,50)"
    p_bins = "(96,0,96)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    wirenBins  = int(w_bins[1:-1].split(',')[0])
    wmin = float(w_bins[1:-1].split(',')[1])
    wmax = float(w_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    LUT = TH2F("LUT","", stripnBins, smin, smax, wirenBins, wmin, wmax)
    
    C = TCanvas("C","",800,600)
# t.Draw(to_draw + ">>LUT", L_cut, "goff")
    t.Draw(to_draw, L_cut, "COLZ")
#    LUT.Draw("same")

    C.Print("%shalfstrip_wire%s"%(output_dir,ext))
     
#______________________________________________________________________________________
def PlothsGEMPad_1():
    input_dir = "files/"
    output_dir = "matching_Combined/"
    ext = ".png"
    trk_eff = "trk_eff_ME1b"
# total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    t = getTree("%s%s"%(input_dir, input_file1), trk_eff)
    
    L_cut = TCut("halfstrip_odd>0 & pad_odd>0 & pt>20")
#  L_cut = AND(L_cut,ok_clct1,ok_gdg1)
    to_draw = "pad_odd:halfstrip_odd"
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    s_bins = "(128,0,128)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    wirenBins  = int(w_bins[1:-1].split(',')[0])
    wmin = float(w_bins[1:-1].split(',')[1])
    wmax = float(w_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    topTitle = " " * 11 + "halfstrip VS Pad for odd chamber in ME1b" + " " *5 + "CMS Simulation Preliminary"
    LUT = TH2F("LUT", topTitle, stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT.SetXTitle("half-strip in ME1b")
    LUT.SetYTitle("Pad")
    
    C = TCanvas("C","",800,600)
    t.Draw(to_draw + ">>LUT", L_cut, "COLZ")
#    t.Draw(to_draw, L_cut, "COLZ")
#    LUT.Draw("same")

    C.Print("%shalfstrip_pad_ME1b_1%s"%(output_dir,ext))
     
#______________________________________________________________________________________
def PlothsGEMPad_2():
    input_dir = "files/"
    output_dir = "matching_Combined/"
    ext = ".png"
    
    trk_eff = "trk_eff_ME1b"
# total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    t = getTree("%s%s"%(input_dir, input_file1), trk_eff)
    L_cut = TCut("halfstrip_even>0 & pad_even>0 & pt>20")
    to_draw = "pad_even:halfstrip_even"
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    s_bins = "(128,0,128)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    wirenBins  = int(w_bins[1:-1].split(',')[0])
    wmin = float(w_bins[1:-1].split(',')[1])
    wmax = float(w_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    topTitle = " " * 11 + "halfstrip VS Pad for even chamber in ME1b" + " " *5 + "CMS Simulation Preliminary"
    LUT = TH2F("LUT", topTitle, stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT.SetXTitle("half-strip in ME1b")
    LUT.SetYTitle("Pad")
     
    C = TCanvas("C","",800,600)
    t.Draw(to_draw + ">>LUT", L_cut, "COLZ")
#    t.Draw(to_draw, L_cut, "COLZ")
#    LUT.Draw("same")

    C.Print("%shalfstrip_pad_ME1b_2%s"%(output_dir,ext))

#______________________________________________________________________________________
def PlothsAndmeanhs_1b():
    input_dir = "files/"
    output_dir = "matching_Combined/"
    ext = ".png"
    trk_eff = "trk_eff_ME1b"
# total numerator selection cut 
    ## the extra brackets around the extra_num_cut are necessary !!
    t = getTree("%s%s"%(input_dir, input_file1), trk_eff)
    
    L_cut = TCut("halfstrip_odd>0 & pad_odd>0 & pt>20")
#  L_cut = AND(L_cut,ok_clct1,ok_gdg1)
    to_draw = "pad_odd:halfstrip_odd"
 
    ## PyROOT works a little different than ROOT when you are plotting 
    ## histograms directly from tree. Hence, this work-around
    s_bins = "(128,0,128)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    wirenBins  = int(w_bins[1:-1].split(',')[0])
    wmin = float(w_bins[1:-1].split(',')[1])
    wmax = float(w_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    topTitle = " " * 11 + "halfstrip VS Mean halfstrip from simhits in ME1b" + " " *5 + "CMS Simulation Preliminary"
    LUT = TH2F("LUT", topTitle, stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT.SetXTitle("half-strip in ME1b")
    LUT.SetYTitle("mean hs")
    
    C = TCanvas("C","",800,600)
    t.Draw(to_draw + ">>LUT", L_cut, "COLZ")
#    t.Draw(to_draw, L_cut, "COLZ")
#    LUT.Draw("same")

    C.Print("%shalfstrip_meanhs_ME1b_1%s"%(output_dir,ext))
     

#______________________________________________________________________________________________
def HSToGEMPadLUT_1():

    f_1 = open('LUT/ME1b_HsToGemPad_1.txt','r')
    k = 0
    s_bins = "(128,0,128)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    LUT_low1 = TH2F("LUT","halfstrip VS pad LUT for odd chamber", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high1 = TH2F("LUT","halfstrip VS pad LUT for odd chamber", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_1:
    	HS_1 = int(line[1:-1].split()[2])
	LowPad_1 = int(line[1:-1].split()[6])
	HighPad_1 = int(line[1:-1].split()[10])
        LUT_low1.Fill(HS_1, LowPad_1)
    	LUT_high1.Fill(HS_1, HighPad_1)
    
    LUT_low1.SetXTitle("halfstrip in ME1b")
    LUT_low1.SetYTitle("Pad number")
    LUT_low1.SetMarkerColor(kRed)
    LUT_high1.SetMarkerColor(kBlue)
    LUT_low1.SetMarkerStyle(21)
    LUT_high1.SetMarkerStyle(20)
    """
#    print HS_1
    f_2 = open('LUT/HSToGEMPad_2.txt','r')
    LUT_low2 = TH2F("LUT","halfstrip VS pad LUT", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high2 = TH2F("LUT","halfstrip VS pad LUT", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_2:
    	HS_2 =int(line[1:-1].split()[2])
	LowPad_2 = int(line[1:-1].split()[6])
	HighPad_2 = int(line[1:-1].split()[10])
        LUT_low2.Fill(HS_2, LowPad_2)
    	LUT_high2.Fill(HS_2, HighPad_2)
    """
    C = TCanvas("C","C",800,600)
    LUT_low1.Draw()
    LUT_high1.Draw("same")

    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.AddEntry(LUT_low1, "HS to Low Pad","pl")
    leg.AddEntry(LUT_high1, "HS to High Pad","pl")
    leg.Draw("same")
    
    C.Print("%shalfstrip_to_GEMPad_ME1b_LUT1%s"%(output_dir,ext))


#______________________________________________________________________________________________
def HSToGEMPadLUT_2():

    f_2 = open('LUT/ME1b_HsToGemPad_2.txt','r')
    s_bins = "(128,0,128)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])

    LUT_low2 = TH2F("LUT","halfstrip VS pad LUT for even Chamber ", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high2 = TH2F("LUT","halfstrip VS pad LUT for even Chamber", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_2:
    	HS_2 =int(line[1:-1].split()[2])
	LowPad_2 = int(line[1:-1].split()[6])
	HighPad_2 = int(line[1:-1].split()[10])
        LUT_low2.Fill(HS_2, LowPad_2)
    	LUT_high2.Fill(HS_2, HighPad_2)


    LUT_low2.SetXTitle("halfstrip in ME1b")
    LUT_low2.SetYTitle("Pad number")
    LUT_low2.SetMarkerColor(kRed)
    LUT_high2.SetMarkerColor(kBlue)
    LUT_low2.SetMarkerStyle(21)
    LUT_high2.SetMarkerStyle(20)

    C = TCanvas("C","C",800,600)
    LUT_low2.Draw()
    LUT_high2.Draw("same")

    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.AddEntry(LUT_low2, "HS to Low Pad","pl")
    leg.AddEntry(LUT_high2, "HS to High Pad","pl")
    leg.Draw("same")
    
    C.Print("%shalfstrip_to_GEMPad_ME1b_LUT2%s"%(output_dir,ext))


#______________________________________________________________________________________________
def HSToGEMPadLUT_1_ME1a():

    f_1 = open('LUT/HsToGEMPad_1_E1_ME1a.txt','r')
    k = 0
    s_bins = "(96,0,96)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])
    
    LUT_low1 = TH2F("LUT","halfstrip VS pad LUT for ME1a of odd chamber for EndCap1", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high1 = TH2F("LUT","halfstrip VS pad LUT for ME1a of odd chamber for EndCap1", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_1:
    	HS_1 = int(line[1:-1].split()[2])
	LowPad_1 = int(line[1:-1].split()[6])
	HighPad_1 = int(line[1:-1].split()[10])
        LUT_low1.Fill(HS_1, LowPad_1)
    	LUT_high1.Fill(HS_1, HighPad_1)
    
    LUT_low1.SetXTitle("halfstrip in ME1a")
    LUT_low1.SetYTitle("Pad number")
    LUT_low1.SetMarkerColor(kRed)
    LUT_high1.SetMarkerColor(kBlue)
    LUT_low1.SetMarkerStyle(21)
    LUT_high1.SetMarkerStyle(20)
    """
#    print HS_1
    f_2 = open('LUT/HSToGEMPad_2.txt','r')
    LUT_low2 = TH2F("LUT","halfstrip VS pad LUT", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high2 = TH2F("LUT","halfstrip VS pad LUT", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_2:
    	HS_2 =int(line[1:-1].split()[2])
	LowPad_2 = int(line[1:-1].split()[6])
	HighPad_2 = int(line[1:-1].split()[10])
        LUT_low2.Fill(HS_2, LowPad_2)
    	LUT_high2.Fill(HS_2, HighPad_2)
    """
    C = TCanvas("C","C",800,600)
    LUT_low1.Draw()
    LUT_high1.Draw("same")

    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.AddEntry(LUT_low1, "HS to Low Pad","pl")
    leg.AddEntry(LUT_high1, "HS to High Pad","pl")
    leg.Draw("same")
    
    C.Print("%shalfstrip_to_GEMPad_E1_ME1a_LUT1%s"%(output_dir,ext))


#______________________________________________________________________________________________
def HSToGEMPadLUT_2_ME1a():

    f_2 = open('LUT/HsToGEMPad_2_E1_ME1a.txt','r')
    s_bins = "(96,0,96)"
    w_bins = "(50,0,50)"
    p_bins = "(200,0,200)"

    stripnBins  = int(s_bins[1:-1].split(',')[0])
    smin = float(s_bins[1:-1].split(',')[1])
    smax = float(s_bins[1:-1].split(',')[2])

    padnBins  = int(p_bins[1:-1].split(',')[0])
    pmin = float(p_bins[1:-1].split(',')[1])
    pmax = float(p_bins[1:-1].split(',')[2])

    LUT_low2 = TH2F("LUT","halfstrip VS pad LUT for ME1a of even Chamber for EndCap1", stripnBins, smin, smax, padnBins, pmin, pmax)
    LUT_high2 = TH2F("LUT","halfstrip VS pad LUT for ME1a of even chamber for EndCap1", stripnBins, smin, smax, padnBins, pmin, pmax)
    for line in f_2:
    	HS_2 =int(line[1:-1].split()[2])
	LowPad_2 = int(line[1:-1].split()[6])
	HighPad_2 = int(line[1:-1].split()[10])
        LUT_low2.Fill(HS_2, LowPad_2)
    	LUT_high2.Fill(HS_2, HighPad_2)


    LUT_low2.SetXTitle("halfstrip in ME1a")
    LUT_low2.SetYTitle("Pad number")
    LUT_low2.SetMarkerColor(kRed)
    LUT_high2.SetMarkerColor(kBlue)
    LUT_low2.SetMarkerStyle(21)
    LUT_high2.SetMarkerStyle(20)

    C = TCanvas("C","C",800,600)
    LUT_low2.Draw()
    LUT_high2.Draw("same")

    leg = TLegend(0.20,0.14,.67,0.4, "", "NDC");
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetTextFont(62)
    leg.AddEntry(LUT_low2, "HS to Low Pad","pl")
    leg.AddEntry(LUT_high2, "HS to High Pad","pl")
    leg.Draw("same")
    
    C.Print("%shalfstrip_to_GEMPad_E1_ME1a_LUT2%s"%(output_dir,ext))


    
#_______________________________________________________________________________________________
if __name__ == "__main__":
 
        input_file1 = "unflipped_Ana_1.root"
        input_dir = "files/"
        output_dir = "matching_Combined/"
        ext = ".png"

	if not os.path.exists(output_dir):
    		os.makedirs(output_dir)
#Plothseta()
#	PlothsGEMPad_1()	
#	PlothsGEMPad_2()
	PlothsAndmeanhs_1b()
        HSToGEMPadLUT_1()
        HSToGEMPadLUT_2()
#        HSToGEMPadLUT_1_ME1a()
#        HSToGEMPadLUT_2_ME1a()


    
