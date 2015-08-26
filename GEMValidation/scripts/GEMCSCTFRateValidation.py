import ROOT
import random
import os
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetStatW(0.07)
ROOT.gStyle.SetStatH(0.06)

ROOT.gStyle.SetOptStat(0)



ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

ROOT.gStyle.SetMarkerStyle(1)

c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()
ptbin = [2.0,   2.5,   3.0,   3.5,   4.0, 4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0, 16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0, 50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 200.0]
myptbin = np.asarray(ptbin)
#___________________________________________________
def getRatecount(tree,todraw,cut):
    
    htemp = ROOT.TH1F("htemp"," ",29,myptbin)
    tree.Draw(todraw+">>htemp",cut)
    return htemp.GetEntries()

#___________________________________________________
def getEventsNum(tree):
    events =0
    hnum = ROOT.TH1F("hnum"," ",3,-1,2)
    tree.Draw("firstl1mu>>hnum","firstl1mu>0")
    events = hnum.GetEntries()
    print " the entries of firstl1mu (firstl1mu>0) ",events
    return events
#____________________________________________________
def getRate(tree, cut):
    #f = ROOT.TFile(file)
    #t = f.Get(dir)
    h = ROOT.TH1F("h"," ",29,myptbin)
    n=1
    for x in ptbin:
#	print "cut ",cut+" && pt>=%f"%x
	content = getRatecount(tree,"pt",cut+"&& pt>=%f"%x)
    	#content = tree.GetEntries(cut+"&& pt>=%f"%x)
	#print "bin n ",n,"pt ",x ,"  content ",content
	h.SetBinContent(n, content)
	n= n+1
    h.Sumw2()
    print "before scale "
    #h.Print("all")
    ntotalEvents = getEventsNum(tree)
    h.Scale(40000./ntotalEvents/3.*0.795)
#    print "after scale "
 #   h.Print("all")
    return h


b1 = ROOT.TH1F("b1","b1",29,myptbin)
b1.GetYaxis().SetRangeUser(0.01,10000)
b1.GetYaxis().SetTitleOffset(1.2)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("Trigger Rate [kHz]")
b1.GetXaxis().SetTitle("L1 muon p_{T} [GeV]")
b1.SetTitle("CMS Simulation Preliminary"+"  "*24 +" PU140 ")
b1.SetStats(0)

treename = "L1TTriggerRate/evtree"
#tfile1 = ROOT.TFile("/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC25_patch1/src/GEMCode/GEMValidation/test/Rate_Neutrino_SLHC25_2019withoutGEM_PU140_0706.root")
#tfile2 = ROOT.TFile("/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC25_patch1/src/GEMCode/GEMValidation/test/Rate_Neutrino_SLHC25_2019withGEM_PU140_0706.root")
#tfile3 = ROOT.TFile("/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC25_patch1/src/GEMCode/GEMValidation/test/Rate_Neutrino_SLHC25_2023NoRPC_PU140_0706.root")
#tfile4 = ROOT.TFile("/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC25_patch1/src/GEMCode/GEMValidation/test/Rate_Neutrino_SLHC25_PU140_0706.root")
tfile1 = ROOT.TFile("Rate_Neutrino_SLHC25_2019withoutGEM_PU140_0706.root")
tfile2 = ROOT.TFile("Rate_Neutrino_SLHC25_2019withGEM_PU140_0706.root")
tfile3 = ROOT.TFile("Rate_Neutrino_SLHC25_2023NoRPC_PU140_0706.root")
tfile4 = ROOT.TFile("Rate_Neutrino_SLHC25_PU140_0706.root")
tree1 = tfile1.Get(treename)
tree2 = tfile2.Get(treename)
tree3 = tfile3.Get(treename)
tree4 = tfile4.Get(treename)
#tree.Print()
cut1 = "nl1tracks>0 && eta<2.14 && eta>1.64 && nstubs>2 && hasME1"
cut2 = "nl1tracks>0 && eta<2.14 && eta>1.64 && nstubs>2 && passGE11"
e1 = getRate(tree1, cut1)
e2 = getRate(tree2, cut2)
e3 = getRate(tree3, cut2)
e4 = getRate(tree4, cut2)
#e4 = getRate("/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC25_patch1/src/GEMCode/GEMValidation/test/Rate_Neutrino_SLHC25_PU140_0706.root",treename, cut)

#e0 = getAllEff("/eos/uscms/store/user/tahuang/SLHC25_patch1_2023Muon_1M_Ana_PU140_Pt2_50_ME11_step0/",treename, den, num)

e4.SetFillColor(ROOT.kBlue+1)
e3.SetFillColor(ROOT.kRed-4)
e2.SetFillColor(ROOT.kMagenta+2)
e1.SetFillColor(ROOT.kGreen+2)

b1.Draw()
e1.Draw("e3same")
e2.Draw("e3same")
e3.Draw("e3same")
e4.Draw("e3same")
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()


legend = ROOT.TLegend(0.5,0.55,0.93,0.92)
legend.SetFillStyle(0)
#legend.SetFillColor(ROOT.kWhite)
legend.SetTextFont(62)
#legend.SetTextSize()
#legend.SetHeader("p_{T}^{sim}>10,2.14>|#eta|>1.64, has at least 3stubs and hasME1")
legend.AddEntry(e1,"#splitline{PhaseI, 3+stubs and}{ one stub in ME11}","f")
legend.AddEntry(e2,"#splitline{PhaseII with GE11 only}{ 3+stubs and YE1/1 bending angle}","f")
legend.AddEntry(e3,"#splitline{PhaseII with GE11GE21 only}{3+stubs and YE1/1 bending angle}","f")
legend.AddEntry(e4,"#splitline{Full PhaseII, 3+stubs }{and YE1/1 bending angle}","f")
legend.Draw("same")

tex = ROOT.TLatex(0.15,0.65,"1.64<|#eta|<2.14")
tex.SetTextSize(0.05)
tex.SetTextFont(42)
tex.SetNDC()
tex.Draw("same")


c1.SaveAs("TFtrack_rate_4scenarios_3stubs_eta2_GE11_PU140_eta_pendcap_0810.pdf")
c1.SaveAs("TFtrack_rate_4scenarios_3stubs_eta2_GE11_PU140_eta_pendcap_0810.png")
