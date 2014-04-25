import ROOT

c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()

def getEff(file,dir,den,num):
    f = ROOT.TFile(file)
    t = f.Get(dir)
    h1 = ROOT.TH1F("h1","h1",40,1.5,2.5)
    t.Draw("abs(eta) >> h1",den)
    h2 = ROOT.TH1F("h2","h2",40,1.5,2.5)
    t.Draw("abs(eta) >> h2",num)
    e = ROOT.TEfficiency(h2,h1)
    return e

b1 = ROOT.TH1F("b1","b1",40,1.5,2.5)
b1.GetYaxis().SetRangeUser(0.60,1.02)
b1.GetYaxis().SetTitleOffset(1.2)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("CLCT reconstruction efficiency")
b1.GetXaxis().SetTitle("#eta of simulated muon track")
b1.SetTitle(" "*10+"CLCT Reco in ME21"+" "*34 + "CMS Simulation Preliminary")
b1.SetStats(0)

treename = "GEMCSCAnalyzer/trk_eff_ME21"
den = "has_csc_sh>0 "
num = "has_csc_sh>0 && has_clct>0"

e1 = getEff("GEMCSC_Ana_Test1_CSC4.root",treename,den,num)
e2 = getEff("GEMCSC_Ana_Test1_step0.root",treename,den,num)
e3 = getEff("GEMCSC_Ana_Test1_step1.root",treename,den,num)
e4 = getEff("GEMCSC_Ana_Test1_all.root",treename,den,num)
#e5 = getEff("GSA_GEMCSC_Step2_PU0.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_PU0.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_PU0.root",treename,den,num)

e1.SetLineColor(ROOT.kBlack)
e1.SetLineWidth(2)
e2.SetLineColor(ROOT.kRed)
e2.SetLineWidth(4)
e4.SetLineColor(ROOT.kAzure-1)
e4.SetLineWidth(2)
e3.SetLineColor(ROOT.kMagenta)
e3.SetLineWidth(2)
"""
e6.SetLineColor(ROOT.kViolet-2)
e6.SetLineWidth(2)
e7.SetLineColor(ROOT.kPink+1)
e7.SetLineWidth(2)
"""
b1.Draw()
e1.Draw("same")
e2.Draw("same")
e3.Draw("same")
e4.Draw("same")
#e5.Draw("same")
#e6.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.23,0.13,0.52,0.40)
legend.SetFillColor(ROOT.kWhite)
legend.SetHeader("PU0")
legend.AddEntry(e1,"CSC SLHC Algorithm (4 hits)","l")
legend.AddEntry(e2,"CSC SLCH Algorithm (3 hits)","l")
legend.AddEntry(e3,"GEM-CSC Algorithm including CLCT from lowQ cclct","l")
legend.AddEntry(e4,"GEM-CSC Algorithm with all improvement","l")
#legend.AddEntry(e5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")

c1.SaveAs("CLCT_ME21_reco_eff_PU0.pdf")
c1.SaveAs("CLCT_ME21_reco_eff_PU0.png")
