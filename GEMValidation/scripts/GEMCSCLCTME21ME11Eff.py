import ROOT
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

def getEff(file,dir,den,num):
    f = ROOT.TFile(file)
    t = f.Get(dir)
    h1 = ROOT.TH1F("h1","h1",100,1.5,2.5)
    t.Draw("abs(eta) >> h1",den)
    h2 = ROOT.TH1F("h2","h2",100,1.5,2.5)
    t.Draw("abs(eta) >> h2",num)
    e = ROOT.TEfficiency(h2,h1)
    return e

b1 = ROOT.TH1F("b1","b1",100,1.5,2.5)
b1.GetYaxis().SetRangeUser(0.0,1.02)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("Efficiency")
b1.GetXaxis().SetTitle("Simulated muon #eta")
b1.SetTitle(" "*12 +"YE1/1 stub reconstruction"+" "*10 + "CMS Phase-II Simulation Preliminary")
b1.SetStats(0)

treename = "GEMCSCAnalyzer/trk_eff_ME11"
den = "has_csc_sh>0 && pt >10"
num = "has_csc_sh>0 && has_lct>0 && pt>10"
num1 = den+" && "+"(has_alct>0 || has_clct>0)"
#e1 = getEff("GEMCSC_Ana_Test1_CSC4.root",treename,den,num)
e3 = getEff("PU140_100k_2019withoutGEM_GEMCSCAna.root",treename,den,num)
#e4 = getEff("PU140_100k_2019withGEM_GEMCSCAna.root",treename,den,num)
#e3 = getEff("PU140_Jason_100k_2023_delta_fixdt_GEMCSCAna.root",treename,den,num)
e4 = getEff("PU140_100k_2023_fixeven_GEMCSCAna.root",treename,den,num)
e1 = getEff("PU140_100k_2023_fixeven_GEMCSCAna.root",treename,den,num1)

#e5 = getEff("GSA_GEMCSC_Step2_Com_PU140.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_Com_PU140.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_Com_PU140.root",treename,den,num)

e1.SetFillColor(ROOT.kBlack)
#e1.SetLineWidth(2)
#e2.SetLineColor(ROOT.kRed)
#e2.SetLineWidth(2)
e3.SetFillColor(ROOT.kRed)
#e3.SetLineWidth(2)
e4.SetFillColor(ROOT.kAzure-1)
#e1.SetFillColor(ROOT.kRed)
#e4.SetLineWidth(2)
"""
e6.SetLineColor(ROOT.kViolet-2)
e6.SetLineWidth(2)
e7.SetLineColor(ROOT.kPink+1)
e7.SetLineWidth(2)
"""
b1.Draw()
#e1.Draw("e3same")
#e2.Draw("same")
e3.Draw("e3same")
e4.Draw("e3same")
#e5.Draw("same")
#e6.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.20,0.15,.90,0.40, "", "brNDC")
legend.SetBorderSize(0)
#legend.SetFillStyle(0)
legend.SetFillColor(ROOT.kWhite)
#legend.SetTextSize(0.05)
legend.SetHeader("PU140, P_{T} > 10GeV")
#legend.AddEntry(e2,"CSC only SLHC Algorithm in 2019","l")
legend.AddEntry(e3,"Phase I detector:CSC only local trigger","f")
#legend.AddEntry(e3,"SLHC TMB:CSC only local trigger","l")
#legend.AddEntry(e3,"GEM-CSC Algorithm in 2023 Jason","l")
legend.AddEntry(e4,"Phase II detector:GEM-CSC local trigger","f")
#legend.AddEntry(e1,"Phase II detector:potentially achievable","f")
#legend.AddEntry(e5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")

c1.SaveAs("LCT_100k_ME11_reco_eff_com_TP_PU140.pdf")
c1.SaveAs("LCT_100k_ME11_reco_eff_com_TP_PU140.png")