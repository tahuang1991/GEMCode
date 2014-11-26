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
    h1 = ROOT.TH1F("h1","h1",50,0,50)
    t.Draw("abs(pt) >> h1",den)
    h2 = ROOT.TH1F("h2","h2",50,0,50)
    t.Draw("abs(pt) >> h2",num)
    e = ROOT.TEfficiency(h2,h1)
    return e

b1 = ROOT.TH1F("b1","b1",50,0,50)
b1.GetYaxis().SetRangeUser(0.0,1)
b1.GetYaxis().SetTitleOffset(1.2)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("efficiency")
b1.GetXaxis().SetTitle("P_{T} of #mu track [GeV/c]")
b1.SetTitle(" "*8 +"LCT reco eff with different dphi cut"+" "*10 + "CMS Simulation Preliminary")
b1.SetStats(0)

treename = "GEMCSCAnalyzer/trk_eff_ME1b"
den = "(has_lct&1)>0 && abs(eta)<2.1 && (has_gem_pad&1)>0"
#den = "(has_lct&2)>0 && abs(eta)<2.1 && (has_gem_pad&2)>0"
num1 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (has_gem_pad&1)>0"
#num1 = "abs(eta)<2.1 && (has_lct&2)>0 && (has_lct&2)>0 && (has_gem_pad&2)>0"
#num2 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (pt>10)"#10GeV
#num2 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00596349)"#10GeV
num2 = "(has_gem_pad&1)>0 && abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.01023299)"#10GeV me11
#num2 = "(has_gem_pad&2)>0 && abs(eta)<2.1 && (has_lct&2)>0 && (has_lct&2)>0 && (abs(dphi_lct_even)<0.00458796)"#10GeV me11
#num3 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (pt>20)"#20GeV 
#num3 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00435298)"#20GeV
num3 = "(has_gem_pad&1)>0 && abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00535176)"#20GeV me11
#num3 = "(has_gem_pad&2)>0 && abs(eta)<2.1 && (has_lct&2)>0 && (has_lct&2)>0 && (abs(dphi_lct_even)<0.00276152)"#20GeV me11
#num4 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (pt>30)"#30GeV 
#num4 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00372145)"#30GeV
num4 = "(has_gem_pad&1)>0 && abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00389050)"#30GeV me11
#num4 = "(has_gem_pad&2)>0 && abs(eta)<2.1 && (has_lct&2)>0 && (has_lct&2)>0 && (abs(dphi_lct_even)<0.00204670)"#30GeV me11
#num5 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (pt>30)"#40GeV 
#num5 = "abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00372145)"#40GeV
num5 = "(has_gem_pad&1)>0 && abs(eta)<2.1 && (has_lct&1)>0 && (has_lct&1)>0 && (abs(dphi_lct_odd)<0.00310539)"#40GeV me11
#num5 = "(has_gem_pad&2)>0 && abs(eta)<2.1 && (has_lct&2)>0 && (has_lct&2)>0 && (abs(dphi_lct_even)<0.00170670)"#40GeV me11

e1 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna_fixnoclct.root",treename,den,num1)
e2 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna_fixnoclct.root",treename,den,num2)
e3 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna_fixnoclct.root",treename,den,num3)
e4 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna_fixnoclct.root",treename,den,num4)
e5 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna_fixnoclct.root",treename,den,num5)
#e5 = getEff("PU140_100k_2023_fixodd_GEMCSCAna.root",treename,den,num)
#e6 = getEff("PU140_100k_2023_fixodd_GEMCSCAna.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_PU140.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_PU140.root",treename,den,num)

e1.SetLineColor(ROOT.kRed)
e1.SetLineWidth(2)
e2.SetLineColor(ROOT.kOrange-3)
e2.SetLineWidth(2)
e4.SetLineColor(ROOT.kSpring-1)
e4.SetLineWidth(2)
e3.SetLineColor(ROOT.kAzure+1)
e3.SetLineWidth(2)
#e6.SetLineColor(ROOT.kViolet-2)
#e6.SetLineWidth(2)
e5.SetLineColor(ROOT.kPink+1)
e5.SetLineWidth(2)

b1.Draw()
e1.Draw("same")
e2.Draw("same")
e3.Draw("same")
e4.Draw("same")
e5.Draw("same")
#e6.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.26,0.16,0.76,0.46,"","brNDC")
#legend.SetFillColor(ROOT.kWhite)
legend.SetFillStyle(0)
#legend.SetBorderSize(0)
legend.SetTextSize(0.05)
#legend.SetHeader("PU140,Phase II detector: GEM-CSC local trigger")
legend.AddEntry(e1,"No dphi cut","l")
legend.AddEntry(e2,"10 Gev dphi cut","l")
legend.AddEntry(e3,"20 Gev dphi cut","l")
legend.AddEntry(e4,"30 Gev dphi cut","l")
legend.AddEntry(e5,"40 Gev dphi cut","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")
tex = ROOT.TLatex(.20,.5,"#splitline{PU140,odd chamber, abs(eta)<2.1}{Phase II detector:GEM-CSC local trigger}")
tex.SetTextSize(0.05)
tex.SetNDC()
tex.Draw("same")

c1.SaveAs("ME11_odd_LCT_reco_eff_pt_dphicut_Tao_PU140.pdf")
c1.SaveAs("ME11_odd_LCT_reco_eff_pt_dphicut_Tao_PU140.png")
