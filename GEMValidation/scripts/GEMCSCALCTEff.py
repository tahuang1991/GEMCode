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
    t.Draw("(eta) >> h1",den)
    h2 = ROOT.TH1F("h2","h2",40,1.5,2.5)
    t.Draw("(eta) >> h2",num)
    e = ROOT.TEfficiency(h2,h1)
    return e

b1 = ROOT.TH1F("b1","b1",40,1.5,2.5)
b1.GetYaxis().SetRangeUser(0.60,1.02)
b1.GetYaxis().SetTitleOffset(1.2)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("ALCT reconstruction efficiency")
b1.GetXaxis().SetTitle("#eta of simulated muon track")
b1.SetTitle("ALCT Reco in ME41(Positive Endcap)"+" "*10 + "CMS Simulation Preliminary")
b1.SetStats(0)

treename = "GEMCSCAnalyzer/trk_eff_ME41"
den = "has_csc_sh>0 "
num = "has_csc_sh>0 && has_alct>0"

#e1 = getEff("GEMCSC_Ana_Test1_CSC4.root",treename,den,num)
#e2 = getEff("GEMCSC_Ana_Test1_Wg3.root",treename,den,num)
#e3 = getEff("GEMCSC_Ana_Test1_Wg3fixed.root",treename,den,num)
e4 = getEff("PU140_100k_2023_delta_fixdt_GEMCSCAna.root",treename,den,num)
#e5 = getEff("GSA_GEMCSC_Step2_PU140.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_PU140.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_PU140.root",treename,den,num)
e4.SetLineColor(ROOT.kAzure-1)
e4.SetLineWidth(2)
"""
e6.SetLineColor(ROOT.kViolet-2)
e6.SetLineWidth(2)
e7.SetLineColor(ROOT.kPink+1)
e7.SetLineWidth(2)
"""
b1.Draw()
#e1.Draw("same")
#e2.Draw("same")
e4.Draw("same")
#e5.Draw("same")
#e6.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.23,0.13,0.62,0.32)
legend.SetFillColor(ROOT.kWhite)
legend.SetHeader("PU140")
legend.AddEntry(e4,"requirement for ALCT >=4 hits","l")
#legend.AddEntry(e2,"requirement for ALCT >=3 hits before fixed quality","l")
#legend.AddEntry(e3,"requirement for ALCT >=3 hits","l")
#legend.AddEntry(e4,"GEM-CSC Algorithm with all improvement","l")
#legend.AddEntry(e5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")

c1.SaveAs("ALCT_100k_ME41_reco_eff_Positive_new_PU140.pdf")
c1.SaveAs("ALCT_100k_ME41_reco_eff_Positive_new_PU140.png")
