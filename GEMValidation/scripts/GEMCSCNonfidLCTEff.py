import ROOT

ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.04)

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
b1.GetYaxis().SetRangeUser(0.0,0.15)
#b1.GetYaxis().SetTitleOffset(1.0)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("inefficiency")
b1.GetXaxis().SetTitle("#eta of simulated muon track")
b1.SetTitle(" "*10+ "Non-Fiducial  Inefficiency due to noalct and noclct in ME21"+" "*20 + "CMS Simulation Preliminary")
b1.SetStats(0)

treename = "GEMCSCAnalyzer/trk_eff_ME21"
den = "nlayers_csc_sh_even>=3 && pt>10"
#num = "nlayers_csc_sh_even>=3 && ((has_alct&2)==0 || ((has_clct&2)==0 && (has_rpc_dg&2)==0) || ((has_clct&2)==0 && quality_alct_even==0 && (has_rpc_dg&2)==1) ||((quality_alct_even==0 || quality_clct_even==3) && (has_rpc_dg&2)==0) )"
#num = "nlayers_csc_sh_even>=3 && ((quality_alct_even==0 || quality_clct_even==3)&&((has_rpc_dg&2)==0))"

#num = "nlayers_csc_sh_even>=3 && ((has_clct&2)==0 && quality_alct_even==0 && (has_rpc_dg&2)==1)"
#num = "nlayers_csc_sh_even>=3 && ((has_clct&2)==0 && (has_rpc_dg&2)==0)"
#num = "nlayers_csc_sh_even>=3 && (has_alct&2)==0"
#num = den+" && (has_lct&2)>0"
num1 = "((has_alct&2)==0 && (has_clct&2)==0)"
num = den + "&&" + num1
print num


#e1 = getEff("GEMCSC_Ana_Test1_CSC4.root",treename,den,num)
#e2 = getEff("GEMCSC_Ana_Test1_Wg3.root",treename,den,num)
#e3 = getEff("100k_Ana_PU140_UpgradeME11ME21.root",treename,den,num)
e4 = getEff("Nlayers_PU140_100k_2023_fixeven_GEMCSCAna.root",treename,den,num)
#e5 = getEff("GSA_GEMCSC_Step2_PU140.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_PU140.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_PU140.root",treename,den,num)

#e1.SetLineColor(ROOT.kBlack)
#e1.SetLineWidth(4)
#e2.SetLineColor(ROOT.kRed)
#e2.SetLineWidth(2)
e4.SetLineColor(ROOT.kAzure-1)
e4.SetLineWidth(2)
#e3.SetLineColor(ROOT.kMagenta)
#e3.SetLineWidth(2)
"""
e6.SetLineColor(ROOT.kViolet-2)
e6.SetLineWidth(2)
e7.SetLineColor(ROOT.kPink+1)
e7.SetLineWidth(2)

"""
#c1.SetLogy()
b1.Draw()
#e1.Draw("same")
#e2.Draw("same")
#e3.Draw("same")
e4.Draw("same")
#e5.Draw("same")
#e6.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.23,0.70,0.72,0.80)
legend.SetFillColor(ROOT.kWhite)
legend.SetHeader("PU140, even chamber")
#legend.AddEntry(e4,"requirement for Anode digi >= 3 layers","l")
#legend.AddEntry(e2,"requirement for ALCT >3 hits before fixed quality","l")
#legend.AddEntry(e3,"requirement for ALCT >= 3 hits","l")
#legend.AddEntry(e4,"GEM-CSC Algorithm with all improvement","l")
#legend.AddEntry(e5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
#legend.Draw("same")
tex = ROOT.TLatex(.25,.5,"PU140, even chamber, P_{T}>10GeV")
tex.SetTextSize(0.05)
tex.SetNDC()
tex.Draw("same")

c1.SaveAs("Nonfid_100k_ME21_ineff_min_even_PU140.pdf")
c1.SaveAs("Nonfid_100k_ME21_ineff_min_even_PU140.png")
