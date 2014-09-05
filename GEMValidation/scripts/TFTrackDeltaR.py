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

#____________________________________________________________________________
def getHistInTree(file,dir,cut):
    f = ROOT.TFile(file)
    t = f.Get(dir)
    h = ROOT.TH1F("h","h",100,0,0.5)
    t.Draw("deltaR >> h",cut)
    print "file: ",file
    print "h type ", type(h)," t type ", type(t)
#ROOT.SetOwnership(h, False)
    c_temp = ROOT.TCanvas()
    h.Draw()
    c_temp.SaveAs("temp.png")
    return h

b1 = ROOT.TH1F("b1","b1",100,0,0.5)
b1.GetYaxis().SetNdivisions(520)
b1.GetYaxis().SetTitle("counts")
b1.GetXaxis().SetTitle("deltaR")
b1.SetTitle(" "*36 + "CMS Simulation Preliminary")
b1.SetStats(0)



treename = "GEMCSCAnalyzer/trk_eff_ALL"
num = " (has_tfTrack)>0 && pt>10"
#print "treename ",treename,"  type ",type(treename)

e1 = getHistInTree("PU140_200k_Pt2-50_2023_correction_GEMCSC.root",treename,num)
#e2 = getEff("PU140_100k_2019withoutGEM_GEMCSCAna.root",treename,den,num)
e4 = getHistInTree("PU140_200k_Pt2-50_2023_lct_GEMCSC.root",treename,num)
#e3 = getEff("PU140_Jason_100k_2023_delta_fixdt_GEMCSCAna.root",treename,den,num)
e2 = getHistInTree("PU140_10k_Pt20_2023_update_GEMCSC.root",treename,num)

print "e1", type(e1)
#e5 = getEff("GSA_GEMCSC_Step2_Com_PU140.root",treename,den,num)
#e6 = getEff("GSA_GEMCSC_Step3_Com_PU140.root",treename,den,num)
#e7 = getEff("GSA_GEMCSC_Step4_Com_PU140.root",treename,den,num)
ROOT.gPad.SetLogy()
b1.Draw()
e1.Draw("same")
e4.Draw("same")
e2.Draw("same")
#e7.Draw("same")

legend = ROOT.TLegend(0.20,0.23,.75,0.5, "", "brNDC")
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.05)
legend.SetHeader("PU140,  P_{T} > 10 GeV")
legend.AddEntry(e1,"simeta, propagated phi","f")
legend.AddEntry(e2,"Jason's ana","f")
#legend.AddEntry(e3,"simeta, simphi","f")
legend.AddEntry(e4,"eta,phi from stubs","f")
#legend.AddEntry(e4,"Upgrade TMB:CSC-RPC local trigger","f")
#legend.AddEntry(e5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(e6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(e7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")

#tex = ROOT.TLatex(0.25,0.5,"PU140, P\_T > 10Gev")
#tex.SetTextSize(0.05)
#tex.SetNDC()
#tex.Draw("same")

c1.SaveAs("TFTrack_200k_deltaR_PU140.pdf")
c1.SaveAs("TFTrack_200k_deltaR_PU140.png")
