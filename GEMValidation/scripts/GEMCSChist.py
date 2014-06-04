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
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

ROOT.gStyle.SetMarkerStyle(1)

c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()

def getTree(file,dir):
    f = ROOT.TFile(file)
    t = f.Get(dir)
    return t

	    
treename = "trk_eff_ME41"	    
tree = getTree("PU140_100k_2023_fixeven_GEMCSCAna.root",treename)
h1 = ROOT.TH1F("h1","h1",40,1.5,2.5)
tree.Draw("(-eta) >> h1","nlayers_csc_sh_odd>=3")
h2 = ROOT.TH1F("h2","h2",40,1.5,2.5)
tree.Draw("(-eta) >> h2","nlayers_csc_sh_odd>3")
h2 = ROOT.TH1F("h2","h2",40,1.5,2.5)
tree.Draw("(-eta) >> h2","nlayers_csc_sh_odd>=3 && nlayers_wg_dg_odd ==3")

h1.Draw()
h2.Draw("same")


legend = ROOT.TLegend(0.23,0.13,0.72,0.45)
legend.SetFillColor(ROOT.kWhite)
legend.SetHeader("PU140, with pt>10")
#legend.AddEntry(h1,"CSC SLHC Algorithm (4 hits)","l")
#legend.AddEntry(h2,"CSC only SLHC Algorithm in 2019","l")
#legend.AddEntry(h3,"GEM-CSC Algorithm in 2019+GE11","l")
#legend.AddEntry(h3,"GEM-CSC Algorithm in 2023 Jason","l")
#legend.AddEntry(h4,"GEM-CSC Algorithm in 2023 after bugfix","l")
#legend.AddEntry(h5,"GEM-CSC Algorithm (step 2)","l")
#legend.AddEntry(h6,"GEM-CSC Algorithm (step 3)","l")
#legend.AddEntry(h7,"GEM-CSC Algorithm (step 4)","l")
legend.Draw("same")

c1.SaveAs("100k_ME41_hist_Negative_fixeven_PU140.pdf")
c1.SaveAs("100k_ME41_hist_Negative_fixeven_PU140.png")
