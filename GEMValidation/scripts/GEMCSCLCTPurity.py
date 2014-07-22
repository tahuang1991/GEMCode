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
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

    
c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()

def dphiCut(h, fractionToKeep):
    """Get the dPhi cut corresponding to the fractionToKeep [95,98,99]"""

    if h.GetEntries() == 0:
        return 2.0000
    
    ax = h.GetXaxis()
    total = h.Integral()
    bin = 1
    fractionToKeep = fractionToKeep/100.
    for b in range(ax.GetNbins()):
        if (h.Integral(0,b)/total > fractionToKeep):
            bin = b - 1
            break

    ## interpolate
    x1 = ax.GetBinUpEdge(bin)
    x2 = ax.GetBinUpEdge(bin + 1)
    y1 = h.Integral(0, bin)/total
    y2 = h.Integral(0, bin + 1)/total
    x = x1 + (fractionToKeep - y1)/(y2-y1)*(x2-x1)
    return x

#_______________________________________________________
def getPurity(file,dir,den,num,xaxis,h_bins):
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",nBins,minBin,maxBin)
    b1.GetYaxis().SetRangeUser(0.6,1.02)
    b1.GetYaxis().SetTitleOffset(1.2)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("LCT Efficiency")
    b1.GetXaxis().SetTitle("%s of track"%xaxis)
    b1.SetTitle("LCT reco"+" "*16 + "CMS Simulation Preliminary")
    b1.SetStats(0)

    h1 = ROOT.TH1F("h1","h1",nBins,minBin,maxBin)
    t.Draw("abs(%s) >> h1"%xaxis,den)
    h2 = ROOT.TH1F("h2","h2",nBins,minBin,maxBin)
    t.Draw("abs(%s) >> h2"%xaxis,num)
    e = ROOT.TEfficiency(h2,h1)
    
    b1.Draw()
    e.Draw("same")
    legend = ROOT.TLegend(0.23,0.20,0.80,0.32)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetHeader("PU140,simPt>40Gev,dphicut on simPt")
#legend.AddEntry(e1,"","l")
    legend.Draw("same")
	
    c1.SaveAs("LCT_reco_ME21_dphicut_%s_Pt40_PU140.pdf"%xaxis)
    c1.SaveAs("LCT_reco_ME21_dphicut_%s_Pt40_PU140.png"%xaxis)
#____________________________________________________________________
def draw2D(file,dir,num,xaxis,yaxis,x_bins,y_bins):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b1.GetYaxis().SetTitle("%s"%yaxis)
    b1.GetXaxis().SetTitle("%s"%xaxis)
    b1.SetTitle("LCT reco in ME21, PU140"+" "*16 + "CMS Simulation Preliminary")
    b1.SetStats(0)
    todraw = "%s"%yaxis+":"+"%s>>b1"%xaxis
    t.Draw(todraw,num,"colz")
    b1.SetMaximum(150)
    b1.Draw("colz")
    ROOT.gPad.SetLogz() 
    legend = ROOT.TLegend(0.15,0.56,0.40,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetHeader("PU140")
#legend.AddEntry(e1,"","l")
# legend.Draw("same")
    tex = ROOT.TLatex(0.30,0.30,"Pt>10Gev,even chamber")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")
	
    c1.SaveAs("LCT_reco_ME21_%s"%xaxis+"_VS_%s_match_even_PU140.pdf"%yaxis)
    c1.SaveAs("LCT_reco_ME21_%s"%xaxis+"_VS_%s_match_even_PU140.png"%yaxis)



#___________________________________________
def gemdphi(file,dir,cut,xaxis,x_bins,pt_min,pt_max):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.GetYaxis().SetTitle("counts")
    b1.GetXaxis().SetTitle("%s"%xaxis)
    b1.SetTitle("LCT reco"+" "*16 + "CMS Simulation Preliminary")
    b1.SetStats(0)
    todraw = "%s>>b1"%xaxis
    t.Draw(todraw,cut,"colz")
    b1.Draw()
    b1.SetBinContent(0,0) 
    b1.SetBinContent(xBins+1,0) 
    dphi_cut = dphiCut(b1,98)
    print "dphicut for simtrack pt %s, 98 "%pt_min, dphi_cut
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
    legend.SetHeader("PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max)
#legend.AddEntry(e1,"","l")
#    legend.Draw("same")
    line1 = "PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max
    line2 = "98 dphicut %f"%dphi_cut
    tex = ROOT.TLatex(0.15,0.45,line1)
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")
    tex2 = ROOT.TLatex(0.15,0.35,line2)
    tex2.SetTextSize(0.05)
    tex2.SetNDC()
    tex2.Draw("same")
	
    c1.SaveAs("LCT_reco_pt_%s"%xaxis+"_%s"%pt_min+"_%s_PU140.pdf"%pt_max)
    c1.SaveAs("LCT_reco_pt_%s"%xaxis+"_%s"%pt_min+"_%s_PU140.png"%pt_max)



#_______________________________________________________________________________
if __name__ == "__main__":

    file = "PU140_100k_2023_LCTTree_GEMCSC.root"
#file = "PU140_100k_2023_FlippedME1a_GEMCSC.root"
    dir = "GEMCSCAnalyzer/LCTTree"
    den = "station==2&& ring==1 && pt>40 &&match>0 "
    num = "station==2&& ring==1 && pt>40 &&match>0 && passdphicut>0 "

    xaxis1 = "pt"
    Pt_bins = "(50,0,50)"
# getPurity(file,dir,den,num,xaxis1,Pt_bins)

    xaxis2 = "eta"
    Eta_bins = "(40,1.5,2.5)"
#   getPurity(file,dir,den,num,xaxis2,Eta_bins)
    
    Eta_bins1 = "(40,-2.5,-1.5)"
    DPhi_bins = "(100,-0.25,0.25)"
    WG_bins = "(50,0,50)"
    WG_bins2 = "(96,0,96)"
    Hs_bins = "(128,0,128)"
    Hs_bins2 = "(96,0,96)"
    Hs_bins3 = "(160,0,160)"
    xaxis = "hs"
    yaxis = "dphi"
    cut1 = "pt>10 && station==2 && match>0 && ring==1 && (chamber%2)==0"
    draw2D(file,dir,cut1,xaxis,yaxis,Hs_bins3,DPhi_bins)
    
    xaxis_pt = "truept"
    yaxis_pt = "pt"
    cut2 = "eta<0 && match>0"
#    draw2D(file,dir,cut2,xaxis_pt,yaxis_pt,Pt_bins,Pt_bins)
    
    cut3 = "pt>10 && station==1 && abs(eta)>2.05"
    pt_min = "10"
    pt_max = "50"
    gemdphige11 = "dphi"
    DPhi_bins = "(100,-0.5,0.5)"
#    gemdphi(file,dir,cut3,gemdphige11,DPhi_bins,pt_min,pt_max)
