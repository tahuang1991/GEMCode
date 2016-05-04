from ROOT import *
import os


gROOT.SetBatch(1)
#gStyle from TStyle
gStyle.SetStatW(0.17)
gStyle.SetStatH(0.2)
gStyle.SetStatX(0.85)
gStyle.SetStatY(0.90)

gStyle.SetOptStat(111110)

gStyle.SetTitleStyle(0)
gStyle.SetTitleAlign(13) ## coord in top left
gStyle.SetTitleX(0.)
gStyle.SetTitleY(1.)
gStyle.SetTitleW(1)
#ROOT.gStyle.SetTitleTextColor(4)
gStyle.SetTitleXSize(0.05)
gStyle.SetTitleYSize(0.05)
gStyle.SetTitleH(0.058)
gStyle.SetTitleBorderSize(0)

gStyle.SetPadLeftMargin(0.126)
gStyle.SetPadRightMargin(0.14)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadBottomMargin(0.13)



def GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle, text, cuts,picname):
    
    c1 = TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    gStyle.SetOptFit(0111)
    gStyle.SetOptStat(0)

        
    chain = TChain(treename1)
    if os.path.isdir(filedir):
    	  ls = os.listdir(filedir)
    	  for x in ls:
		x = filedir[:]+x
		chain.Add(x)
    elif os.path.isfile(filedir):
	  chain.Add(filedir)
    else:
	  print " it is not file or dir ", filedir

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    b1 = TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("%s"%st_title)
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%xtitle)
    #b1.SetStats(0)

    b1.Sumw2() 
    chain.Draw(xaxis+">>b1",cuts)
   # b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)
    fb2 = TF1("fb3","gaus(0)+landau(3)",-.2,.2)
    #fb2.SetParameter(0,1.0)
    #fb2.SetParameter(1,0.00032)
    #fb2.SetParameter(2,0.00097)
    #fb2.SetParameter(4,-1.704e-4)
    #fb2.SetParameter(5,0.000301)
    b1.Draw() 
    #b1.Fit(fb2)
    b1.Fit("gaus")
    #b1.Fit("landau")

    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f,20<p_{T}<50}"%(st_title,etamin,etamax))
    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f}"%(st_title,etamin,etamax))
    tex1 = TLatex(0.13,.8,"%s"%text)
    tex1.SetNDC()
    tex1.SetTextFont(42)
    #tex1.SetTextColorAlpha(kRed, 0.6)
    tex1.Draw("same")
    #fitting= "p_{0}Exp(-#frac{(x-p_{1})^{2}}{2p_{2}^{2}})+p_{3}Landau(x,p_{4},p_{5})"
    fitting= "p_{0}Exp(#frac{-(x-p_{1})^{2}}{2p_{2}^{2}})"
    tex2 = TLatex(0.16,.25,"%s"%fitting)
    tex2.SetNDC()
    tex2.SetTextFont(42)
    #tex1.SetTextColorAlpha(kRed, 0.6)
    tex2.Draw("same")
    
    c1.SaveAs("%s"%picname+"_Gaus.pdf")
    c1.SaveAs("%s"%picname+"_Gaus.png")

def GEMCSCResolution2D(filedir,treename1,xaxis, yaxis, x_bins, y_bins, st_title, xtitle, ytitle, cuts,picname):
    
    c1 = TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    gStyle.SetOptFit(0111)
    gStyle.SetOptStat(0)

        
    chain = TChain(treename1)
    if os.path.isdir(filedir):
    	  ls = os.listdir(filedir)
    	  for x in ls:
		x = filedir[:]+x
		chain.Add(x)
    elif os.path.isfile(filedir):
	  chain.Add(filedir)
    else:
	  print " it is not file or dir ", filedir

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    b1 = TH2F("b1","b1",xBins,xminBin,xmaxBin, yBins, yminBin, ymaxBin)
    b1.SetTitle("%s"%st_title)
    b1.GetYaxis().SetTitle("%s"%ytitle)
    b1.GetXaxis().SetTitle("%s"%xtitle)
    #b1.SetStats(0)

    b1.Sumw2() 
    chain.Draw(yaxis+":"+xaxis+">>b1",cuts)
   # b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)

    b1.Draw("colz") 
    #b1.Fit("gaus")

    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f,20<p_{T}<50}"%(st_title,etamin,etamax))
    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f}"%(st_title,etamin,etamax))
    tex1 = TLatex(0.35,.9,"%s"%st_title)
    tex1.SetNDC()
    #tex1.SetTextColorAlpha(kRed, 0.6)
    #tex1.Draw("same")
    
    c1.SaveAs("%s"%picname+"_Gaus.pdf")
    c1.SaveAs("%s"%picname+"_Gaus.png")

yaxistitle="#Phi_{#vec{p}} at ME11"
yaxis="phi_sh-ptphi_sh"#phi from momentum
#filedir = "/eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt2_50_20160301/"
filedir = "GEMCSC_outana_0420_ctau1000.root"
#filedir = "GEMCSC_outana_0420_ctau100.root"
#filedir = "GEMCSC_outana_0420_ctau10.root"
filedir = "GEMCSC_outana_0425_prompt.root"
treename1="GEMCSCAnalyzer/trk_eff_CSC_ME11"
x_bins = "(100,-0.01,0.01)"
"""
xaxis="(phi_lct_odd-phi_cscsh_odd)"#phi from momentum
xtitle="(phi_lct_odd-phi_cscsh_odd)"#phi from momentum

picname="GEMCSC_phi_ME11_odd"
cuts ="abs(phi_lct_odd-phi_cscsh_odd)<1"
st_title = "#phi resolution in ME11 odd chamber"
GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle,cuts,picname)
xaxis = "(dphi_lct_odd)-(phi_cscsh_odd-phi_gemsh_odd)"
xtitle = "#Delta#Phi_{bending}"
picname="GEMCSC_deltaphi_Pt2_50_ME11_odd"
cuts ="abs((dphi_lct_odd)-(phi_cscsh_odd-phi_gemsh_odd))<1"
st_title = "GEM-CSC bending angle resolution in ME11 odd chamber"
GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle,cuts,picname)

#phi resolution
CSCMEs = ["ME11","ME12","ME13","ME21","ME22","ME31","ME32","ME41"]
evenodds = ["even","odd"]
for me in CSCMEs:
	for a in evenodds:
		treename1 = "GEMCSCAnalyzer/trk_eff_CSC_%s"%(me)
		xaxis = "phi_lct_%s-phi_cscsh_%s"%(a,a)
		xtitle = "phi_lct_%s-phi_cscsh_%s"%(a,a)
		picname = "GEMCSC_pt2_50_phi_%s_%s"%(me,a)
		st_title = "#phi resolution in %s %s chamber"%(me, a)
		cuts = "abs(phi_lct_%s-phi_cscsh_%s)<1"%(a,a)
		GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle,cuts,picname)
pts = [5,10,20,30,40]
"""
GEMCSCs = ["ME11","ME21"]
evenodds = ["even","odd"]
slope = {
   "ME11even": -0.132,
   "ME11odd": -0.059,
   "ME21even": -0.105,
   "ME21odd": -.064
	}
deltaphi0=-0.022
ddphi = 0.004
D_GEMCSC={
   "ME11even":45.789,
   "ME11odd":16.49,
   "ME21even":48.281,
   "ME21odd":23.58
}
for me in GEMCSCs:
	for a in evenodds:
		dd=-1
		#for pt in pts:
		while dd<10:
    			dd = dd+1
			#if dd!=5:
			#	continue
			#filedir = "/eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt%d_20160301/"%pt
			#filedir = "/eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt2_50_20160301/"
			treename1 = "GEMCSCAnalyzer/trk_eff_CSC_%s"%(me)
			deltaphi1 = deltaphi0+dd*ddphi
			deltaphi2 = deltaphi0+(dd+1)*ddphi
			print "deltaphi1 ",deltaphi1," deltaphi2 ",deltaphi2			
			#xaxis2 = "(phi_lct_%s-phi_pad_%s)-(phi_cscsh_%s-phi_gemsh_%s)"%(a,a,a,a)
  			case= "%s%s"%(me,a)
			xaxis2 = "10*((dphi_pad_%s)-(dphi_sh_%s)*(1+%f))/%f"%(a,a,slope[case],D_GEMCSC[case])
			#yaxis = "(dphi_pad_%s)"%(a)
			#xaxis = "(dphi_lct_%s)"%(a)
			yaxis = "(dphi_sh_%s)"%(a)
			xtitle = "#Delta#Phi_{GEM-CSC}^{L1}"
			xtitle2 = "#Delta#Phi_{GEM-CSC}^{L1} - #Delta#Phi_{GEM-CSC}^{SIM}"
			ytitle = "#Delta#Phi_{GEM-CSC}^{SIM}"
			picname = "DeltaPhiRes_Slice_ctau0_0420/GEMCSC1D_distanceCorr_L1andSim_%s_%s_ctau0_0425_%d"%(me,a,dd)
			picname2 = "DeltaPhiRes_Slice_ctau0_0420/GEMCSC1D_bending_distribution_Sim_%s_%s_ctau0_0425_%d"%(me,a,dd)
			st_title = "GEM-CSC bending angle distribution in %s %s chamber"%(me, a)
			"""
			xaxis = "abs(dphi_lct_%s)"%(a)
			xtitle = "#Delta#Phi_{(GEM,CSC)}, %s %s"%(me, a)
			picname = "DeltaPhi_Pt%d/GEMCSC_pt%d_absdeltaphi_%s_%s"%(pt,pt,me,a)
			st_title = "GEM-CSC bending angle in %s %s chamber, p_{T}=%d"%(me, a,pt)
			"""
			cuts = "abs(phi_lct_%s)<4 && abs(phi_pad_%s)<4 && abs(phi_cscsh_%s)<4 && abs(phi_gemsh_%s)<4"%(a,a,a,a)
			cuts_phi = "(dphi_sh_%s)>%f && (dphi_sh_%s)<=%f"%(a,deltaphi1,a,deltaphi2)
			x_bins = "(50,%f,%f)"%(deltaphi0+(dd-2)*ddphi, deltaphi0+(dd+3)*ddphi)
			y_bins = "(100,-0.005,0.005)"
			text = "prompt muon, %.3f<#Delta#Phi_{GEM-CSC}^{SIM} <%.3f"%(deltaphi1, deltaphi2)
			#text="c#tau=1000"
			GEMCSCResolution(filedir,treename1,xaxis2,y_bins,st_title, xtitle2, text, cuts+"&&"+cuts_phi,picname)
			#GEMCSCResolution(filedir,treename1,yaxis,y_bins,st_title, ytitle, text, cuts+"&&"+cuts_phi,picname2)
    			#dd = dd+1
			#GEMCSCResolution2D(filedir,treename1,xaxis, yaxis, x_bins, y_bins, st_title, xtitle, ytitle, cuts,picname)
			#GEMCSCResolution2D(filedir,treename1,"pt", yaxis,"(20,0,50)", y_bins, st_title, "simulated p_{T}", ytitle, cuts,"DeltaPhiRes2D_Pt2_50/GEMCSC2D_pt2_50_bending_deltaVssimPt_%s_%s"%(me,a))
