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



def GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle, cuts,text,picname, fitting):
    
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
    		if os.path.isdir(x):
			continue
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

    chain.Draw(xaxis+">>b1",cuts)
    
    #b1.Sumw2() 
    #b1.Print("ALL")
   # b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)

    b1.Draw() 
    g1 = TF1("g1","gaus",-2*b1.GetStdDev(),2*b1.GetStdDev())
    if fitting:
    	b1.Fit(g1,"LRQ")

    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f,20<p_{T}<50}"%(st_title,etamin,etamax))
    #tex1 = TLatex(0.35,.8,"#splitline{%s}{%.1f<|#eta|<%.1f}"%(st_title,text))
    tex1 = TLatex(0.2,.8,"%s"%text)
    tex1.SetNDC()
    #tex1.SetTextColorAlpha(kRed, 0.6)
    tex1.SetTextColor(kRed)
    tex1.Draw("same")
    
    c1.SaveAs("%s"%picname+".pdf")
    c1.SaveAs("%s"%picname+".png")
    c1.SaveAs("%s"%picname+".C")
    #x1 = xaxis.split('-').[0]
    #c2 = TCanvas()
    x1 = xaxis.split('-')[0]
    x2 = xaxis.split('-')[1]
    xtitle1 = xtitle.split('-')[0]
    xtitle2 = xtitle.split('-')[1]
    
    #b2 = TH1F("b2","b2",xBins,xminBin,xmaxBin)
    #b2.SetTitle("%s distribution "%xtitle2)
    #b2.GetYaxis().SetTitle("Events")
    #b2.GetXaxis().SetTitle("%s"%xtitle2)
    #b1.SetStats(0)

    #chain.Draw(x2+">>b2",cuts)

    #b2.Draw()
    #tex1.Draw("same")
    #c2.SaveAs("%s_%s_distribution"%(picname,x2)+".png")
    
yaxistitle="#Phi_{#vec{p}} at ME11"
yaxis="phi_sh-ptphi_sh"#phi from momentum
filedir = "out_ana.root"
filedir = "/eos/uscms/store/user/tahuang/SLHC26_patch1_2023Muon_1M_Ana_PU0_Pt20_0721/"
filedir = "/eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/GEMCSCAna_ctau0_Pt2_50_0816_v5sim/160816_222241/0000/"
filedir = "/eos/uscms/store/user/tahuang/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/GEMCSCAna_ctau0_Pt2_50_0827/160828_002000/0000/"
#filedir2 = "/eos/uscms/store/user/tahuang/DarkSUSY_MH-125_MGammaD-20000_ctau1000_14TeV_madgraph-pythia6-tauola/GEMCSCAna_DarkSUSY_ctau1000_0816_v5sim/160816_222334/0000/"
filedir2_old = "/eos/uscms/store/user/tahuang/DarkSUSY_MH-125_MGammaD-20000_ctau1000_14TeV_madgraph-pythia6-tauola/GEMCSCAna_DarkSUSY_ctau1000_0826/160826_215302/0000/"
#filedir2_old = "/eos/uscms/store/user/tahuang/DarkSUSY_MH-125_MGammaD-20000_ctau1000_14TeV_madgraph-pythia6-tauola/GEMCSCAna_DarkSUSY_ctau1000_0827/160828_001931/0000/"
filedir2 = "/eos/uscms/store/user/tahuang/DarkSUSY_MH-125_MGammaD-20000_ctau1000_14TeV_madgraph-pythia6-tauola/GEMCSCAna_DarkSUSY_ctau1000_0827/160828_001931/0000/"
treename1="GEMCSCAnalyzer/trk_eff_CSC_ME11"
x_bins = "(400,-0.005,0.005)"
"""
xaxis="(phi_lct_odd-phi_cscsh_odd)"#phi from momentum
xtitle="(phi_lct_odd-phi_cscsh_odd)"#phi from momentum

picname="GEMCSC_phi_ME11_odd"
cuts ="abs(phi_lct_odd-phi_cscsh_odd)<1"
st_title = "#phi resolution in ME11 odd chamber"
GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle,cuts,picname)
"""
xaxis = "(dphi_lct_odd)-(phi_cscsh_odd-phi_gemsh_odd)"
xtitle = "#Delta#Phi_{bending}"
picname="GEMCSC_deltaphi_ME11_odd"
cuts ="abs((dphi_lct_odd)-(phi_cscsh_odd-phi_gemsh_odd))<1"
st_title = "GEM-CSC bending angle resolution in ME11 odd chamber, p_{T}=2-50"
#GEMCSCResolution(filedir,treename1,xaxis,x_bins,st_title, xtitle,cuts,picname)
#GEMCSCResolution(filedir,treename1,"abs(dphi_lct_odd)",x_bins,"GEM-CSC bending angle distribution", "#Delta#Phi_{GEM,CSC}","abs(dphi_lct_odd)<.5","GEMCSC_pt40_deltaphi_ME11_odd")

#GEMCSCResolution(filedir,treename1,"ptphi_sh_even-phiM_gemcsc_even","(100,-0.05,0.05)","differences of direction from simhits and from gemcsc bending angle and position,ME11, even", "#Phi_{simhits}^{momentum}-#Phi_{gemcsc}^{momentum}","ptphi_sh_even>-9 && phiM_gemcsc_even>-9","phiM_Resolution_ME11_even")
#GEMCSCResolution(filedir,treename1,"ptphi_sh_odd-phiM_gemcsc_odd","(100,-0.05,0.05)","differences of direction from simhits and from gemcsc bending angle and position,ME11, odd", "#Phi_{simhits}^{momentum}-#Phi_{gemcsc}^{momentum}","ptphi_sh_odd>-9 && phiM_gemcsc_odd>-9","phiM_Resolution_ME11_odd")
#CSCMEs = ["ME11","ME12","ME13","ME21","ME22","ME31","ME32","ME41"]
#CSCMEs = ["ME1b","ME12","ME21","ME22"]
CSCMEs = ["ME1b", "ME21"]
CSCMEs = ["ME31","ME32","ME41"]
GEMs = ["GE11","GE21"]
evenodds = ["even","odd"]
#etas = [1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4]
etas = [1.2,1.4,1.6,1.8,2.0,2.2]
etas = [1.6,1.8,2.0,2.2,2.4]
bending = [-0.007,-0.005,-0.003,-0.001,0.001,0.003,0.005,0.007]
for m in range(len(CSCMEs)):
	for a in evenodds:
	    for neta in range(len(etas)-1):
	      #for nbend in range(len(bending)-1):
		me = CSCMEs[m]
		"""
		if (m==0 or m==2) and neta<2:
			continue
		if (m==1 or m==3) and neta>1:
			continue
		"""
		if m==0 or m==1:
			x_bins = "(100,-0.002,0.002)"
		if m==2 or m==3:
			x_bins = "(100,-0.003,0.003)"
		#ge = GEMs[m]
		x_bins = "(100,-1.0,1.0)"
		treename1 = "GEMCSCAnalyzer/trk_eff_CSC_%s"%(me)
		xaxis = "phi_layer3_fit_%s-phi_cscsh_%s"%(a,a)
		xaxis2 = "phi_lct_%s-phi_cscsh_%s"%(a,a)
		#xtitle = "dphi_lct_%s-dphi_sh_%s"%(a,a)
		#xtitle = "dphi_lct_%s-dphi_sh_%s"%(a,a)
		#xtitle = "#Delta#phi_{GEMCSC}^{L1}-#Delta#phi_{GEMCSC}^{SIM}"
		xtitle = "#phi_{CSC}^{L1}(fit)-#phi_{CSC}^{SIM}"
		xtitle2 = "#phi_{CSC}^{L1}(lct)-#phi_{CSC}^{SIM}"
		picname = "PhiRes_pt2_50_0827_compared/GEMCSC_ctau1000_pt5_50_CSCphi_fit_%s_%s_neta%d"%(me,a,neta)
		picname2 = "PhiRes_pt2_50_0827_compared/GEMCSC_ctau1000_pt5_50_CSCphi_lct_%s_%s_neta%d"%(me,a, neta)
		picname3 = "PhiRes_pt2_50_0827/GEMCSC_ctau1000_pt5_50_GMEphi_2strippad_%s_%s_neta%d"%(me,a, neta)
		picname4 = "PhiRes_pt2_50_0827/GEMCSC_ctau1000_pt5_50_GMEphi_4strippad_%s_%s_neta%d"%(me,a, neta)
		picname5 = "PhiRes_pt2_50_0827/GEMCSC_ctau1000_pt5_50_GMEphi_1strippad_%s_%s_neta%d"%(me,a, neta)
		st_title = "#phi_{CSC} resolution in %s chamber"%(me)
		#st_title2 = "#phi_{GEM} resolution in %s  chamber"%(ge)
		#text = "%.1f<|#eta|<%.1f, %.3f<#Delta#Phi^{SIM}<%.3f"%(etas[neta],etas[neta+1],bending[nbend],bending[nbend+1])
		text = "#splitline{%.1f<|#eta|<%.1f}{p_{T}>5 GeV, c#tau=1000}"%(etas[neta],etas[neta+1])
		#text = ""
		cuts = "eta<0 && pt>5 && fabs(phi_cscsh_%s)<4 &&  abs(phi_lct_%s)<4 && %f<fabs(eta) && fabs(eta)<%f"%(a,a,etas[neta],etas[neta+1])
		#cuts = "pt>5 && fabs(phi_cscsh_%s)<4 &&  abs(phi_lct_%s)<4 && %f<fabs(eta) && fabs(eta)<%f && fabs(phi_gemsh_%s)<4"%(a,a,etas[neta],etas[neta+1], a)
		#cuts = "fabs(dphi_sh_%s)>0&& abs(dphi_sh_%s)<1 && fabs(dphi_lct_%s)>0&& abs(dphi_lct_%s)<1"%(a,a,a,a)
		#cuts = "fabs(dphi_lct_%s)>0 &&fabs(dphi_lct_%s)<9&& dphi_sh_%s>%f && dphi_sh_%s<%f && %f<fabs(eta) && fabs(eta)<%f"%(a,a,a,bending[nbend],a,bending[nbend+1],etas[neta],etas[neta+1])
		GEMCSCResolution(filedir2_old,treename1,xaxis,x_bins,st_title, xtitle,cuts,text,picname+"_0826",True)
		GEMCSCResolution(filedir2,treename1,xaxis,x_bins,st_title, xtitle,cuts,text,picname+"_0828",True)
		#GEMCSCResolution(filedir1,treename1,xaxis2,x_bins,st_title, xtitle2,cuts,text,picname2,True)
		#GEMCSCResolution(filedir2,treename1,"phi_2strippad_%s-phi_gemsh_%s"%(a,a),x_bins,st_title2, "#phi_{GEMPad}(2strippad)-#phi_{GEM}^{SIM}",cuts,text,picname3,True)
		#GEMCSCResolution(filedir2,treename1,"phi_4strippad_%s-phi_gemsh_%s"%(a,a),x_bins,st_title2, "#phi_{GEMPad}(4strippad)-#phi_{GEM}^{SIM}",cuts,text,picname4,True)
		#GEMCSCResolution(filedir2,treename1,"phi_1strippad_%s-phi_gemsh_%s"%(a,a),x_bins,st_title2, "#phi_{GEMPad}(1strippad)-#phi_{GEM}^{SIM}",cuts,text,picname5,True)
