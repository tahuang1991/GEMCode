void ME31()
{
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  // example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
//  example_plot( iPeriod, 11 );  // left-aligned
//  example_plot( iPeriod, 33 );  // right-aligned

  writeExtraText = false;       // remove Preliminary
  
  // example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

//  example_plot( iPeriod, 11 );  // default: left-aligned
//  example_plot( iPeriod, 22 );  // centered
  example_plot( iPeriod, 33 );  // right-aligned


  
}

TCanvas* example_plot( int iPeriod, int iPos )
{ 
  int W = 800;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TString canvName = "ME31_TP_2_";
  canvName += W;
  canvName += "-";
  canvName += H;
  canvName += "_";  
  canvName += iPeriod;
  if( writeExtraText ) canvName += "-prelim";
  if( iPos%10==0 ) canvName += "-out";
  else if( iPos%10==1 ) canvName += "-left";
  else if( iPos%10==2 )  canvName += "-center";
  else if( iPos%10==3 )  canvName += "-right";

  TCanvas* canv = new TCanvas(canvName,canvName,10,10,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(1);
  canv->SetTicky(1);



  {
    // Observed data
//=========Macro generated from canvas: c1/c1
//=========  (Mon Aug  4 17:10:03 2014) by ROOT version5.34/07
  /* TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   gStyle->SetOptStat(0);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.126);
   c1->SetRightMargin(0.04);
   c1->SetTopMargin(0.06);
   c1->SetBottomMargin(0.13);
   c1->SetFrameBorderMode(0);*/
   
   TH1F *b1 = new TH1F("b1","",36,1.6,2.5);
   b1->SetMinimum(0.0);
   b1->SetMaximum(1.1);
   b1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   b1->SetLineColor(ci);
   b1->GetXaxis()->SetTitle("Simulated muon #eta");
   //b1->GetXaxis()->SetLabelFont(42);
   //b1->GetXaxis()->SetLabelSize(0.035);
   //b1->GetXaxis()->SetTitleSize(0.05);
   //b1->GetXaxis()->SetTitleFont(42);
   b1->GetYaxis()->SetTitle("Efficiency");
   b1->GetYaxis()->SetTitleOffset(0.8);
   //b1->GetYaxis()->SetNdivisions(520);
   //b1->GetYaxis()->SetLabelFont(42);
   //b1->GetYaxis()->SetLabelSize(0.035);
   //b1->GetYaxis()->SetTitleSize(0.05);
   //b1->GetYaxis()->SetTitleFont(42);
   //b1->GetZaxis()->SetLabelFont(42);
   //b1->GetZaxis()->SetLabelSize(0.035);
   //b1->GetZaxis()->SetTitleSize(0.035);
   //b1->GetZaxis()->SetTitleFont(42);
   b1->Draw();
   Double_t xAxis1[32] = {1.65, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.45}; 
   
   TEfficiency * h1_clone1 = new TEfficiency("h1_clone","",31,xAxis1);
   h1_clone1->SetMarkerSize(0);

   h1_clone1->SetConfidenceLevel(0.6826895);
   h1_clone1->SetBetaAlpha(1);
   h1_clone1->SetBetaBeta(1);
   h1_clone1->SetWeight(1);
   h1_clone1->SetStatisticOption(0);
   h1_clone1->SetPosteriorMode(0);
   h1_clone1->SetShortestInterval(0);
   h1_clone1->SetTotalEvents(0,20698);
   h1_clone1->SetPassedEvents(0,18153);
   h1_clone1->SetTotalEvents(1,49);
   h1_clone1->SetPassedEvents(1,27);
   h1_clone1->SetTotalEvents(2,381);
   h1_clone1->SetPassedEvents(2,310);
   h1_clone1->SetTotalEvents(3,688);
   h1_clone1->SetPassedEvents(3,607);
   h1_clone1->SetTotalEvents(4,743);
   h1_clone1->SetPassedEvents(4,695);
   h1_clone1->SetTotalEvents(5,767);
   h1_clone1->SetPassedEvents(5,718);
   h1_clone1->SetTotalEvents(6,717);
   h1_clone1->SetPassedEvents(6,668);
   h1_clone1->SetTotalEvents(7,755);
   h1_clone1->SetPassedEvents(7,718);
   h1_clone1->SetTotalEvents(8,772);
   h1_clone1->SetPassedEvents(8,708);
   h1_clone1->SetTotalEvents(9,788);
   h1_clone1->SetPassedEvents(9,569);
   h1_clone1->SetTotalEvents(10,722);
   h1_clone1->SetPassedEvents(10,497);
   h1_clone1->SetTotalEvents(11,734);
   h1_clone1->SetPassedEvents(11,642);
   h1_clone1->SetTotalEvents(12,725);
   h1_clone1->SetPassedEvents(12,659);
   h1_clone1->SetTotalEvents(13,829);
   h1_clone1->SetPassedEvents(13,768);
   h1_clone1->SetTotalEvents(14,769);
   h1_clone1->SetPassedEvents(14,716);
   h1_clone1->SetTotalEvents(15,766);
   h1_clone1->SetPassedEvents(15,701);
   h1_clone1->SetTotalEvents(16,775);
   h1_clone1->SetPassedEvents(16,715);
   h1_clone1->SetTotalEvents(17,747);
   h1_clone1->SetPassedEvents(17,682);
   h1_clone1->SetTotalEvents(18,756);
   h1_clone1->SetPassedEvents(18,641);
   h1_clone1->SetTotalEvents(19,695);
   h1_clone1->SetPassedEvents(19,474);
   h1_clone1->SetTotalEvents(20,711);
   h1_clone1->SetPassedEvents(20,584);
   h1_clone1->SetTotalEvents(21,783);
   h1_clone1->SetPassedEvents(21,705);
   h1_clone1->SetTotalEvents(22,772);
   h1_clone1->SetPassedEvents(22,700);
   h1_clone1->SetTotalEvents(23,765);
   h1_clone1->SetPassedEvents(23,696);
   h1_clone1->SetTotalEvents(24,765);
   h1_clone1->SetPassedEvents(24,688);
   h1_clone1->SetTotalEvents(25,801);
   h1_clone1->SetPassedEvents(25,733);
   h1_clone1->SetTotalEvents(26,766);
   h1_clone1->SetPassedEvents(26,684);
   h1_clone1->SetTotalEvents(27,746);
   h1_clone1->SetPassedEvents(27,668);
   h1_clone1->SetTotalEvents(28,757);
   h1_clone1->SetPassedEvents(28,655);
   h1_clone1->SetTotalEvents(29,499);
   h1_clone1->SetPassedEvents(29,386);
   h1_clone1->SetTotalEvents(30,160);
   h1_clone1->SetPassedEvents(30,89);
   h1_clone1->SetTotalEvents(31,0);
   h1_clone1->SetPassedEvents(31,0);
   h1_clone1->SetTotalEvents(32,4);
   h1_clone1->SetPassedEvents(32,3);


   ci = TColor::GetColor("#990000");
   h1_clone1->SetFillColor(ci);
   h1_clone1->Draw("e3samep");
   Double_t xAxis2[32] = {1.65, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.45}; 
   
   TEfficiency * h1_clone2 = new TEfficiency("h1_clone","",31,xAxis2);
   h1_clone2->SetMarkerSize(0); 
  
   h1_clone2->SetConfidenceLevel(0.6826895);
   h1_clone2->SetBetaAlpha(1);
   h1_clone2->SetBetaBeta(1);
   h1_clone2->SetWeight(1);
   h1_clone2->SetStatisticOption(0);
   h1_clone2->SetPosteriorMode(0);
   h1_clone2->SetShortestInterval(0);
   h1_clone2->SetTotalEvents(0,42875);
   h1_clone2->SetPassedEvents(0,38320);
   h1_clone2->SetTotalEvents(1,134);
   h1_clone2->SetPassedEvents(1,103);
   h1_clone2->SetTotalEvents(2,765);
   h1_clone2->SetPassedEvents(2,678);
   h1_clone2->SetTotalEvents(3,1419);
   h1_clone2->SetPassedEvents(3,1335);
   h1_clone2->SetTotalEvents(4,1602);
   h1_clone2->SetPassedEvents(4,1519);
   h1_clone2->SetTotalEvents(5,1524);
   h1_clone2->SetPassedEvents(5,1466);
   h1_clone2->SetTotalEvents(6,1583);
   h1_clone2->SetPassedEvents(6,1575);
   h1_clone2->SetTotalEvents(7,1627);
   h1_clone2->SetPassedEvents(7,1624);
   h1_clone2->SetTotalEvents(8,1593);
   h1_clone2->SetPassedEvents(8,1588);
   h1_clone2->SetTotalEvents(9,1585);
   h1_clone2->SetPassedEvents(9,1548);
   h1_clone2->SetTotalEvents(10,1483);
   h1_clone2->SetPassedEvents(10,1445);
   h1_clone2->SetTotalEvents(11,1626);
   h1_clone2->SetPassedEvents(11,1616);
   h1_clone2->SetTotalEvents(12,1578);
   h1_clone2->SetPassedEvents(12,1572);
   h1_clone2->SetTotalEvents(13,1571);
   h1_clone2->SetPassedEvents(13,1566);
   h1_clone2->SetTotalEvents(14,1636);
   h1_clone2->SetPassedEvents(14,1622);
   h1_clone2->SetTotalEvents(15,1618);
   h1_clone2->SetPassedEvents(15,1605);
   h1_clone2->SetTotalEvents(16,1612);
   h1_clone2->SetPassedEvents(16,1608);
   h1_clone2->SetTotalEvents(17,1546);
   h1_clone2->SetPassedEvents(17,1542);
   h1_clone2->SetTotalEvents(18,1590);
   h1_clone2->SetPassedEvents(18,1569);
   h1_clone2->SetTotalEvents(19,1371);
   h1_clone2->SetPassedEvents(19,1327);
   h1_clone2->SetTotalEvents(20,1467);
   h1_clone2->SetPassedEvents(20,1449);
   h1_clone2->SetTotalEvents(21,1581);
   h1_clone2->SetPassedEvents(21,1572);
   h1_clone2->SetTotalEvents(22,1576);
   h1_clone2->SetPassedEvents(22,1570);
   h1_clone2->SetTotalEvents(23,1576);
   h1_clone2->SetPassedEvents(23,1572);
   h1_clone2->SetTotalEvents(24,1671);
   h1_clone2->SetPassedEvents(24,1668);
   h1_clone2->SetTotalEvents(25,1519);
   h1_clone2->SetPassedEvents(25,1503);
   h1_clone2->SetTotalEvents(26,1593);
   h1_clone2->SetPassedEvents(26,1576);
   h1_clone2->SetTotalEvents(27,1544);
   h1_clone2->SetPassedEvents(27,1536);
   h1_clone2->SetTotalEvents(28,1503);
   h1_clone2->SetPassedEvents(28,1490);
   h1_clone2->SetTotalEvents(29,1094);
   h1_clone2->SetPassedEvents(29,1072);
   h1_clone2->SetTotalEvents(30,341);
   h1_clone2->SetPassedEvents(30,337);
   h1_clone2->SetTotalEvents(31,0);
   h1_clone2->SetPassedEvents(31,0);
   h1_clone2->SetTotalEvents(32,9);
   h1_clone2->SetPassedEvents(32,8);


   
   ci = TColor::GetColor("#9999ff");
   h1_clone2->SetFillColor(ci);
   h1_clone2->Draw("e3samep");
   
   TLegend *leg = new TLegend(0.22,0.15,.89,.32,"","brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry(h1_clone1,"Phase I detector : CSC Local Trigger","f");
   entry=leg->AddEntry(h1_clone2,"CSC+RPC Integrated Local Trigger","f");
   leg->Draw();
   TLatex *   tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.05);
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->DrawLatex(0.657,0.56,"p_{T}^{sim} > 10 GeV");
   tex->DrawLatex(0.1,0.94,"YE3/1 stub reconstruction");
   
  /* c1->Modified();
   c1->cd();
   c1->SetSelected(c1);*/
}
  // writing the lumi information and the CMS "logo"
  {
    CMS_lumi( canv, iPeriod, iPos );
  }

  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();

  canv->Print(canvName+"_0_11.pdf",".pdf");
  canv->Print(canvName+"_0_11.png",".png");
  canv->Print(canvName+"_0_11.eps",".eps");

  TFile *f = new TFile("ME31.root", "RECREATE");
  canv->Write();
  f->Close();
  return canv;
}

