void ME41()
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

  TString canvName = "ME41_TP_2_";
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
  canv->SetTickx(0);
  canv->SetTicky(0);


    // Observed data

{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Aug  4 18:05:38 2014) by ROOT version5.34/07
 /*  TCanvas *c1 = new TCanvas("c1", "c1",11,52,700,500);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetFrameBorderMode(0);*/
   
   TH1F *b1 = new TH1F("b1","",32,1.7,2.5);
   b1->SetMinimum(0.0);
   b1->SetMaximum(1.10);
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
   //b1->GetYaxis()->SetLabelFont(42);
   //b1->GetYaxis()->SetLabelSize(0.035);
   //b1->GetYaxis()->SetTitleSize(0.05);
   //b1->GetYaxis()->SetTitleFont(42);
   b1->GetZaxis()->SetLabelFont(42);
   b1->GetZaxis()->SetLabelSize(0.035);
   b1->GetZaxis()->SetTitleSize(0.035);
   b1->GetZaxis()->SetTitleFont(42);
   b1->Draw("e3");
   Double_t xAxis1[28] = {1.75, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.45}; 
   
   TEfficiency * h1_clone1 = new TEfficiency("h1_clone","",27,xAxis1);
   h1_clone1->SetMarkerSize(0); 
   h1_clone1->SetConfidenceLevel(0.6826895);
   h1_clone1->SetBetaAlpha(1);
   h1_clone1->SetBetaBeta(1);
   h1_clone1->SetWeight(1);
   h1_clone1->SetStatisticOption(0);
   h1_clone1->SetPosteriorMode(0);
   h1_clone1->SetShortestInterval(0);


   h1_clone1->SetTotalEvents(0,17382);
   h1_clone1->SetPassedEvents(0,15345);
   h1_clone1->SetTotalEvents(1,198);
   h1_clone1->SetPassedEvents(1,146);
   h1_clone1->SetTotalEvents(2,519);
   h1_clone1->SetPassedEvents(2,442);
   h1_clone1->SetTotalEvents(3,750);
   h1_clone1->SetPassedEvents(3,691);
   h1_clone1->SetTotalEvents(4,772);
   h1_clone1->SetPassedEvents(4,706);
   h1_clone1->SetTotalEvents(5,788);
   h1_clone1->SetPassedEvents(5,726);
   h1_clone1->SetTotalEvents(6,723);
   h1_clone1->SetPassedEvents(6,670);
   h1_clone1->SetTotalEvents(7,722);
   h1_clone1->SetPassedEvents(7,567);
   h1_clone1->SetTotalEvents(8,701);
   h1_clone1->SetPassedEvents(8,505);
   h1_clone1->SetTotalEvents(9,821);
   h1_clone1->SetPassedEvents(9,725);
   h1_clone1->SetTotalEvents(10,769);
   h1_clone1->SetPassedEvents(10,727);
   h1_clone1->SetTotalEvents(11,766);
   h1_clone1->SetPassedEvents(11,710);
   h1_clone1->SetTotalEvents(12,775);
   h1_clone1->SetPassedEvents(12,703);
   h1_clone1->SetTotalEvents(13,747);
   h1_clone1->SetPassedEvents(13,695);
   h1_clone1->SetTotalEvents(14,777);
   h1_clone1->SetPassedEvents(14,683);
   h1_clone1->SetTotalEvents(15,681);
   h1_clone1->SetPassedEvents(15,510);
   h1_clone1->SetTotalEvents(16,705);
   h1_clone1->SetPassedEvents(16,588);
   h1_clone1->SetTotalEvents(17,783);
   h1_clone1->SetPassedEvents(17,708);
   h1_clone1->SetTotalEvents(18,772);
   h1_clone1->SetPassedEvents(18,725);
   h1_clone1->SetTotalEvents(19,765);
   h1_clone1->SetPassedEvents(19,704);
   h1_clone1->SetTotalEvents(20,765);
   h1_clone1->SetPassedEvents(20,700);
   h1_clone1->SetTotalEvents(21,801);
   h1_clone1->SetPassedEvents(21,751);
   h1_clone1->SetTotalEvents(22,766);
   h1_clone1->SetPassedEvents(22,701);
   h1_clone1->SetTotalEvents(23,739);
   h1_clone1->SetPassedEvents(23,664);
   h1_clone1->SetTotalEvents(24,560);
   h1_clone1->SetPassedEvents(24,455);
   h1_clone1->SetTotalEvents(25,174);
   h1_clone1->SetPassedEvents(25,117);
   h1_clone1->SetTotalEvents(26,2);
   h1_clone1->SetPassedEvents(26,2);
   h1_clone1->SetTotalEvents(27,0);
   h1_clone1->SetPassedEvents(27,0);
   h1_clone1->SetTotalEvents(28,0);
   h1_clone1->SetPassedEvents(28,0);


   ci = TColor::GetColor("#990000");
   h1_clone1->SetFillColor(ci);
   h1_clone1->Draw("e3samep");
   Double_t xAxis2[28] = {1.75, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.45}; 
   
   TEfficiency * h1_clone2 = new TEfficiency("h1_clone","",27,xAxis2);
   
   h1_clone2->SetMarkerSize(0);
   h1_clone2->SetConfidenceLevel(0.6826895);
   h1_clone2->SetBetaAlpha(1);
   h1_clone2->SetBetaBeta(1);
   h1_clone2->SetWeight(1);
   h1_clone2->SetStatisticOption(0);
   h1_clone2->SetPosteriorMode(0);
   h1_clone2->SetShortestInterval(0);

   h1_clone2->SetTotalEvents(0,36080);
   h1_clone2->SetPassedEvents(0,32171);
   h1_clone2->SetTotalEvents(1,373);
   h1_clone2->SetPassedEvents(1,307);
   h1_clone2->SetTotalEvents(2,1198);
   h1_clone2->SetPassedEvents(2,1105);
   h1_clone2->SetTotalEvents(3,1617);
   h1_clone2->SetPassedEvents(3,1541);
   h1_clone2->SetTotalEvents(4,1593);
   h1_clone2->SetPassedEvents(4,1547);
   h1_clone2->SetTotalEvents(5,1586);
   h1_clone2->SetPassedEvents(5,1575);
   h1_clone2->SetTotalEvents(6,1488);
   h1_clone2->SetPassedEvents(6,1478);
   h1_clone2->SetTotalEvents(7,1608);
   h1_clone2->SetPassedEvents(7,1527);
   h1_clone2->SetTotalEvents(8,1541);
   h1_clone2->SetPassedEvents(8,1408);
   h1_clone2->SetTotalEvents(9,1557);
   h1_clone2->SetPassedEvents(9,1517);
   h1_clone2->SetTotalEvents(10,1636);
   h1_clone2->SetPassedEvents(10,1629);
   h1_clone2->SetTotalEvents(11,1618);
   h1_clone2->SetPassedEvents(11,1608);
   h1_clone2->SetTotalEvents(12,1612);
   h1_clone2->SetPassedEvents(12,1603);
   h1_clone2->SetTotalEvents(13,1548);
   h1_clone2->SetPassedEvents(13,1542);
   h1_clone2->SetTotalEvents(14,1624);
   h1_clone2->SetPassedEvents(14,1607);
   h1_clone2->SetTotalEvents(15,1370);
   h1_clone2->SetPassedEvents(15,1329);
   h1_clone2->SetTotalEvents(16,1451);
   h1_clone2->SetPassedEvents(16,1418);
   h1_clone2->SetTotalEvents(17,1581);
   h1_clone2->SetPassedEvents(17,1574);
   h1_clone2->SetTotalEvents(18,1576);
   h1_clone2->SetPassedEvents(18,1570);
   h1_clone2->SetTotalEvents(19,1576);
   h1_clone2->SetPassedEvents(19,1570);
   h1_clone2->SetTotalEvents(20,1671);
   h1_clone2->SetPassedEvents(20,1664);
   h1_clone2->SetTotalEvents(21,1519);
   h1_clone2->SetPassedEvents(21,1511);
   h1_clone2->SetTotalEvents(22,1593);
   h1_clone2->SetPassedEvents(22,1578);
   h1_clone2->SetTotalEvents(23,1531);
   h1_clone2->SetPassedEvents(23,1525);
   h1_clone2->SetTotalEvents(24,1155);
   h1_clone2->SetPassedEvents(24,1137);
   h1_clone2->SetTotalEvents(25,381);
   h1_clone2->SetPassedEvents(25,372);
   h1_clone2->SetTotalEvents(26,3);
   h1_clone2->SetPassedEvents(26,3);
   h1_clone2->SetTotalEvents(27,0);
   h1_clone2->SetPassedEvents(27,0);
   h1_clone2->SetTotalEvents(28,0);
   h1_clone2->SetPassedEvents(28,0);



   ci = TColor::GetColor("#9999ff");
   h1_clone2->SetFillColor(ci);
   h1_clone2->Draw("e3samep");
   
   TLegend *leg = new TLegend(0.2,0.15,.89,.32,"","brNDC");
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
   tex->DrawLatex(0.62,0.60,"p_{T}^{sim} > 10 GeV");
   tex->DrawLatex(0.1,0.94,"YE4/1 stub reconstruction");
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

  TFile *f = new TFile("ME41.root", "RECREATE");
  canv->Write();
  f->Close();

  return canv;
}

