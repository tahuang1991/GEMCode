void ME11()
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

  TString canvName = "ME11_TP_2_";
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
//=========  (Mon Aug  4 18:04:36 2014) by ROOT version5.34/07
/*   TCanvas *c1 = new TCanvas("c1", "c1",11,52,700,500);
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
   
   TH1F *b1 = new TH1F("b1","",40,1.5,2.5);
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
   //b1->GetYaxis()->SetLabelFont(42);
   //b1->GetYaxis()->SetLabelSize(0.035);
   //b1->GetYaxis()->SetTitleSize(0.05);
   //b1->GetYaxis()->SetTitleFont(42);
   b1->GetZaxis()->SetLabelFont(42);
   b1->GetZaxis()->SetLabelSize(0.035);
   b1->GetZaxis()->SetTitleSize(0.035);
   b1->GetZaxis()->SetTitleFont(42);
   b1->Draw("");
   Double_t xAxis1[37] = {1.5, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.5}; 
   
   TEfficiency * h1_clone1 = new TEfficiency("h1_clone","",36,xAxis1);
   
   h1_clone1->SetMarkerSize(0);
   h1_clone1->SetConfidenceLevel(0.6826895);
   h1_clone1->SetBetaAlpha(1);
   h1_clone1->SetBetaBeta(1);
   h1_clone1->SetWeight(1);
   h1_clone1->SetStatisticOption(0);
   h1_clone1->SetPosteriorMode(0);
   h1_clone1->SetShortestInterval(0);
   
   h1_clone1->SetTotalEvents(0,0);
   h1_clone1->SetPassedEvents(0,0);
   h1_clone1->SetTotalEvents(1,1204);
   h1_clone1->SetPassedEvents(1,1017);
   h1_clone1->SetTotalEvents(2,1142);
   h1_clone1->SetPassedEvents(2,1027);
   h1_clone1->SetTotalEvents(3,1447);
   h1_clone1->SetPassedEvents(3,1418);
   h1_clone1->SetTotalEvents(4,1584);
   h1_clone1->SetPassedEvents(4,1555);
   h1_clone1->SetTotalEvents(5,1544);
   h1_clone1->SetPassedEvents(5,1507);
   h1_clone1->SetTotalEvents(6,1544);
   h1_clone1->SetPassedEvents(6,1508);
   h1_clone1->SetTotalEvents(7,1564);
   h1_clone1->SetPassedEvents(7,1531);
   h1_clone1->SetTotalEvents(8,1532);
   h1_clone1->SetPassedEvents(8,1503);
   h1_clone1->SetTotalEvents(9,1492);
   h1_clone1->SetPassedEvents(9,1469);
   h1_clone1->SetTotalEvents(10,1558);
   h1_clone1->SetPassedEvents(10,1521);
   h1_clone1->SetTotalEvents(11,1454);
   h1_clone1->SetPassedEvents(11,1420);
   h1_clone1->SetTotalEvents(12,1514);
   h1_clone1->SetPassedEvents(12,1458);
   h1_clone1->SetTotalEvents(13,1576);
   h1_clone1->SetPassedEvents(13,1524);
   h1_clone1->SetTotalEvents(14,1584);
   h1_clone1->SetPassedEvents(14,1535);
   h1_clone1->SetTotalEvents(15,1456);
   h1_clone1->SetPassedEvents(15,1404);
   h1_clone1->SetTotalEvents(16,1488);
   h1_clone1->SetPassedEvents(16,1433);
   h1_clone1->SetTotalEvents(17,1460);
   h1_clone1->SetPassedEvents(17,1415);
   h1_clone1->SetTotalEvents(18,1678);
   h1_clone1->SetPassedEvents(18,1624);
   h1_clone1->SetTotalEvents(19,1554);
   h1_clone1->SetPassedEvents(19,1497);
   h1_clone1->SetTotalEvents(20,1556);
   h1_clone1->SetPassedEvents(20,1383);
   h1_clone1->SetTotalEvents(21,1566);
   h1_clone1->SetPassedEvents(21,1430);
   h1_clone1->SetTotalEvents(22,1514);
   h1_clone1->SetPassedEvents(22,1324);
   h1_clone1->SetTotalEvents(23,1594);
   h1_clone1->SetPassedEvents(23,1460);
   h1_clone1->SetTotalEvents(24,1568);
   h1_clone1->SetPassedEvents(24,1500);
   h1_clone1->SetTotalEvents(25,1583);
   h1_clone1->SetPassedEvents(25,1519);
   h1_clone1->SetTotalEvents(26,1586);
   h1_clone1->SetPassedEvents(26,1522);
   h1_clone1->SetTotalEvents(27,1550);
   h1_clone1->SetPassedEvents(27,1486);
   h1_clone1->SetTotalEvents(28,1553);
   h1_clone1->SetPassedEvents(28,1480);
   h1_clone1->SetTotalEvents(29,1544);
   h1_clone1->SetPassedEvents(29,1463);
   h1_clone1->SetTotalEvents(30,1620);
   h1_clone1->SetPassedEvents(30,1538);
   h1_clone1->SetTotalEvents(31,1542);
   h1_clone1->SetPassedEvents(31,1460);
   h1_clone1->SetTotalEvents(32,1506);
   h1_clone1->SetPassedEvents(32,1405);
   h1_clone1->SetTotalEvents(33,1493);
   h1_clone1->SetPassedEvents(33,1319);
   h1_clone1->SetTotalEvents(34,1009);
   h1_clone1->SetPassedEvents(34,858);
   h1_clone1->SetTotalEvents(35,661);
   h1_clone1->SetPassedEvents(35,526);
   h1_clone1->SetTotalEvents(36,205);
   h1_clone1->SetPassedEvents(36,118);
   h1_clone1->SetTotalEvents(37,0);
   h1_clone1->SetPassedEvents(37,0);



   ci = TColor::GetColor("#990000");
   h1_clone1->SetFillColor(ci);
   h1_clone1->Draw("e3samep");
   Double_t xAxis2[37] = {1.5, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2, 2.025, 2.05, 2.075, 2.1, 2.125, 2.15, 2.175, 2.2, 2.225, 2.25, 2.275, 2.3, 2.325, 2.35, 2.375, 2.4, 2.425, 2.45, 2.5}; 
   
   TEfficiency * h1_clone2 = new TEfficiency("h1_clone","",36,xAxis2);
   
   h1_clone2->SetMarkerSize(0);
   h1_clone2->SetConfidenceLevel(0.6826895);
   h1_clone2->SetBetaAlpha(1);
   h1_clone2->SetBetaBeta(1);
   h1_clone2->SetWeight(1);
   h1_clone2->SetStatisticOption(0);
   h1_clone2->SetPosteriorMode(0);
   h1_clone2->SetShortestInterval(0);
   
   h1_clone2->SetTotalEvents(0,0);
   h1_clone2->SetPassedEvents(0,0);
   h1_clone2->SetTotalEvents(1,2286);
   h1_clone2->SetPassedEvents(1,1984);
   h1_clone2->SetTotalEvents(2,2460);
   h1_clone2->SetPassedEvents(2,2256);
   h1_clone2->SetTotalEvents(3,3134);
   h1_clone2->SetPassedEvents(3,3083);
   h1_clone2->SetTotalEvents(4,3138);
   h1_clone2->SetPassedEvents(4,3121);
   h1_clone2->SetTotalEvents(5,3132);
   h1_clone2->SetPassedEvents(5,3116);
   h1_clone2->SetTotalEvents(6,3092);
   h1_clone2->SetPassedEvents(6,3072);
   h1_clone2->SetTotalEvents(7,3074);
   h1_clone2->SetPassedEvents(7,3041);
   h1_clone2->SetTotalEvents(8,3098);
   h1_clone2->SetPassedEvents(8,3081);
   h1_clone2->SetTotalEvents(9,3204);
   h1_clone2->SetPassedEvents(9,3169);
   h1_clone2->SetTotalEvents(10,3048);
   h1_clone2->SetPassedEvents(10,3025);
   h1_clone2->SetTotalEvents(11,3168);
   h1_clone2->SetPassedEvents(11,3122);
   h1_clone2->SetTotalEvents(12,3254);
   h1_clone2->SetPassedEvents(12,3223);
   h1_clone2->SetTotalEvents(13,3186);
   h1_clone2->SetPassedEvents(13,3161);
   h1_clone2->SetTotalEvents(14,3172);
   h1_clone2->SetPassedEvents(14,3137);
   h1_clone2->SetTotalEvents(15,2976);
   h1_clone2->SetPassedEvents(15,2943);
   h1_clone2->SetTotalEvents(16,3256);
   h1_clone2->SetPassedEvents(16,3219);
   h1_clone2->SetTotalEvents(17,3156);
   h1_clone2->SetPassedEvents(17,3114);
   h1_clone2->SetTotalEvents(18,3142);
   h1_clone2->SetPassedEvents(18,3109);
   h1_clone2->SetTotalEvents(19,3272);
   h1_clone2->SetPassedEvents(19,3222);
   h1_clone2->SetTotalEvents(20,3236);
   h1_clone2->SetPassedEvents(20,3160);
   h1_clone2->SetTotalEvents(21,3224);
   h1_clone2->SetPassedEvents(21,3139);
   h1_clone2->SetTotalEvents(22,3096);
   h1_clone2->SetPassedEvents(22,3025);
   h1_clone2->SetTotalEvents(23,3321);
   h1_clone2->SetPassedEvents(23,3251);
   h1_clone2->SetTotalEvents(24,3104);
   h1_clone2->SetPassedEvents(24,3040);
   h1_clone2->SetTotalEvents(25,3204);
   h1_clone2->SetPassedEvents(25,3147);
   h1_clone2->SetTotalEvents(26,3168);
   h1_clone2->SetPassedEvents(26,3109);
   h1_clone2->SetTotalEvents(27,3152);
   h1_clone2->SetPassedEvents(27,3089);
   h1_clone2->SetTotalEvents(28,3152);
   h1_clone2->SetPassedEvents(28,3063);
   h1_clone2->SetTotalEvents(29,3342);
   h1_clone2->SetPassedEvents(29,3258);
   h1_clone2->SetTotalEvents(30,3036);
   h1_clone2->SetPassedEvents(30,2942);
   h1_clone2->SetTotalEvents(31,3185);
   h1_clone2->SetPassedEvents(31,3089);
   h1_clone2->SetTotalEvents(32,3087);
   h1_clone2->SetPassedEvents(32,2989);
   h1_clone2->SetTotalEvents(33,2920);
   h1_clone2->SetPassedEvents(33,2713);
   h1_clone2->SetTotalEvents(34,2067);
   h1_clone2->SetPassedEvents(34,1834);
   h1_clone2->SetTotalEvents(35,1384);
   h1_clone2->SetPassedEvents(35,1201);
   h1_clone2->SetTotalEvents(36,409);
   h1_clone2->SetPassedEvents(36,248);
   h1_clone2->SetTotalEvents(37,0);
   h1_clone2->SetPassedEvents(37,0);


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
   entry=leg->AddEntry(h1_clone2,"GEM+CSC Integrated Local Trigger","f");
   leg->Draw();
   TLatex *   tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.05);
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->DrawLatex(0.65,0.65,"p_{T}^{sim} > 10 GeV");
   tex->DrawLatex(0.1,0.94,"YE1/1 stub reconstruction");
   
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
  
  TFile *f = new TFile("ME11.root", "RECREATE");
  canv->Write();
  f->Close();

  return canv;
}
