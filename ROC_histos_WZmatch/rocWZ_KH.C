#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <TChain.h>
#include "TMath.h"
#include "TGraphErrors.h"

void rocWZ_KH(){
  TFile * file_inputdistributions1;
  TFile * file_inputdistributions2;
  TFile * file_output;
  TString filename = "ROCplot_WJ1_TotalBkg";     //file name to be saved
  TString name;
  
  file_inputdistributions1 = new TFile("TChiWZ_600_100_MC2018.root", "READ");
  file_inputdistributions2 = new TFile("Total_MC2018.root", "READ");
  // based on J1 
  ///* 
  TH1F * h_sig = (TH1F*) file_inputdistributions1->Get("wAK8J1Tau21");
  TH1F * h_bkg = (TH1F*) file_inputdistributions2->Get("wAK8J1Tau21");
  TH1F * h_sigZDisc = (TH1F*) file_inputdistributions1->Get("wAK8J1zDis");
  TH1F * h_bkgZDisc = (TH1F*) file_inputdistributions2->Get("wAK8J1zDis");
  // TH1F * h_sigZHDisc = (TH1F*) file_inputdistributions1->Get("AK8J1zhDisMD");
  // TH1F * h_bkgZHDisc = (TH1F*) file_inputdistributions2->Get("AK8J1zhDisMD");
  TH1F * h_sigdoubleBDisc = (TH1F*) file_inputdistributions1->Get("wAK8J1wDis");
  TH1F * h_bkgdoubleBDisc = (TH1F*) file_inputdistributions2->Get("wAK8J1wDis");
  TH1F * h_sigdeepdoubleBDisc = (TH1F*) file_inputdistributions1->Get("wAK8J1wDisDC");
  TH1F * h_bkgdeepdoubleBDisc = (TH1F*) file_inputdistributions2->Get("wAK8J1wDisDC");
  //*/
  

  
  // based on J2
  /*
  TH1F * h_sig = (TH1F*) file_inputdistributions1->Get("AK8J2Tau21");
  TH1F * h_bkg = (TH1F*) file_inputdistributions2->Get("AK8J2Tau21");
  TH1F * h_sigZDisc = (TH1F*) file_inputdistributions1->Get("AK8J2zDis");
  TH1F * h_bkgZDisc = (TH1F*) file_inputdistributions2->Get("AK8J2zDis");
  TH1F * h_sigdoubleBDisc = (TH1F*) file_inputdistributions1->Get("AK8J2doubleBDis");
  TH1F * h_bkgdoubleBDisc = (TH1F*) file_inputdistributions2->Get("AK8J2doubleBDis");
  TH1F * h_sigdeepdoubleBDisc = (TH1F*) file_inputdistributions1->Get("AK8J2deepdoubleBDis");
  TH1F * h_bkgdeepdoubleBDisc = (TH1F*) file_inputdistributions2->Get("AK8J2deepdoubleBDis");
  */
  cout<<h_sig->GetNbinsX()<<"\t"<<h_bkg->GetNbinsX()<<endl;
  cout<<h_sigZDisc->GetNbinsX()<<"\t"<<h_bkgZDisc->GetNbinsX()<<endl;
  cout<<"********"<<endl;

  Double_t s[101],s_ZDisc[101],s_ZHDisc[101],s_doubleBDisc[101],s_deepdoubleBDisc[101] ={0};
  Double_t b[101],b_ZDisc[101],b_ZHDisc[101],b_doubleBDisc[101],b_deepdoubleBDisc[101] ={0};
  Double_t s_total,s_zDiscTotal=0,s_zHDiscTotal=0, s_doubleBTotal=0, s_deepdoubleBTotal=0;
  Double_t b_total,b_zDiscTotal=0,b_zHDiscTotal=0, b_doubleBTotal=0, b_deepdoubleBTotal=0;
  Double_t s1[101],s1_ZDisc[101],s1_ZHDisc[101],s1_doubleBDisc[101],s1_deepdoubleBDisc[101] = {0};
  Double_t b1[101],b1_ZDisc[101],b1_ZHDisc[101],b1_doubleBDisc[101],b1_deepdoubleBDisc[101] = {0};


  TCanvas *c1 = new TCanvas("c1","Individual plot", 200,10,700,500);
  c1->cd();
  h_sig->Draw();
  h_bkg->Draw("sames");
  
  s_total=  h_sig->Integral(0, 100);
  b_total=  h_bkg->Integral(0, 100);
  s_zDiscTotal=  h_sigZDisc->Integral(0, 100);
  b_zDiscTotal=  h_bkgZDisc->Integral(0, 100);
  // s_zHDiscTotal=  h_sigZHDisc->Integral(0, 100);
  // b_zHDiscTotal=  h_bkgZHDisc->Integral(0, 100);
  s_doubleBTotal = h_sigdoubleBDisc->Integral(0, 100);
  b_doubleBTotal = h_bkgdoubleBDisc->Integral(0, 100);
  s_deepdoubleBTotal = h_sigdeepdoubleBDisc->Integral(0, 100);
  b_deepdoubleBTotal = h_bkgdeepdoubleBDisc->Integral(0, 100);

  cout<<"sig total :"<<s_total<<endl;
  cout<<"bkg total :"<<b_total<<endl;
  int i;
  int j;
  for( i= 0 , j = 100; i <= 100; i++,j--){
    s[i] = h_sig->Integral(0, i);
    b[i] = h_bkg->Integral(0, i);
    s_ZDisc[j] = h_sigZDisc->Integral(j, 100);
    b_ZDisc[j] = h_bkgZDisc->Integral(j, 100);
    // s_ZHDisc[j] = h_sigZHDisc->Integral(j, 100);
    // b_ZHDisc[j] = h_bkgZHDisc->Integral(j, 100);
    s_doubleBDisc[j] = h_sigdoubleBDisc->Integral(j, 100);
    b_doubleBDisc[j] = h_bkgdoubleBDisc->Integral(j, 100);
    s_deepdoubleBDisc[j] = h_sigdeepdoubleBDisc->Integral(j, 100);
    b_deepdoubleBDisc[j] = h_bkgdeepdoubleBDisc->Integral(j, 100);

    cout<<"sig bin cont in "<<i<<" is "<<s[i]<<endl;
    cout<<"bkg bin cont in "<<i<<" is "<<b[i]<<endl;
    cout<<"***********************************"<<endl;
    cout<<"sigdeepDisc bin cont in "<<j<<" is "<<s_ZDisc[j]<<endl;
    cout<<"bkgdeepDisc bin cont in "<<j<<" is "<<b_ZDisc[j]<<endl;
    cout<<"**********Ratios********"<<endl;  


    s1[i] = s[i]/s_total;
    b1[i] = b[i]/b_total;
    s1_ZDisc[i] = s_ZDisc[j]/s_zDiscTotal;
    b1_ZDisc[i] = b_ZDisc[j]/b_zDiscTotal;
    // s1_ZHDisc[i] = s_ZHDisc[j]/s_zHDiscTotal;
    // b1_ZHDisc[i] = b_ZHDisc[j]/b_zHDiscTotal;
    s1_doubleBDisc[i] = s_doubleBDisc[j]/s_doubleBTotal;
    b1_doubleBDisc[i] = b_doubleBDisc[j]/b_doubleBTotal;
    s1_deepdoubleBDisc[i] = s_deepdoubleBDisc[j]/s_deepdoubleBTotal;
    b1_deepdoubleBDisc[i] = b_deepdoubleBDisc[j]/b_deepdoubleBTotal;

    //when j=30 doubleB =0.3, i=70
    if (i==35){cout<<"for tau21 = 0.35 X: "<< b1[i] <<" Y "<< s1[i]<<endl;}
    if (i==8){cout<<"for Wdis = 0.92 X: " << b1_doubleBDisc[i]<<" Y "<< s1_doubleBDisc[i]<<endl;}
    // if (i==70){cout<<"for doubleB=0.3 X,Y: "<<b1_doubleBDisc[i]<<" Y "<< s1_doubleBDisc[i]<<endl;} // X 0.0841232 Y 0.89172 
    // if (b1_deepdoubleBDisc[i] > 0.082 && b1_deepdoubleBDisc[i] < 0.086) {cout<<i<<endl;}// i =57 (0.084634, 0.961783) 
    // if (i==57) {cout<<"deepdouble X "<<b1_deepdoubleBDisc[i]<<" Y "<< s1_deepdoubleBDisc[i]<<endl;}// i=57 means j=43
    //cout<<"Ratio TY  Sig/Total tau21 "<<i<<" is "<<s1[i]<<endl;
    //   cout<<"Ratio DY SigDisc/Total deepDisc "<<j<<" is "<<s1_ZDisc[j]<<endl;
    //cout<<"Ratio TX Bkg/Total tau21 "<<i<<" is "<<b1[i]<<endl;
    // cout<<"Ratio DX BkgDisc/Total deepDisc8"<<j<<" is "<<b1_ZDisc[j]<<endl;
    cout<<"***************"<<endl;
  }

  //
  // Creating TGraph for ROC curves 
  TGraph *roc_graph = new TGraph(101, b1, s1); // tau21
  TGraph *roc_graph1 = new TGraph(101, b1_ZDisc, s1_ZDisc); // deepAK8Zdisc.
  TGraph *roc_graph2 = new TGraph(101, b1_doubleBDisc, s1_doubleBDisc); // double B disc.
  TGraph *roc_graph3 = new TGraph(101, b1_deepdoubleBDisc, s1_deepdoubleBDisc); // deep double B disc.
  //  TGraph *roc_graph4 = new TGraph(101, b1_ZHDisc, s1_ZHDisc); // deepAK8ZH DIS. MD
  
  roc_graph->SetLineColor(2);
  roc_graph->SetLineWidth(2);
  // roc_graph->SetMarkerColor(2);
  // roc_graph->SetMarkerSize(1.5);
  // roc_graph->SetMarkerStyle(22);
  
  roc_graph1->SetLineColor(4);
  roc_graph1->SetLineWidth(2);
  // roc_graph1->SetMarkerColor(4);
  // roc_graph1->SetMarkerSize(1.5);
  // roc_graph1->SetMarkerStyle(21);

  roc_graph2->SetLineColor(6);
  roc_graph2->SetLineWidth(2);
  if (i==70){
   roc_graph2->SetMarkerColor(6);
   roc_graph2->SetMarkerSize(1.5);
   roc_graph2->SetMarkerStyle(23);
  }

  roc_graph3->SetLineColor(1);
  roc_graph3->SetLineWidth(2);
  if (i==57){
  roc_graph3->SetMarkerColor(8);
  roc_graph3->SetMarkerSize(1.5);
  roc_graph3->SetMarkerStyle(24);
  }

  // roc_graph4->SetLineColor(3);
  // roc_graph4->SetLineWidth(2);

  //
  // Drawing on canvas
  TCanvas *c2 = new TCanvas("c2","ROC plot", 200,10,700,500);
  //c2->SetBorderMode(0);
  c2->SetGridx();
  c2->SetGridy();
  c2->Modified();
  c2->Update();
  name = filename+".pdf";
  //c2->SaveAs(name);
  TMultiGraph *mg = new TMultiGraph();
 
  mg->Add(roc_graph, "LP");
  //mg->Add(roc_graph1,"LP");
  mg->Add(roc_graph2,"LP");
  mg->Add(roc_graph3,"LP");
  //  mg->Add(roc_graph4,"LP");
  mg->Draw("ALP");
  mg->SetTitle(" ROC curve J1  Sig:TChiWZ Bkg: Total Bkg");
  mg->GetXaxis()->SetTitle("Bkg selection");
  mg->GetYaxis()->SetTitle("Signal selection");
  //  gPad->SetLogx();
  auto legend = new TLegend(0.45,0.4,0.9,0.7);
  //  legend->AddEntry(roc_graph1,"AK8J1 Deep Z Disc.","l");
  legend->AddEntry(roc_graph,"AK8J1 Tau21","l");
  legend->AddEntry(roc_graph2,"AK8J1 W Disc.","l");
  legend->AddEntry(roc_graph3,"AK8J1 W Disc. MD","l");
  //  legend->AddEntry(roc_graph4,"AK8J1 Deep ZH Disc. MD","l");
  legend->Draw();
 
  // TF1 f1("f",[&](double *s1, double *){ return roc_graph->Eval(s1[0]); },0,1,0);
  // double integral = f1.Integral(0,1);
  // cout<<"integral"<<integral<<endl; // gives area under the curve

  // double AOC_roc_graph = roc_graph->Integral(0, -1);
  // cout<<"AOC_roc_graph"<<AOC_roc_graph<<endl;
  //roc_graph->Draw("A*L");
  //roc_graph1->Draw("sames");
  //file_inputdistributions1->Close();
  //file_inputdistributions2->Close();

  // Save drawn plots to a file
  c2->SaveAs(name);
  
}


