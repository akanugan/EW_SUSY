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

void roc_KH(){
  TFile * file_inputdistributions1;
  TFile * file_inputdistributions2;
  TFile * file_output;
  TString filename = "ROC_plotWisJ1_ZjetsBkg";     //file name to be saved
  TString name;
  
  file_inputdistributions1 = new TFile("SIG_WZ_notau21_jet1isW/TChiWZ_600_100_MC2018.root", "READ");
  file_inputdistributions2 = new TFile("BKG_NOtau21_binned1/ZJetsToNuNu_HT_MC2018.root", "READ");

  // based on J1 
  /*
  TH1F * h_sigDisc = (TH1F*) file_inputdistributions1->Get("AK8J1wDis");
  TH1F * h_bkgDisc = (TH1F*) file_inputdistributions2->Get("AK8J1wDis");
  TH1F * h_sig = (TH1F*) file_inputdistributions1->Get("AK8Tau21");
  TH1F * h_bkg = (TH1F*) file_inputdistributions2->Get("AK8Tau21");
  */

  // based on J2
  TH1F * h_sigDisc = (TH1F*) file_inputdistributions1->Get("AK8J1wDis");
  TH1F * h_bkgDisc = (TH1F*) file_inputdistributions2->Get("AK8J1wDis");
  TH1F * h_sig = (TH1F*) file_inputdistributions1->Get("AK8Tau21");
  TH1F * h_bkg = (TH1F*) file_inputdistributions2->Get("AK8Tau21");
  

  int n_bin_s, n_bin_b;
  n_bin_s = h_sig->GetNbinsX();
  n_bin_b = h_bkg->GetNbinsX();
  
  //h_sig->Scale(1/h_sig->Integral(0, n_bin_s));
  //h_bkg->Scale(1/h_bkg->Integral(0, n_bin_b));
  
  cout<<n_bin_s<<"\t"<<n_bin_b<<endl;
  cout<<h_sigDisc->GetNbinsX()<<"\t"<<h_bkgDisc->GetNbinsX()<<endl;
  cout<<"********"<<endl;

  Double_t s[101],s_Disc[101];
  Double_t b[101],b_Disc[101];
  Double_t s_total,s_totalDisc=0;
  Double_t b_total,b_totalDisc=0;
  Double_t s1[101],s1_Disc[101] = {0};
  Double_t b1[101],b1_Disc[101] = {0};


  TCanvas *c1 = new TCanvas("c1","Individual plot", 200,10,700,500);
  c1->cd();
  h_sig->Draw();
  h_bkg->Draw("sames");
  
  s_total=  h_sig->Integral(0, 100);
  b_total=  h_bkg->Integral(0, 100);
  s_totalDisc=  h_sigDisc->Integral(0, 100);
  b_totalDisc=  h_bkgDisc->Integral(0, 100);
  cout<<"sig total :"<<s_total<<endl;
  cout<<"bkg total :"<<b_total<<endl;
  int i;
  int j;
  for( i= 0 , j = 100; i <= 100; i++,j--){
    s[i] = h_sig->Integral(0, i);
    b[i] = h_bkg->Integral(0, i);
    s_Disc[j] = h_sigDisc->Integral(j, 100);
    b_Disc[j] = h_bkgDisc->Integral(j, 100);
    cout<<"sig bin cont in "<<i<<" is "<<s[i]<<endl;
    cout<<"bkg bin cont in "<<i<<" is "<<b[i]<<endl;
    cout<<"***********************************"<<endl;
    cout<<"sigdeepDisc bin cont in "<<j<<" is "<<s_Disc[j]<<endl;
    cout<<"bkgdeepDisc bin cont in "<<j<<" is "<<b_Disc[j]<<endl;
    cout<<"**********Ratios********"<<endl;  


    s1[i] = s[i]/s_total;
    b1[i] = b[i]/b_total;
    s1_Disc[j] = s_Disc[j]/s_totalDisc;
    b1_Disc[j] = b_Disc[j]/b_totalDisc;
    cout<<"Ratio TY  Sig/Total tau21 "<<i<<" is "<<s1[i]<<endl;
    cout<<"Ratio DY SigDisc/Total deepDisc "<<j<<" is "<<s1_Disc[j]<<endl;
    cout<<"Ratio TX Bkg/Total tau21 "<<i<<" is "<<b1[i]<<endl;
    cout<<"Ratio DX BkgDisc/Total deepDisc8"<<j<<" is "<<b1_Disc[j]<<endl;
    cout<<"***************"<<endl;
  }

  //
  // Creating TGraph for ROC curves 
  TGraph *roc_graph = new TGraph(101, b1, s1); // tau21
  TGraph *roc_graph1 = new TGraph(101, b1_Disc, s1_Disc); // deepAK8disc.
  roc_graph->SetLineColor(2);
  roc_graph->SetLineWidth(2);
  roc_graph->SetMarkerColor(2);
  roc_graph->SetMarkerSize(1.5);
  roc_graph->SetMarkerStyle(22);
  
  roc_graph1->SetLineColor(4);
  roc_graph1->SetLineWidth(2);
  roc_graph1->SetMarkerColor(4);
  roc_graph1->SetMarkerSize(1.5);
  roc_graph1->SetMarkerStyle(21);
  
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
  mg->Add(roc_graph1,"LP");
  mg->Draw("ALP");
  mg->SetTitle(" ROC curve J1  Sig:TChiWZ Bkg: Z+jets Bkg");
  mg->GetXaxis()->SetTitle("Bkg selection");
  mg->GetYaxis()->SetTitle("Signal selection");

  auto legend = new TLegend(0.1,0.8,0.48,0.9);
  legend->AddEntry(roc_graph1,"DeepAk8J1 W Disc.","l");
  legend->AddEntry(roc_graph,"AK8J1 Tau21","l");
  legend->Draw();
 
  TF1 f1("f",[&](double *s1, double *){ return roc_graph->Eval(s1[0]); },0,1,0);
  double integral = f1.Integral(0,1);
  cout<<"integral"<<integral<<endl; // gives area under the curve

  double AOC_roc_graph = roc_graph->Integral(0, -1);
  cout<<"AOC_roc_graph"<<AOC_roc_graph<<endl;
  //roc_graph->Draw("A*L");
  //roc_graph1->Draw("sames");
  //file_inputdistributions1->Close();
  //file_inputdistributions2->Close();

  // Save drawn plots to a file
  c2->SaveAs(name);
  
}


