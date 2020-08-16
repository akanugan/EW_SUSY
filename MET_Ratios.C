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

void tmp(){
  TFile * file_inputdistributions;
  TFile * file_inputdistributions2;
  TFile * file_inputdistributions3;
  TFile * file_output;
  TString filename = "Ratio_highSRWandZ";     //file name to be saved
  TString name;
  

  file_inputdistributions = new TFile("wandzbkg.root", "READ");
  //file_inputdistributions = new TFile("TTJets_MC2018.root", "READ");
  //BKG_WX_conditions1a/ Histos_deep2b_mtbcut_massSB/
  //ST__MC2018.root //     TTJets_MC2018.root
  // file_inputdistributions1 = new TFile("BKG_baseline+Msd45to65_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_NOtau21_binned1/ZJetsToNuNu_HT_MC2018.root", "READ");
  //file_inputdistributions3 = new TFile("BKG_baseline+Msd100to120_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions1 = new TFile("BKG_baseline+T21lessthan0.35_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_baseline+0.35T21to0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions3 = new TFile("BKG_baseline+T21greaterthan0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");

  TH1F * h_Bkg1 = (TH1F*) file_inputdistributions->Get("whMETvBin"); //wzMETvBin
  TH1F * h_Bkg2 = (TH1F*) file_inputdistributions->Get("whMETvBin_RegC");  // wzMETvBin_RegC WHAK8J1MassNo2bTag      //WHAK8J1Massnobtag
  TH1F * h_Bkg3 = (TH1F*) file_inputdistributions->Get("whMETvBin_RegB");  //WHAK8J1MassNo2bTag      //WHAK8J1Massnobtag
  TH1F * h_Bkg4 = (TH1F*) file_inputdistributions->Get("whMETvBin_RegD");  //WHAK8J1MassNo2bTag      //WHAK8J1Massnobtag
  // TH1F * h_Bkg3 = (TH1F*) file_inputdistributions->Get("HWAK8J1MassSB");
  // TH1F * h_Bkg4 = (TH1F*) file_inputdistributions->Get("HWAK8J1MassNo2bTag");

  // h_Bkg1->Add(h_Bkg3); // higgs mass tagged
  // h_Bkg2->Add(h_Bkg4); // higgs mass no tag

  //  const Int_t NBINS = 4;
  //Double_t edges[NBINS + 1] = {50.0, 85.0, 135.0, 160.0, 190.0, 250.0};
  //  Double_t edges[NBINS + 1] = {30.0, 65.0, 135., 200.0, 250.0};
  TH1F *h_varbin = (TH1F*)h_Bkg1->Rebin(1);
  TH1F *h_varbin2 = (TH1F*)h_Bkg2->Rebin(1);
  h_varbin->Divide(h_varbin2);
  h_Bkg3->Divide(h_Bkg4);
  h_varbin->GetXaxis()->SetRange(0,300);
  h_varbin->GetYaxis()->SetRangeUser(0.,3.0);
  h_varbin->SetTitle("");
  h_varbin->GetXaxis()->SetTitle("MET");
  h_varbin->GetYaxis()->SetTitle("Ratios  {pass/fail}");
  h_varbin->SetStats(0);
  h_varbin->GetXaxis()->SetTitleSize(.7);
  h_varbin->GetYaxis()->SetTitleSize(.7);
  h_varbin->GetXaxis()->SetTitleSize(20);
  h_varbin->GetYaxis()->SetTitleSize(20);
  h_varbin->GetXaxis()->SetTitleFont(43); 
  h_varbin->GetYaxis()->SetTitleFont(43); 
  h_varbin->SetLineColor(kRed);
  h_Bkg3->SetLineColor(kBlue);


 
  auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  //c1->SetFillColor(42);
  c1->SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);

  h_varbin->Draw("E");
  h_Bkg3->Draw("E SAME");

  TLegend *legend1=new TLegend(0.6, 0.7, 0.9, 0.88);
  legend1->AddEntry(h_varbin,"W+Z Bkg","l");
  // legend1->AddEntry(h_Bkg3,"SB TTbar+Rare Bkg","l");
  legend1->Draw();

  c1->Print(filename+"wh.png");
  c1->Update();
  //c1->SaveAs(filename+".pdf");
}

  
