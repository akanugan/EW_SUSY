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

void ABCDratios(){
  TFile * file_inputdistributions;
  TFile * file_inputdistributions2;
  TFile * file_inputdistributions3;
  TFile * file_output;
  TString filename = "TTJets";     //file name to be saved
  TString name;
  

  file_inputdistributions = new TFile("Histos_WH_info5/TTJets_MC2018.root", "READ");
  //file_inputdistributions = new TFile("TTJets_MC2018.root", "READ");
  //BKG_WX_conditions1a/ Histos_deep2b_mtbcut_massSB/
  //ST__MC2018.root //     TTJets_MC2018.root
  // file_inputdistributions1 = new TFile("BKG_baseline+Msd45to65_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_NOtau21_binned1/ZJetsToNuNu_HT_MC2018.root", "READ");
  //file_inputdistributions3 = new TFile("BKG_baseline+Msd100to120_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions1 = new TFile("BKG_baseline+T21lessthan0.35_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_baseline+0.35T21to0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions3 = new TFile("BKG_baseline+T21greaterthan0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");

  TH1F * h_Bkg1 = (TH1F*) file_inputdistributions->Get("WHAK8J2MassSB");
  TH1F * h_Bkg2 = (TH1F*) file_inputdistributions->Get("WHAK8J2MassNo2bTag");  //WHAK8J1MassNo2bTag      //WHAK8J1Massnobtag
  TH1F * h_Bkg3 = (TH1F*) file_inputdistributions->Get("HWAK8J1MassSB");
  TH1F * h_Bkg4 = (TH1F*) file_inputdistributions->Get("HWAK8J1MassNo2bTag");

  h_Bkg1->Add(h_Bkg3); // higgs mass tagged
  h_Bkg2->Add(h_Bkg4); // higgs mass no tag

  const Int_t NBINS = 4;
  //Double_t edges[NBINS + 1] = {50.0, 85.0, 135.0, 160.0, 190.0, 250.0};
  Double_t edges[NBINS + 1] = {50.0, 85.0, 135.0, 200., 250.0};
  TH1F *h_varbin = (TH1F*)h_Bkg1->Rebin(NBINS,"hvarbin",edges);
  TH1F *h_varbin2 = (TH1F*)h_Bkg2->Rebin(NBINS,"hvarbin2",edges);
  h_varbin->Divide(h_varbin2);
  h_varbin->GetXaxis()->SetRange(0,300);
  h_varbin->GetYaxis()->SetRangeUser(0.,0.7);
  h_varbin->SetTitle("");
  h_varbin->GetXaxis()->SetTitle("Soft drop mass (GeV)");
  h_varbin->GetYaxis()->SetTitle("Ratio {pass/fail}");
  h_varbin->SetStats(0);
  h_varbin->GetXaxis()->SetTitleSize(.7);
  h_varbin->GetYaxis()->SetTitleSize(.7);
  h_varbin->GetXaxis()->SetTitleSize(20);
  h_varbin->GetYaxis()->SetTitleSize(20);
  h_varbin->GetXaxis()->SetTitleFont(43); 
  h_varbin->GetYaxis()->SetTitleFont(43); 

  int rebin = 2;
  h_Bkg1->Rebin(rebin);
  h_Bkg2->Rebin(rebin);
  
  h_Bkg1->Divide(h_Bkg2); // bin by bin
  h_Bkg1->GetXaxis()->SetRange(0,300);
  h_Bkg1->GetYaxis()->SetRangeUser(0.,0.7);
  h_Bkg1->SetTitle(""); // gStyle->SetOptTitle(0);
  h_Bkg1->GetXaxis()->SetTitle("Soft drop mass (GeV)");
  h_Bkg1->GetYaxis()->SetTitle("R_{pass/fail}");
  h_Bkg1->SetStats(0);
  h_Bkg1->GetXaxis()->SetTitleSize(20);
  h_Bkg1->GetYaxis()->SetTitleSize(20);
  h_Bkg1->GetXaxis()->SetTitleFont(43); 
  h_Bkg1->GetYaxis()->SetTitleFont(43); 

  
  

  auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  //c1->SetFillColor(42);
  c1->SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  
  h_Bkg1->Draw("E");
  c1->Print(filename+"i_deepbb.png");
  h_varbin->Draw("E");
  c1->Print(filename+"_deepbb.png");
  c1->Update();
  //c1->SaveAs(filename+".pdf");
}

  
