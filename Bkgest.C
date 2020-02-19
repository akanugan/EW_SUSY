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

void Bkgest(){
  TFile * file_inputdistributions;
  TFile * file_inputdistributions2;
  TFile * file_inputdistributions3;
  TFile * file_output;
  TString filename = "TotalBkg";     //file name to be saved
  TString name;
  

  file_inputdistributions = new TFile("BKG_WX_conditions1/Total_MC2018.root", "READ");

  // file_inputdistributions1 = new TFile("BKG_baseline+Msd45to65_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_NOtau21_binned1/ZJetsToNuNu_HT_MC2018.root", "READ");
  //file_inputdistributions3 = new TFile("BKG_baseline+Msd100to120_NoT21/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions1 = new TFile("BKG_baseline+T21lessthan0.35_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions2 = new TFile("BKG_baseline+0.35T21to0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");
  // file_inputdistributions3 = new TFile("BKG_baseline+T21greaterthan0.5_Nomasscut/ZJetsToNuNu_HT_MC2018.root", "READ");

  TH1F * h_Bkg1 = (TH1F*) file_inputdistributions->Get("WHAK8J2Mass");
  TH1F * h_Bkg2 = (TH1F*) file_inputdistributions->Get("WHAK8J2Massnobtag");
  TH1F * h_Bkg3 = (TH1F*) file_inputdistributions->Get("HWAK8J1Mass");
  TH1F * h_Bkg4 = (TH1F*) file_inputdistributions->Get("HWAK8J1Massnobtag");

  h_Bkg1->Add(h_Bkg3);
  h_Bkg2->Add(h_Bkg4);

  int n_bin_1, n_bin_2, n_bin_3, n_bin_4;
  n_bin_1 = h_Bkg1->GetNbinsX();
  n_bin_2 = h_Bkg2->GetNbinsX();
  n_bin_3 = h_Bkg3->GetNbinsX();
  n_bin_4 = h_Bkg4->GetNbinsX();
  
  //h_sig->Scale(1/h_sig->Integral(0, n_bin_s));
  //h_bkg->Scale(1/h_bkg->Integral(0, n_bin_b));
  
  cout<<n_bin_1<<"\t"<<n_bin_2<<"\t"<<n_bin_3<<"\t"<<n_bin_4<<endl;
  Double_t Bkg1_A = 0,Bkg1_B =0 ,Bkg1_C =0 , Bkg2_A =0 ,Bkg2_B =0 ,Bkg2_C =0 , Bkg3_A =0 ,Bkg3_B =0 ,Bkg3_C =0 , Bkg4_A =0 ,Bkg4_B =0 ,Bkg4_C =0, R1,R2,R3,P1,P2,P3;
  Double_t Bkg1_A_T = 0,Bkg1_B_T =0 ,Bkg1_C_T =0 , Bkg2_A_T =0 ,Bkg2_B_T =0 ,Bkg2_C_T =0 , Bkg3_A_T =0 ,Bkg3_B_T =0 ,Bkg3_C_T =0 , Bkg4_A_T =0 ,Bkg4_B_T =0 ,Bkg4_C_T =0;
  Double_t Bkg1_A_err2_T =0 , Bkg1A_error =0, Bkg2_A_err2_T =0 , Bkg2A_error =0,  Bkg3_A_err2_T =0 , Bkg3A_error =0,  Bkg4_A_err2_T =0 , Bkg4A_error =0;
  Double_t Bkg1_B_err2_T =0 , Bkg1B_error =0, Bkg2_B_err2_T =0 , Bkg2B_error =0,  Bkg3_B_err2_T =0 , Bkg3B_error =0,  Bkg4_B_err2_T =0 , Bkg4B_error =0;
  Double_t Bkg1_C_err2_T =0 , Bkg1C_error =0, Bkg2_C_err2_T =0 , Bkg2C_error =0,  Bkg3_C_err2_T =0 , Bkg3C_error =0,  Bkg4_C_err2_T =0 , Bkg4C_error =0;

  for(double bin = 0; bin <= 17; bin ++){ //Region A
    
    Bkg1_A=  h_Bkg1->GetBinContent(bin); 
    Bkg2_A=  h_Bkg2->GetBinContent(bin); 
    // Bkg3_A=  h_Bkg3->GetBinContent(bin); 
    // Bkg4_A=  h_Bkg4->GetBinContent(bin); 

    Bkg1_A_err2_T +=  TMath::Power(h_Bkg1->GetBinError(bin),2); 
    Bkg2_A_err2_T +=  TMath::Power(h_Bkg2->GetBinError(bin),2); 
    // Bkg3_A_err2_T +=  TMath::Power(h_Bkg3->GetBinError(bin),2); 
    // Bkg4_A_err2_T +=  TMath::Power(h_Bkg4->GetBinError(bin),2); 
    
    Bkg1_A_T += Bkg1_A;
    Bkg2_A_T += Bkg2_A;
    // Bkg3_A_T += Bkg3_A;
    // Bkg4_A_T += Bkg4_A;
    
  }
  Bkg1A_error= TMath::Sqrt(Bkg1_A_err2_T);
  Bkg2A_error= TMath::Sqrt(Bkg2_A_err2_T);
  // Bkg3A_error= TMath::Sqrt(Bkg3_A_err2_T);
  // Bkg4A_error= TMath::Sqrt(Bkg4_A_err2_T);
      
  cout << "Bkg 1 error: "<< Bkg1A_error <<endl;
    
  for(double bin =18; bin <= 27; bin ++){ //Region B
    
    Bkg1_B=  h_Bkg1->GetBinContent(bin); 
    Bkg2_B=  h_Bkg2->GetBinContent(bin); 
    // Bkg3_B=  h_Bkg3->GetBinContent(bin);  
    // Bkg4_B=  h_Bkg4->GetBinContent(bin); 
  
    Bkg1_B_err2_T +=  TMath::Power(h_Bkg1->GetBinError(bin),2); 
    Bkg2_B_err2_T +=  TMath::Power(h_Bkg2->GetBinError(bin),2); 
    // Bkg3_B_err2_T +=  TMath::Power(h_Bkg3->GetBinError(bin),2); 
    // Bkg4_B_err2_T +=  TMath::Power(h_Bkg4->GetBinError(bin),2); 
   
    Bkg1_B_T += Bkg1_B;
    Bkg2_B_T += Bkg2_B;
    // Bkg3_B_T += Bkg3_B;
    // Bkg4_B_T += Bkg4_B;
   
  }
  
  Bkg1B_error= TMath::Sqrt(Bkg1_B_err2_T);
  Bkg2B_error= TMath::Sqrt(Bkg2_B_err2_T);
  // Bkg3B_error= TMath::Sqrt(Bkg3_B_err2_T);
  // Bkg4B_error= TMath::Sqrt(Bkg4_B_err2_T);
  
  for(double bin = 28; bin <= 60; bin ++){ //Region C
    
    Bkg1_C=  h_Bkg1->GetBinContent(bin);   
    Bkg2_C=  h_Bkg2->GetBinContent(bin); 
    // Bkg3_C=  h_Bkg3->GetBinContent(bin); 
    // Bkg4_C=  h_Bkg4->GetBinContent(bin); 
    
    Bkg1_C_err2_T +=  TMath::Power(h_Bkg1->GetBinError(bin),2); 
    Bkg2_C_err2_T +=  TMath::Power(h_Bkg2->GetBinError(bin),2); 
    // Bkg3_C_err2_T +=  TMath::Power(h_Bkg3->GetBinError(bin),2); 
    // Bkg4_C_err2_T +=  TMath::Power(h_Bkg4->GetBinError(bin),2); 
  
    Bkg1_C_T += Bkg1_C;
    Bkg2_C_T += Bkg2_C;
    // Bkg3_C_T += Bkg3_C;
    // Bkg4_C_T += Bkg4_C;
   
  }
  
  Bkg1C_error= TMath::Sqrt(Bkg1_C_err2_T);
  Bkg2C_error= TMath::Sqrt(Bkg2_C_err2_T);
  // Bkg3C_error= TMath::Sqrt(Bkg3_C_err2_T);
  // Bkg4C_error= TMath::Sqrt(Bkg4_C_err2_T);
  

  Double_t R1_error =0, R2_error =0, R3_error =0, P1_error =0, P2_error =0, P3_error =0;
  R1 = Bkg1_A_T/Bkg2_A_T;
  R1_error = abs(R1)*TMath::Sqrt(TMath::Power(Bkg2A_error/Bkg2_A_T,2) + TMath::Power(Bkg1A_error/Bkg1_A_T,2));    // if C=A/B;  C_error=C*Sqrt((deltaA/A)^2+(deltaB/B)^2)
  // cout<<"r1 error: "<<R1_error<<endl;
  R2 = Bkg1_B_T/Bkg2_B_T;
  R2_error = abs(R2)*TMath::Sqrt(TMath::Power(Bkg2B_error/Bkg2_B_T,2) + TMath::Power(Bkg1B_error/Bkg1_B_T,2));   
  
  R3 = Bkg1_C_T/Bkg2_C_T;
  R3_error = abs(R3)*TMath::Sqrt(TMath::Power(Bkg2C_error/Bkg2_C_T,2) + TMath::Power(Bkg1C_error/Bkg1_C_T,2));   
 
  // P1 = Bkg3_A_T/Bkg4_A_T;
  // P1_error = abs(P1)*TMath::Sqrt(TMath::Power(Bkg4A_error/Bkg4_A_T,2) + TMath::Power(Bkg3A_error/Bkg3_A_T,2));   
  // P2 = Bkg3_B_T/Bkg4_B_T;
  // P2_error = abs(P2)*TMath::Sqrt(TMath::Power(Bkg4B_error/Bkg4_B_T,2) + TMath::Power(Bkg3B_error/Bkg3_B_T,2));   
  // P3 = Bkg3_C_T/Bkg4_C_T;
  // P3_error = abs(P3)*TMath::Sqrt(TMath::Power(Bkg4C_error/Bkg4_C_T,2) + TMath::Power(Bkg3C_error/Bkg3_C_T,2));   

  cout<< " B1 = "<< Bkg1_A_T<<"\u00b1"<<Bkg1A_error << " | A1 = " << Bkg1_B_T<<"\u00b1"<<Bkg1B_error << "| B1 = " << Bkg1_C_T<<"\u00b1"<<Bkg1C_error <<"\n"<<endl;
  cout<< " D = "<< Bkg2_A_T<<"\u00b1"<<Bkg2A_error << " | C = " << Bkg2_B_T<<"\u00b1"<<Bkg2B_error << "| D = " << Bkg2_C_T<<"\u00b1"<<Bkg2C_error<<"\n"<<endl;
  cout<< " B1/D = "<< R1 <<"\u00b1"<< R1_error << " | A1/C = " << R2 <<"\u00b1"<< R2_error  << " | B1/D = " << R3 <<"\u00b1"<<R3_error<<"\n"<<endl;

  // cout<< " B1 = "<< Bkg3_A_T<<"\u00b1"<<Bkg3A_error << " | A1 = " << Bkg3_B_T<<"\u00b1"<<Bkg3B_error << "| B1 = " << Bkg3_C_T<<"\u00b1"<<Bkg3C_error<<"\n"<<endl;
  // cout<< " D = "<< Bkg4_A_T<<"\u00b1"<<Bkg4A_error << " | C = " << Bkg4_B_T<<"\u00b1"<<Bkg4B_error << "| D = " << Bkg4_C_T<<"\u00b1"<<Bkg4C_error<<"\n"<<endl;
  // cout<< " B1/D = "<< P1 <<"\u00b1"<< P1_error<< " | A1/C = " << P2 <<"\u00b1"<< P2_error << " | B1/D = " << P3<<"\u00b1"<<P3_error<<"\n"<<endl;


  // cout<<"Bkg1 Tau21<0.35 = "<< Bkg1_A << " || Bkg1 Tau21>0.35 = " << Bkg1_B << " || Bkg1 (Tau21<0.35)/(Tau21>0.35) =  "<< Bkg1_R <<endl;
  // cout<<"Bkg2 Tau21<0.35 = "<< Bkg2_A << " || Bkg2 Tau21>0.35 = " << Bkg2_B << " || Bkg2 (Tau21<0.35)/(Tau21>0.35) = "<< Bkg2_R <<endl;
  // cout<<"Bkg3 Tau21<0.35 = "<< Bkg3_A << " || Bkg3 Tau21>0.35 = " << Bkg3_B << " || Bkg3 (Tau21<0.35)/(Tau21>0.35) = "<< Bkg3_R <<endl;

  // auto c1 = new TCanvas("c1","Profile histogram",200,10,700,500);
  // double massbins[4] ={50.,85.,135.,250.};
  
  // TGraphErrors* hprof  = new TGraphErrors("hprof","Profile of ratios",3,massbins);
  // //void SetHistFillColor(Color_t color = 1);
  // //void SetHistFillStyle(Style_t styl = 0);
  // hprof->SetBinEntries(1,1);
  // hprof->SetBinContent(1,R1);
  // hprof->SetBinError(  1,R1_error);
  // cout<<"r1 is: "<<R1<<endl;
  // cout<<"r1 error is: "<<R1_error<<endl;
  // hprof->SetBinEntries(2,1);
  // hprof->SetBinContent(2,R2);
  // hprof->SetBinError(  2,R2_error);
  // hprof->SetBinEntries(3,1);
  // hprof->SetBinContent(3,R3);
  // hprof->SetBinError(  3,R3_error);
  // //hprof->SetBinContent(0,100000);
  // std::cout<<"bin 1 content is "<<hprof->GetBinContent(1)<<std::endl;
  // std::cout<<"bin 0 error is "<<hprof->GetBinError(1)<<std::endl;

  // hprof->Draw("L");

  
  auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  //c1->SetFillColor(42);
  c1->SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  const Int_t n = 3;
  Double_t xbins[n+1]  = {50.,85.,135.,250.};
  Double_t x[n]  = {(xbins[0]+xbins[1])/2,(xbins[1]+xbins[2])/2,(xbins[2]+xbins[3])/2};
  //Double_t x[n]  = {67.5,110.,192.5};
  
  Double_t y[n]  = {R1,R2,R3};
  //  Double_t ex[n] = {17.5,25.,57.5};
  Double_t ex[n] = {(xbins[1]-xbins[0])/2,(xbins[2]-xbins[1])/2,(xbins[3]-xbins[2])/2};
  Double_t ey[n] = {R1_error,R2_error,R3_error};
  auto gr = new TGraphErrors(n,x,y,ex,ey);

  gr->SetTitle(filename);
  gr->GetYaxis()->SetTitle("R_{pass/fail}"); 
  gr->GetXaxis()->SetTitle("Softdrop Jet Mass (GeV)"); 
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");
  
  TLatex l;
  l.SetTextSize(0.050);
  //l.SetTextAngle(30.);
  // for (int i=0; i<n; i++) {
  //   l.DrawLatex(x[i],y[i],Form("(%g,%g)",x[i],y[i]));
  //  }
  l.DrawLatex(x[0]+0.01,y[0]+0.01,"B1/D");
  l.DrawLatex(x[1]+0.01,y[1]+0.01,"A1/C");
  l.DrawLatex(x[2]+0.01,y[2]+0.01,"B1/D");
  c1->SaveAs(filename+".pdf");
}

  
