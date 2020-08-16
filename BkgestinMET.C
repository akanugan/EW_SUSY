#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"
#include"TMath.h"

//gROOT->ProcessLine(".L Bkgest_conditions1METbinned.C");

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=7,nBG=6;    //Specify no. of files
TFile *f[nfiles];
int col[14]={kPink+1,kMagenta,kGreen,kBlue,kRed,kYellow,kOrange,kCyan,kBlue+2,kPink+1,kOrange+1,kGreen+3,kBlue+7};  //Specify Colors b's

TCanvas *c_cA=new TCanvas("kinVar","",1500,900);
//TCanvas *test=new TCanvas("test","",1500,900);

void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
TLegend *legend1=new TLegend(0.4501, 0.65,  0.88, 0.88);
//TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);
//TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);

//TLegend *legend2=new TLegend(0.7, 0.9,  0.90, 0.65);
//TLegend *legend2=new TLegend(0.6, 0.90,  0.98, 0.45);
void setLastBinAsOverFlow(TH1D*);
void BkgestinMET(){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);

  TString varName =  "wzMET";
  TString varName2 = "wzMETvBin_RegB";
  TString varName3 = "wzMETvBin_RegC";
  TString varName4 = "wzMETvBin_RegD";

  // TString varName5 = "WHMET_RegC";
  // TString varName6 = "HWMET_RegC";
  // TString varName7 = "WHMET_RegD";
  // TString varName8 = "HWMET_RegD";
  TString xLabel = " MET (GeV)";
  int rebin=1;

  f[0] = new TFile("ST__MC2018.root");
  f[1] = new TFile("Rare_MC2018.root");
  f[2] = new TFile("TTJets_MC2018.root");
  f[3] = new TFile("QCD_HT_MC2018.root");
  f[4] = new TFile("WJetsToLNu_HT_MC2018.root");
  f[5] = new TFile("ZJetsToNuNu_HT_MC2018.root");
  f[6] = new TFile("Total_MC2018.root");
 

  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","");
  //TH1D *h_R;
  TH1D *h_MET_R[nfiles];
  for(int i=0;i<nfiles;i++){
    sprintf(name,"hist_file%i",i);
    h_MET_R[i]=new TH1D(name,name,21,0.5,21.5);
  }
  vector<double> Bcnt;
  double intLumi=137;
  TLatex tl1;

  static const int numMETbins = 6;
  Double_t METbins[numMETbins+1] = {250.,350.,500.,600.,700., 800.,2000.};
  
  // //rebin total
  TH1D *h_METTotal=(TH1D*)f[6]->FindObjectAny(varName);
  TH1D *h_METTotal2=(TH1D*)f[6]->FindObjectAny(varName2);
  TH1D *h_METTotal3=(TH1D*)f[6]->FindObjectAny(varName3);
  TH1D *h_METTotal4=(TH1D*)f[6]->FindObjectAny(varName4);
  
  // TH1D *h_METTotal5=(TH1D*)f[6]->FindObjectAny(varName5);
  // TH1D *h_METTotal6=(TH1D*)f[6]->FindObjectAny(varName6);
  // TH1D *h_METTotal7=(TH1D*)f[6]->FindObjectAny(varName7);
  // TH1D *h_METTotal8=(TH1D*)f[6]->FindObjectAny(varName8);


  // h_METTotal->Add(h_METTotal2); // region a
  // h_METTotal3->Add(h_METTotal4); //region b
  // h_METTotal5->Add(h_METTotal6); // region c
  // h_METTotal7->Add(h_METTotal8); //region d
  
  TH1D* h_A = (TH1D*)h_METTotal->Rebin(numMETbins,"h_A",METbins);  
  TH1D* h_B = (TH1D*)h_METTotal2->Rebin(numMETbins,"h_B",METbins);  
  TH1D* h_C = (TH1D*)h_METTotal3->Rebin(numMETbins,"h_C",METbins);  
  TH1D* h_D = (TH1D*)h_METTotal4->Rebin(numMETbins,"h_D",METbins);  

  TH1D* overlay = new TH1D("overlay","Predicted BC/D",numMETbins,METbins);
 
  for (int bin=1; bin < numMETbins+1; bin++){  //total
    //for (int bin=0; bin < numMETbins; bin++){
    Double_t A,B,C,D,Pred_A =0;
    Double_t B_err,C_err,D_err,Pred_A_err =0;
    B = h_B->GetBinContent(bin); // B 
    printf("%d\n",bin);
    cout<<"B is: "<<B<<endl;
    B_err = h_B->GetBinError(bin); // B error 
    cout<<"B err is: "<<B_err<<endl;
    C = h_C->GetBinContent(bin); // C
    cout<<"C is: "<<C<<endl; 
    C_err = h_C->GetBinError(bin); // C error 
    cout<<"C err is: "<<C_err<<endl;
    D = h_D->GetBinContent(bin); // D 
    cout<<"D is: "<<D<<endl;
    D_err = h_D->GetBinError(bin); // D error 
    cout<<"D err is: "<<D_err<<endl;
    Pred_A = abs((B*C)/D);
    cout<<"pred A is "<<Pred_A<<endl;
    Pred_A_err = Pred_A * TMath::Sqrt( TMath::Power(B_err/B,2) + TMath::Power(C_err/C,2) + TMath::Power(D_err/D,2));
    cout<<"pred A error  is "<<Pred_A_err<<endl;

    A = h_A->GetBinContent(bin); // A  
    cout<<"Observed A is "<< A <<endl;
    cout <<"Difference : " << abs(Pred_A - A) <<endl; 
      
    overlay->SetBinContent(bin,Pred_A);
    overlay->SetBinError(bin,Pred_A_err);
    TAxis *xaxis = overlay->GetXaxis();
    Double_t binCenter = xaxis->GetBinCenter(bin);
    cout<<"bin center is "<<binCenter<<endl;
    cout<<"*************"<<endl;
    
  }
  
  overlay->Draw("HIST");

  
  TH1D *h_total;
  //individual
  for(int i=0; i < nBG; i++){ //0 to 5
    
    //    printf("%d\n",i);
    
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    TH1D *h_MET2=(TH1D*)f[i]->FindObjectAny(varName2);
    TH1D *h_MET3=(TH1D*)f[i]->FindObjectAny(varName3);
    TH1D *h_MET4=(TH1D*)f[i]->FindObjectAny(varName4);

    // TH1D *h_MET5=(TH1D*)f[i]->FindObjectAny(varName5);
    // TH1D *h_MET6=(TH1D*)f[i]->FindObjectAny(varName6);
    // TH1D *h_MET7=(TH1D*)f[i]->FindObjectAny(varName7);
    // TH1D *h_MET8=(TH1D*)f[i]->FindObjectAny(varName8);

    // h_MET->Add(h_MET2); // region a
    
    // h_MET3->Add(h_MET4); //region b
    // h_MET5->Add(h_MET6); // region c
    // h_MET7->Add(h_MET8); //region d
        
    TH1D* h_MET_binned = (TH1D*)h_MET->Rebin(numMETbins,"h_MET_binned",METbins);
 
    TH1D* h_MET2_binned = (TH1D*)h_MET2->Rebin(numMETbins,"h_MET2_binned",METbins);
    TH1D* h_MET3_binned = (TH1D*)h_MET3->Rebin(numMETbins,"h_MET3_binned",METbins);
    TH1D* h_MET4_binned = (TH1D*)h_MET4->Rebin(numMETbins,"h_MET4_binned",METbins);
 
    cout<<"Bkg:  "<< f[i]->GetName() <<endl;
    cout <<"A in 1st bin is: " <<h_MET_binned->GetBinContent(1)<<" +/- "<<h_MET_binned->GetBinError(1)<<endl; 
    cout <<"B in 1st bin is: " <<h_MET2_binned->GetBinContent(1)<<" +/- "<<h_MET2_binned->GetBinError(1)<<endl; 
    cout <<"C in 1st bin is: " <<h_MET3_binned->GetBinContent(1)<<" +/- "<<h_MET3_binned->GetBinError(1)<<endl; 
    cout <<"D in 1st bin is: " <<h_MET4_binned->GetBinContent(1)<<" +/- "<<h_MET4_binned->GetBinError(1)<<endl; 
    cout<<" \n  "<<endl;
    cout <<"A in 2nd bin is: " <<h_MET_binned->GetBinContent(2)<<" +/- "<<h_MET_binned->GetBinError(2)<<endl; 
    cout <<"B in 2nd bin is: " <<h_MET2_binned->GetBinContent(2)<<" +/- "<<h_MET2_binned->GetBinError(2)<<endl; 
    cout <<"C in 2nd bin is: " <<h_MET3_binned->GetBinContent(2)<<" +/- "<<h_MET3_binned->GetBinError(2)<<endl; 
    cout <<"D in 2nd bin is: " <<h_MET4_binned->GetBinContent(2)<<" +/- "<<h_MET4_binned->GetBinError(2)<<endl; 
    cout<<" ***  \n  "<<endl;

 
    if(i<=(nBG-1)){
      hs_var->Add(h_MET_binned);
      if(i==0) h_total = (TH1D*)h_MET->Clone("totalHist");
      else h_total->Add(h_MET);
    }
    TH1D* h_total_rebinned = (TH1D*)h_total->Rebin(numMETbins,"h_total_rebinned",METbins);  
    decorate(h_MET_binned,i,f[i]->GetName());
    
    if(i==nBG-1) {
      c_cA->cd();
      hs_var->Draw("BAR HIST");
      //hs_var->Draw("colz");
      hs_var->Draw("HIST");
      // hs_var->SetMinimum(0.8);
      hs_var->SetMaximum(hs_var->GetMaximum()*10);
      h_total_rebinned->SetFillStyle(3017);
      h_total_rebinned->SetFillColor(kBlack);
      h_total_rebinned->Draw("e2 same");

      //   decorate(hs_var,i,f[i]->GetName()); 
      //hs_var->GetYaxis()->SetRangeUser(100.5,20000);
    }
    if(i>=nBG){ 
      //overlay BC/D here
      // c_cA->cd();
      // overlay->SetMarkerStyle(20);
      // overlay->SetMarkerColor(col[7]);
      // overlay->SetLineColor(col[8]);
      // overlay->SetLineWidth(3);
      // // overlay->Draw("HIST");
      // overlay->Draw("e2 hist same");      
      // overlay->Draw("HIST");
      // overlay->Draw("SAME");
      
    }
    drawlegend(h_MET_binned,i,f[i]->GetName());
    if(i==nBG-1){
      //hs_var->GetXaxis()->SetRangeUser(0,1500); 
      //hs_var->GetXaxis()->SetTitleOffset(.90);
      hs_var->GetXaxis()->SetTitle("MET");
      hs_var->GetYaxis()->SetTitle("Events");
      hs_var->SetTitle(0);
    }
  }
  
  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  c_cA->cd();
  //  overlay->SetMarkerStyle(4);
  overlay->SetMarkerColor(col[10]);
  overlay->SetLineColor(col[10]);
  overlay->SetLineWidth(3);
  // overlay->Draw("HIST");
  overlay->Draw("E SAME ");      
  // overlay->Draw("HIST");
  // overlay->Draw("SAME");
  drawlegend(overlay,6,f[6]->GetName());
  gPad->SetLogy();
  legend1->Draw();
  //  gPad->RedrawAxis();
  //  hs_var->GetXaxis()->SetTitle(xLabel);
  
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.68,0.91,name2);

  name3 = varName+".pdf";
  c_cA->SaveAs(name3);
  cout<<"*****************************************************"<<endl;
  cout<<"Int Lumi(inv.fb) for file1:"<<setprecision(4)<<intLumi<<endl;
    
}

void decorate(THStack *hs,int i,const char* fname){
  //  hs->SetMinimum(0.5);
  //hs->SetTitle(0);
  // hs->GetXaxis()->SetLabelSize(.05);
  // hs->GetYaxis()->SetLabelSize(.05);
  // hs->GetXaxis()->SetTitleSize(0.05);
  //hs->GetYaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  //  gPad->Update();
  gStyle->SetOptStat(0);
}
void decorate(TH1D* hist,int i,const char* fname){
  hist->SetLineColor(col[i]);
  if(i<nBG) {
    hist->SetFillColor(col[i]);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
  }
  else hist->SetLineWidth(2);
  hist->SetTitle(0);
  hist->GetXaxis()->SetLabelSize(.06);
  hist->GetYaxis()->SetLabelSize(.06);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.06);
  // drawlegend(hist,i,fname);
  //gPad->Update();
  setLastBinAsOverFlow(hist);
  gStyle->SetOptStat(0);
  //Hlist.Add(hist);
}

void drawlegend(TH1D *hist,int i,const char* fname){
  gStyle->SetLegendBorderSize(0);
 
  TString lName=fname;
  
  if(lName.Contains("ST_")){lName="Single t";}
  else if(lName.Contains("DYJetsToLL")){lName="DY(l^{+}l^{-})";}
  else if(lName.Contains("WJetsToLNu")){lName="W(l#nu)+ jets";}
  else if(lName.Contains("Rare")){lName="Rare";}
  else if(lName.Contains("TTJets")){lName="t #bar{t}";}
  else if(lName.Contains("QCD")){lName="QCD";}
  else if(lName.Contains("WGJetsToLNuG")){lName="W(l#nu)+ #gamma";}
  else if(lName.Contains("ZJetsToNuNu")){lName="Z(#nu#bar{#nu})+ jets";}
  else if(lName.Contains("TTGJets")){lName="t #bar{t}+ #gamma";}
  else if(lName.Contains("GJets")){lName="#gamma +jets";}
  else if(lName.Contains("Run2016")){lName="Data";}
  // else if(lName.Contains("TChiWZ_1000_1")){lName="TChiWZ_1000_1";}
  // else if(lName.Contains("TChiWZ_800_1")){lName="TChiWZ_800_1";}
  // else if(lName.Contains("TChiWZ_600_1")){lName="TChiWZ_600_1";}
  
  else if(lName.Contains("TChiWZ_300_1")){lName="TChiWZ_300_1";}
  
  else if(lName.Contains("TChiWZ_500_1")){lName="TChiWZ_500_1";}
  else if(lName.Contains("TChiWZ_600_100")){lName="TChiWZ(600,100)";}
  //  else if(lName.Contains("TChiWZ_600_100")){lName="TChiWZ_600_100";}
  else if(lName.Contains("TChiWZ_600_200")){lName="TChiWZ_600_200";}
  else if(lName.Contains("Total")){lName="Predicted B*C/D";}


  else if(lName.Contains("TChiWZ_600_400")){lName="TChiWZ_600_400";}
  else if(lName.Contains("TChiWZ_600_500")){lName="TChiWZ_600_500";}
  //  else if(lName.Contains("TChiWZ_800_1")){lName="TChiWZ_800_1";}
  else if(lName.Contains("TChiWZ_800_100")){lName="TChiWZ(800,100)";}
  else if(lName.Contains("TChiWZ_1000_100")){lName="TChiWZ(1000,100)";}
  else if(lName.Contains("TChiWZ_800_300")){lName="TChiWZ_800_300";}
  else if(lName.Contains("TChiWZ_800_500")){lName="TChiWZ_800_500";}
  else if(lName.Contains("TChiWZ_800_700")){lName="TChiWZ_800_700";}
  else if(lName.Contains("TChipmWW_600_100")){lName="TChipmWW_600_100";}
  else if(lName.Contains("TChipmWW_800_100")){lName="TChipmWW_800_100";}
  else if(lName.Contains("TChipmWW_1000_100")){lName="TChipmWW_1000_100";}
  // else if(lName.Contains("TChiWZ_900_100")){lName="TChiWZ_900_100";}
  //else if(lName.Contains("TChiWZ_1000_1")){lName="TChiWZ_1000_1";}

  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="T5bbbbZG(1.6,0.15)";}
  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="#tilde{g}_{1600}#rightarrow b#bar{b}#tilde{#chi}_{1,150}^{0}";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName="T5bbbb_ZG_1550";}

  // const char *l_name=lName.c_str();
  if(i<nBG)legend1->AddEntry(hist,lName,"f");
  else legend1->AddEntry(hist,lName,"l");
  // legend1->SetTextSize(0.04);
}


void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
  
  if(lastBinCt!=0 && overflCt!=0)
    //    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );
    lastBinErr = sqrt( (lastBinErr*lastBinErr) + (overflErr*overflErr) );  

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
    
}
