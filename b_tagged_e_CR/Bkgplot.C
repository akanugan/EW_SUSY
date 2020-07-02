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

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=1,nBG=1;    //Specify no. of files
TFile *f[nfiles];
int col[14]={kBlue,kOrange,kPink+1,kTeal+9,kGreen,kYellow,kCyan,kRed,kBlue+2,kMagenta,kPink+1,kOrange+1,kGreen+3,kBlue+7};  //Specify Colors b's

TCanvas *c_cA=new TCanvas("kinVar","",1500,900);

void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
//TLegend *legend1=new TLegend(0.4501, 0.65,  0.88, 0.88);
TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);
//TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);

//TLegend *legend2=new TLegend(0.7, 0.9,  0.90, 0.65);
//TLegend *legend2=new TLegend(0.6, 0.90,  0.98, 0.45);
void setLastBinAsOverFlow(TH1D*);
void Bkgplot(){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  //TString varName = "AK8J2wDis";
  //  TString varName = "AK8J2wDisnoBtag";
  //  TString varName = "AK8J1wDiscase2";
  TString varName = "wzMETvBin_RegD";
  TString xLabel = "MET in Reg D (GeV)";
  int rebin=1;

  //f[0] = new TFile("BKGhistos_allcuts/ST__MC2018.root");
  // f[1] = new TFile("BKGhistos_allcuts/Rare_MC2018.root");
  // f[3] = new TFile("BKGhistos_allcuts/QCD_HT_MC2018.root");
  // f[2] = new TFile("BKGhistos_allcuts/TTJets_MC2018.root");
  // f[4] = new TFile("BKGhistos_allcuts/WJetsToLNu_HT_MC2018.root");
  // f[5] = new TFile("BKGhistos_allcuts/ZJetsToNuNu_HT_MC2018.root");
  f[0] = new TFile("WJetsToLNu_HT_MC2018.root");

  // f[6] = new TFile("TChiWZ_600_100_MC2018.root");
  // f[7] = new TFile("TChiWZ_800_100_MC2018.root");
  // f[8] = new TFile("TChiWZ_1000_100_MC2018.root");
  //  f[9] = new TFile("TChiWZ_800_500_MC2018.root");
  //f[10] = new TFile("TChiWZ_800_700_MC2018.root");
  //  f[12] = new TFile("TChiWZ_900_100_MC2018.root");
  //f[11] = new TFile("TChiWZ_1000_1_MC2018.root");
  


  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","MET Stacked");
  //TH1D *h_R;
  TH1D *h_MET_R[nfiles];
  for(int i=0;i<nfiles;i++){
    sprintf(name,"hist_file%i",i);
    h_MET_R[i]=new TH1D(name,name,21,0.5,21.5);
  }
  vector<double> Bcnt;
  double intLumi=137;
  TLatex tl1;
  for(int i=0;i<nfiles;i++){
        
    printf("%d\n",i);
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    //    h_MET->Rebin(rebin);
    //    h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    //    h_MET->SetMinimum(100);
    decorate(h_MET,i,f[i]->GetName());
    
    if(i<=(nBG-1))  hs_var->Add(h_MET);

    if(i==nBG-1) {
      c_cA->cd();
      hs_var->Draw("BAR HIST");
      //hs_var->Draw("colz");
      hs_var->Draw("HIST");
      hs_var->SetMinimum(0.8);
      hs_var->SetMaximum(hs_var->GetMaximum()*10);
      decorate(hs_var,i,f[i]->GetName()); 
      //hs_var->GetYaxis()->SetRangeUser(100.5,20000);
    }
    if(i>=nBG){ 
      c_cA->cd(); 
      h_MET->SetMarkerStyle(20);
      h_MET->SetMarkerColor(col[i]);
      h_MET->SetLineColor(col[i]);
      h_MET->SetLineWidth(3);
      h_MET->Draw("hist same");
      //h_MET->Draw("colz");
      //      h_MET->GetYaxis()->SetRangeUser(0.5,20000);
      //      h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    }
    drawlegend(h_MET,i,f[i]->GetName());
    if(i==nfiles-1){
      //hs_var->GetXaxis()->SetRangeUser(0,1500); 
      //hs_var->GetXaxis()->SetTitleOffset(.90);
      hs_var->GetXaxis()->SetTitle(xLabel); hs_var->GetYaxis()->SetTitle("Events");hs_var->SetTitle(0);
    }
  }
  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  c_cA->cd();
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
  //  gPad->Update();
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
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );
  
  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
    
}
