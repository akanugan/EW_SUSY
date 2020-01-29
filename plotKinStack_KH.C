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
#include "TH1D.h"
#include"TMath.h"

using namespace std;

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=9,nBG=6;    //Specify no. of files
TFile *f[nfiles];
int col[11]={kPink+1,kTeal+9,kGreen,kYellow,kOrange,kBlue,kCyan,kRed,kBlue+2,kMagenta,kPink+1};  //Specify Colors b's

TCanvas *c_cA=new TCanvas("kinVar","plot of a kin var",1200,1200);

void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);

TLegend *legend1=new TLegend(0.4501, 0.7,  0.93, 0.88);
//TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);
//TLegend *legend1=new TLegend(0.1,0.7,0.48,0.9);

//TLegend *legend2=new TLegend(0.7, 0.9,  0.90, 0.65);
//TLegend *legend2=new TLegend(0.6, 0.90,  0.98, 0.45);
void setLastBinAsOverFlow(TH1D*);
void plotKinStack_KH(){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  TString varName = "AK8Mass";
  TString xLabel = "Softdrop mass AK8jet1 (GeV) ";
  int rebin=1;

  f[0] = new TFile("BKGhistos_allcuts/ST__MC2018.root");
  f[1] = new TFile("BKGhistos_allcuts/Rare_MC2018.root");
  f[3] = new TFile("BKGhistos_allcuts/QCD_HT_MC2018.root");
  f[2] = new TFile("BKGhistos_allcuts/TTJets_MC2018.root");
  f[4] = new TFile("BKGhistos_allcuts/WJetsToLNu_HT_MC2018.root");
  f[5] = new TFile("BKGhistos_allcuts/ZJetsToNuNu_HT_MC2018.root");
  f[6] = new TFile("SIG_WZ_histos_baseline/TChiWZ_600_100_MC2018.root");  
  f[7] = new TFile("SIG_WZ_histos_baseline/TChiWZ_800_100_MC2018.root");
  f[8] = new TFile("SIG_WZ_histos_baseline/TChiWZ_1000_100_MC2018.root");

  // f[6] = new TFile("SIG_WW_tightmass_vetoB/TChipmWW_600_100_MC2018.root");  
  // f[7] = new TFile("SIG_WW_tightmass_vetoB/TChipmWW_800_100_MC2018.root");
  // f[8] = new TFile("SIG_WW_tightmass_vetoB/TChipmWW_1000_100_MC2018.root");

  TH1D *h_Sig=(TH1D*)f[6]->FindObjectAny(varName);
  TH1D *h_Sig1=(TH1D*)f[7]->FindObjectAny(varName);
  TH1D *h_Sig2=(TH1D*)f[8]->FindObjectAny(varName);

  //
  // Top panel
  //
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.0, 1.0, 1.0);
  pad1->SetTopMargin(0.12);
  pad1->SetBottomMargin(0.30); // Upper and lower plot are joined
  pad1->SetLeftMargin(0.1);
  pad1->SetRightMargin(0.05);
  pad1->Draw();
  pad1->cd();

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
  for(int i=0;i<nfiles;i++){
        
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    
    h_MET->Rebin(rebin);
    //    h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    //    h_MET->SetMinimum(100);
    decorate(h_MET,i,f[i]->GetName());
    
    if(i<=(nBG-1)) hs_var->Add(h_MET);

    if(i==nBG-1) {
      //c_cA->cd();
      hs_var->Draw("BAR HIST");
      hs_var->Draw("HIST");
      hs_var->SetMinimum(0.8);
      hs_var->SetMaximum(hs_var->GetMaximum()*10);
      decorate(hs_var,i,f[i]->GetName()); 
      //hs_var->GetYaxis()->SetRangeUser(100.5,20000);
    }
    if(i>=nBG){ 
      //c_cA->cd(); 
      h_MET->SetMarkerStyle(20);
      h_MET->SetMarkerColor(col[i]);
      h_MET->SetLineColor(col[i]);
      h_MET->SetLineWidth(3);
      h_MET->Draw("hist same");
      //      h_MET->GetYaxis()->SetRangeUser(0.5,20000);
      //      h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    }
    drawlegend(h_MET,i,f[i]->GetName());
    if(i==nfiles-1){
      hs_var->GetXaxis()->SetRangeUser(0,2100); 
      hs_var->GetXaxis()->SetTitleOffset(.90);
      //hs_var->GetXaxis()->SetTitle(xLabel); hs_var->GetYaxis()->SetTitle("Events");hs_var->SetTitle(0);
    }
  }
  
  // AR
  //TH1D *hRatio = (TH1D*)h_Sig ->Clone();
  //TH1* hRatio = new TH1I("h1", "h1 title", 100, 0.0, 4.0);
  TH1 *stack_sum = static_cast<TH1*>(hs_var->GetStack()->Last());
  TH1 *sstack_sum = (TH1D*)h_Sig ->Clone(); // just to have same bin content
  TH1 *sigbkg_sum = (TH1D*)h_Sig ->Clone(); // just to have same bin content
  TH1 *sigbkg_sum1  = (TH1D*)h_Sig1 ->Clone(); // just to have same bin content
  TH1 *sigbkg_sum2  = (TH1D*)h_Sig2 ->Clone(); // just to have same bin content
  if (stack_sum->Integral() !=0){
  for (int bin=0;bin<=stack_sum->GetNcells();++bin) {
    sstack_sum->SetBinContent(bin,sqrt(stack_sum->GetBinContent(bin)));
    sigbkg_sum->SetBinContent(bin,sqrt(stack_sum->GetBinContent(bin)+h_Sig->GetBinContent(bin)));
    sigbkg_sum1->SetBinContent(bin,sqrt(stack_sum->GetBinContent(bin)+h_Sig1->GetBinContent(bin)));
    sigbkg_sum2->SetBinContent(bin,sqrt(stack_sum->GetBinContent(bin)+h_Sig2->GetBinContent(bin)));
  }
  }
  //TRatioPlot *hRatio = new TRatioPlot(hratio,stack_sum,"Ratio");
  //printf("%f\n",static_cast<float>(hRatio->Integral()));

  sigbkg_sum->Add(sstack_sum,-1);
  sigbkg_sum->Scale(2);
  sigbkg_sum1->Add(sstack_sum,-1);
  sigbkg_sum1->Scale(2);
  sigbkg_sum2->Add(sstack_sum,-1);
  sigbkg_sum2->Scale(2);

  //Option1:  Q = 2(sqrt(S+B)-sqrt(B))
  TH1D *hRatio = (TH1D*)sigbkg_sum ->Clone();
  TH1D *hRatio1 = (TH1D*)sigbkg_sum1 ->Clone();
  TH1D *hRatio2 = (TH1D*)sigbkg_sum2 ->Clone();
  //

  //Option2: Q = S/sqrt(B)
  // TH1D *hRatio = (TH1D*)h_Sig ->Clone();
  // TH1D *hRatio1 = (TH1D*)h_Sig1 ->Clone();
  // TH1D *hRatio2 = (TH1D*)h_Sig2 ->Clone();
  // hRatio->Divide(sstack_sum);
  // hRatio1->Divide(sstack_sum);
  // hRatio2->Divide(sstack_sum);
  // 
 
  printf("%f\n",static_cast<float>(hRatio->Integral()));
  //hRatio->Divide(static_cast<TH1*>(hs_var->GetStack()->Last()));
  //printf("%f\n",static_cast<float>(static_cast<TH1*>(hs_var->GetStack()->Last())->Integral()));

  hs_var->Draw("BAR HIST");
  hs_var->Draw("HIST");
  h_Sig->Draw("hist same");
  h_Sig1->Draw("hist same");
  h_Sig2->Draw("hist same");
  //hs_var->SetTitle("");
  hs_var->GetXaxis()->CenterTitle(true);
  hs_var->GetYaxis()->CenterTitle(true);
  hs_var->GetYaxis()->SetTitle("# of events");
  //hs_var->GetXaxis()->SetTitle(xLabel);
  hs_var->GetYaxis()->SetTitleOffset(1.0);
  hs_var->GetYaxis()->SetTitleFont(42);
  hs_var->GetYaxis()->SetTitleSize(0.045);
  hs_var->GetYaxis()->SetLabelSize(0.035);
  hs_var->GetYaxis()->SetLabelFont(42);
      
  hs_var->GetXaxis()->SetTitleOffset(0.65);
  hs_var->GetXaxis()->SetTitleFont(42);
  hs_var->GetXaxis()->SetLabelSize(0.02);
  hs_var->GetXaxis()->SetLabelFont(42);
  hs_var->GetXaxis()->SetTitleSize(0.035);
  pad1->SetLogy();
  //h_v17->SetMarkerColor(kRed);
  //hs_var->SetMarkerColor(kBlue);
  //hs_var->Draw("hist");
  //h_MET->Draw("hist same");

  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  //c_cA->cd();
  gPad->SetLogy();legend1->Draw();
  //  gPad->RedrawAxis();
  //  hs_var->GetXaxis()->SetTitle(xLabel);
  
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.68,0.91,name2);
  
  //
  // Bottom panel
  //
  //// Ratio Plot in lower panel

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.02, 1.0, 0.30);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.4);
  pad2->SetLeftMargin(0.1);
  pad2->SetRightMargin(0.05);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  hRatio->GetXaxis()->SetTitle(xLabel);
  hRatio->GetXaxis()->SetTitleFont(42);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetXaxis()->CenterTitle(true);
  hRatio->GetXaxis()->SetTitleSize(0.19);  
  hRatio->GetYaxis()->CenterTitle(true);
  hRatio->GetYaxis()->SetTitle("2(#sqrt{S+B}-#sqrt{B})");
  hRatio->GetYaxis()->SetTitleOffset(0.25);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetTitleFont(42);
  hRatio->GetYaxis()->SetLabelFont(42);
  hRatio->GetXaxis()->SetRangeUser(0,2100); //set to same as pad 1

  hRatio->SetMarkerSize(1.);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineWidth(3);
  hRatio->SetMarkerStyle(20);
  hRatio->SetAxisRange(0.0,5,"Y");
  hRatio->GetXaxis()->SetLabelSize(0.07);
  hRatio->GetYaxis()->SetLabelSize(0.12);

  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTickLength(0.025);
  hRatio->GetXaxis()->SetTickLength(0.08);
  
  //hRatio->SetLineWidth(5);
  hRatio->SetMaximum(4.9);
  hRatio->Draw("hist");
  hRatio1->Draw("hist same");
  hRatio2->Draw("hist same");
  Double_t Gmax=hRatio->GetXaxis()->GetBinUpEdge(hRatio->GetNbinsX()-5);
  Double_t Gmin=hRatio->GetXaxis()->GetXmin();
  printf("%f\n",Gmax);
  TLine *line = new TLine(Gmin,2,2100,2);
  line->SetLineColor(kBlue);
  line->Draw();
  //

  name3 = varName+".pdf";
  c_cA->SaveAs(name3);
  cout<<"*****************************************************"<<endl;
  cout<<"Int Lumi(inv.fb) for file1:"<<setprecision(4)<<intLumi<<endl;}


void decorate(THStack *hs,int i,const char* fname){
  //  hs->SetMinimum(0.5);
  //hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.05);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
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
  else if(lName.Contains("TChiWZ_600_100")){lName="TChiWZ_600_100";}
  else if(lName.Contains("TChiWZ_800_100")){lName="TChiWZ_800_100";}
  else if(lName.Contains("TChiWZ_1000_100")){lName="TChiWZ_1000_100";}

  else if(lName.Contains("TChiWZ_300_1")){lName="TChiWZ_300_1";}
  else if(lName.Contains("TChiWZ_400_1")){lName="TChiWZ_400_1";}
  else if(lName.Contains("TChiWZ_500_1")){lName="TChiWZ_500_1";}

  else if(lName.Contains("TChipmWW_600_100")){lName="TChiWW_600_100";}
  else if(lName.Contains("TChipmWW_800_100")){lName="TChiWw_800_100";}
  else if(lName.Contains("TChipmWW_1000_100")){lName="TChiWW_1000_100";}

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
