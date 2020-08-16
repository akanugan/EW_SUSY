
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

TLatex textOnTop,intLumiE;
double intLumi=137;
TString name;
bool saveCanvas=0;

void ROCplotter(){
  TFile *fSig = TFile::Open("T2tt_MC2016_fastsim.root");
  //  TFile *fBG  = TFile::Open("ZJetsToNuNu_HT_MC2018.root");
  //TFile *fBG  = TFile::Open("WJetsToLNu_HT_MC2018.root");
  //  TFile *fBG  = TFile::Open("QCD_HT_MC2018.root");
  TFile *fBG  = TFile::Open("ROC_Wmatch_v3/Bkg/QCD_HT_LDP_MC2016.root");

  vector<TString> sigHist = {"wAK8J1wDis","wAK8J1wDisDC","wAK8J1Tau21"};
  vector<TString> bgHist  = {"wAK8J1wDis","wAK8J1wDisDC","wAK8J1Tau21"};
  vector<TString> legName = {"Deep W",  "Deep W (MD)",   "#tau_{21}"};
  vector<int> col = {kRed,kBlue,kCyan};

  vector<double> cutValues = {0.918, 0.704, 0.35};//for each ROC curve, the point corresponding to this values will be shown.
  //  vector<double> cutValues = {0.918, 0.245, 0.35};//for each ROC curve, the point corresponding to this values will be shown.
  //  vector<double> cutValues = {0.918, 0.479, 0.45};//for each ROC curve, the point corresponding to this values will be shown.
  
  TMarker *mark[sigHist.size()];
  
  TH1D *hS[sigHist.size()], *hBG[bgHist.size()];
  vector<float> sigEff[sigHist.size()], bgEff[bgHist.size()];
  
  TCanvas *roc = new TCanvas("rocPlot","ROC plot",1200,800);
  TLegend *leg = new TLegend(0.65, 0.78,  0.85, 0.55);
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);
  
  for(int i=0;i<sigHist.size();i++){
    fSig->GetObject(sigHist[i],hS[i]);
    
    vector<float> vec;
    float total = hS[i]->Integral();
    if(total < 0.000001){cout<<"0 entry error in getting "<<sigHist[i]<<" from file "<<fSig<<endl; return;}
    
    for(int j=1;j<=hS[i]->GetNbinsX()+1;j++){
      if(sigHist[i].Contains("Tau21")) vec.push_back(hS[i]->Integral(1,j)/total);
      else vec.push_back(hS[i]->Integral(j,hS[i]->GetNbinsX()+1)/total);
    }
    sigEff[i] = vec;
  }
  for(int i=0;i<bgHist.size();i++){
    fBG->GetObject(bgHist[i],hBG[i]);
    vector<float> vec;
    float total = hBG[i]->Integral();
    if(total < 0.000001){cout<<"0 entry error in getting "<<bgHist[i]<<" from file "<<fBG<<endl; return;}
    
    for(int j=1;j<=hBG[i]->GetNbinsX()+1;j++){
       if(bgHist[i].Contains("Tau21")) vec.push_back(hBG[i]->Integral(1,j)/total);
       else vec.push_back(hBG[i]->Integral(j,hBG[i]->GetNbinsX()+1)/total);
    }
    bgEff[i] = vec;
  }

  //-----------------
  TGraph *gr[sigHist.size()];
  for(int i=0;i<sigHist.size();i++){
    gr[i] = new TGraph(bgEff[i].size(),&(bgEff[i][0]),&(sigEff[i][0]));
    //    gr[i] = new TGraph(bgEff[i].size(),&(sigEff[i][0]),&(bgEff[i][0]));
    gr[i]->SetLineColor(col[i]);
    gr[i]->SetMarkerColor(col[i]);
    gr[i]->SetLineWidth(2);
    leg->AddEntry(gr[i],legName[i],"l");

    //gr[i]->SetTitle(";BG efficiency (Z(#nu#nu));Signal efficiency (t#bar{t})");
    //    gr[i]->SetTitle(";BG efficiency (QCD);Signal efficiency (t#bar{t})");
    gr[i]->SetTitle(";BG efficiency (t#bar{t});Signal efficiency (t#bar{t})");
    
    if(i==0) gr[i]->Draw();
    else gr[i]->Draw("same");
    double x,y;
    if(sigHist[i].Contains("Tau21")){
      x = hBG[i]->Integral(1,hBG[i]->FindBin(cutValues[i]))/hBG[i]->Integral();
      y = hS[i]->Integral(1,hS[i]->FindBin(cutValues[i]))/hS[i]->Integral();
    }
    else{
      x = hBG[i]->Integral(hBG[i]->FindBin(cutValues[i]),hBG[i]->GetNbinsX()+1)/hBG[i]->Integral();
      y = hS[i]->Integral(hS[i]->FindBin(cutValues[i]),hS[i]->GetNbinsX()+1)/hS[i]->Integral();
    }
    cout<<legName[i]<<" | WP = "<<cutValues[i]<<" | SigEff = "<<y<<" | BGEff = "<<x<<endl;
    mark[i] = new TMarker(x,y,21);
    mark[i]->Draw();
  }
  leg->Draw();
  // for(int i=0;i<sigHist.size();i++){
  //   cout<<sigHist[i]<<endl;
  //   for(int j=0;j<sigEff[i].size();j++)
  //     cout<<sigEff[i][j]<<" ";
  // }
  
}
