#define SignalReg_cxx
#include "SignalReg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  SignalReg ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void SignalReg::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  TTree *newtree = fChain->CloneTree(0);
  
  TFile *fXsec = new TFile("/uscms/home/vhegde/nobackup/Physics/skims_Run2ProductionV12/MakeSkims/SignalRegion/SignalSkimming/TChiWG_MassScan.root");
  TH1D *h1=(TH1D*)fXsec->Get("mEWeakino_N2C1Xsec");
  
  float xsec = 0.0, mass=0.;
  if(s_data.Contains("TChiWZ_1000")) mass = 1000.0;
  else if(s_data.Contains("TChiWZ_900")) mass = 900.0;
  else if(s_data.Contains("TChiWZ_800")) mass = 800.0;
  else if(s_data.Contains("TChiWZ_700")) mass = 700.0;
  else if(s_data.Contains("TChiWZ_600")) mass = 600.0;
  
  xsec = h1->GetBinContent(h1->FindBin(mass));
  cout<<"Assigning xsec as: "<<xsec<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    CrossSection = xsec;
    NumEvents = nentries;
    Weight = xsec/NumEvents;
    newtree->Fill();
  }
}
