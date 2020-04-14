#ifndef SignalReg_H
#define SignalReg_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include <vector>

class SignalReg : public NtupleVariables{

 public:
  SignalReg(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~SignalReg();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  Double_t find_bjets_mtbmin();
  void print(Long64_t);

  //Variables defined
  bool isMC=true;
  double wt=0,lumiInfb=35.815165;
  double deepCSVvalue = 0;
  double massLowW = 65., massHighW = 90.; //65-90, 55-100
  double massLowZ = 65., massHighZ = 90.; //65-90, 55-100
  double massLowH = 85., massHighH = 135.;
  double bbscore = 0.3; //  double b score
  double deepbbscore = 0.7; // deep double b score
  double deepAK8Wscore = 0.9; // deepAK8W score
  vector<TLorentzVector> bjets;

  TH1D *h_filters;
  //  TH1D *h_MET;
  //  TH1D *h_MHT;
  TH1D *h_HT;
  TH1D *h_madHT;

  /* TH1D *h_NJets; */
  /* TH1D *h_BTags; */

  //TH1D *h_MT;
  //TH1D *h_MT2J;
  /* TH1D *h_dPhiMETAK8; */
  /* TH1D *h_dPhiAK8J1J2; */

  /* TH1D *h_AK8J1Pt, *h_AK8J1Mass, *h_AK8J1Eta, *h_AK8J1Tau21, *h_AK8J1wDisDC, *h_AK8J1zhDisDC; */
  /* TH1D *h_AK8J2Pt, *h_AK8J2Mass, *h_AK8J2Eta, *h_AK8J2Tau21, *h_AK8J2wDisDC, *h_AK8J2zhDisDC; */
  /* TH1D *h_AK8J1Mass1, *h_AK8J2Mass1; */
  /* TH1D *h_AK8J1Mass2, *h_AK8J2Mass2; */
  /* TH1D *h_AK8J1Mass3, *h_AK8J2Mass3; */
  /* TH1D *h_AK8J1Mass4, *h_AK8J2Mass4; */

  //WH SR
  TH1D *h_WHAK8J1Pt, *h_WHAK8J1Mass, *h_WHAK8J1Eta, *h_WHAK8J1Tau21, *h_WHAK8J1wDis;
  TH1D *h_WHAK8J2Pt, *h_WHAK8J2Mass, *h_WHAK8J2Eta, *h_WHAK8J2Tau21, *h_WHAK8J2wDis;
  TH1D *h_WHMET;
  TH1D *h_WHMT;
  TH1D *h_WHMT2J;
  TH1D *h_WHMETa;
  TH1D *h_WHMETc;

  TH1D *h_HWAK8J1Pt, *h_HWAK8J1Mass, *h_HWAK8J1Eta, *h_HWAK8J1Tau21, *h_HWAK8J1wDis;
  TH1D *h_HWAK8J2Pt, *h_HWAK8J2Mass, *h_HWAK8J2Eta, *h_HWAK8J2Tau21, *h_HWAK8J2wDis;
  TH1D *h_HWMET;
  TH1D *h_HWMT;
  TH1D *h_HWMT2J;
  TH1D *h_HWMETa;
  TH1D *h_HWMETc;
  
  TH1D *h_WHmtbmin, *h_WHmct;
  TH1D *h_HWmtbmin, *h_HWmct;
  TH2D *h2_WHmtbminMct, *h2_WHmtbminHmass;
  TH2D *h2_HWmtbminMct, *h2_HWmtbminHmass;

  // WW SR
  TH1D *h_WWMET;
  TH1D *h_WWMT;
  TH1D *h_WWMT2J;
  TH1D *h_WWAK8J1Pt;
  TH1D *h_WWAK8J1Eta;
  TH1D *h_WWAK8J1Mass;

  //Bkest_MET
  TH1D *h_WHMET_RegA, *h_WHMET_RegB, *h_WHMET_RegC, *h_WHMET_RegD;
  TH1D *h_WHAK8J2Mass_RegA, *h_WHAK8J2Mass_RegB, *h_WHAK8J2Mass_RegC, *h_WHAK8J2Mass_RegD;
  TH1D *h_HWMET_RegA, *h_HWMET_RegB, *h_HWMET_RegC, *h_HWMET_RegD;
  TH1D *h_HWAK8J1Mass_RegA, *h_HWAK8J1Mass_RegB, *h_HWAK8J1Mass_RegC, *h_HWAK8J1Mass_RegD;  

  // Mass SB
  TH1D *h_WHAK8J1MassSB, *h_WHAK8J1MassNo2bTag;
  TH1D *h_WHAK8J2MassSB, *h_WHAK8J2MassNo2bTag;
  TH1D *h_HWAK8J1MassSB, *h_HWAK8J1MassNo2bTag;
  TH1D *h_HWAK8J2MassSB, *h_HWAK8J2MassNo2bTag;

  // GenReco match
  TH1D *h_AK8J1doubleBDis, *h_AK8J2doubleBDis;
  TH1D *h_AK8J1deepdoubleBDis, *h_AK8J2deepdoubleBDis;
  TH1D *h_AK8J1Tau21, *h_AK8J2Tau21;
  TH1D *h_AK8J1zDis, * h_AK8J2zDis;
  TH1D *h_AK8J1deepdoubleBDisQ, *h_AK8J2deepdoubleBDisQ;
  TH1D *h_AK8J1zhDisMD, *h_AK8J2zhDisMD;

  TH1D *h_mtbmin, *h_mct;
  /* TH2D *h2_AK8J1J2Tau21; */
  /* TH2D *h2_DisxodRAK8J1; */
  /* TH2D *h2_DisdRAK8J2; */
  /* TH2D *h2_Tau21dRAK8J1; */
  /* TH2D *h2_Tau21dRAK8J2; */
  /* TH2D *h2_AK8J1Mass_J1Tau21; */
  /* TH2D *h2_AK8J2Mass_J2Tau21; */
  
  /* TH1D *h_dPhi1; */
  /* TH1D *h_dPhi2; */
  /* TH1D *h_dPhi3; */
  /* TH1D *h_dPhi4; */

  TH1F *h_cutflow;
  TFile *oFile;
  
};
#endif

#ifdef SignalReg_cxx

void SignalReg::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  TString name,title;
 
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);

  h_cutflow = new TH1F("CutFlow","cut flow",25,0,25);
  h_filters = new TH1D("Filters","Filters: Bin1 : all nEvnts, other bins: filter pass/fail",10,0,10);
  
  /* h_MET = new TH1D("MET","MET",200,0,2000); */
  /* h_MHT = new TH1D("MHT","MHT",200,0,2000); */
   h_HT = new TH1D("HT","HT",100,0,5000); 
   h_madHT = new TH1D("madHT","madHT",100,0,5000); 

  /* h_NJets = new TH1D("NJets","NJets with pT > 30, |eta| < 20.4",20,0,20);   */
  /* h_BTags = new TH1D("BTags","BTags with DeepCSV MedWP",10,0,10);   */
  
  //h_MT = new TH1D("mT","mT(MET,AK8J)",200,0,2000);
  // h_MT2J = new TH1D("mT2J","mT(MET,AK8J2)",200,0,2000);  
  //h_dPhiMETAK8 = new TH1D("dPhiMETAK8","dPhi(MET,AK8J)",40,0,4);
  /* h_dRGenAK8J1 = new TH1D("dRGenAK8J1","deltaR(WGen,AK8J1)",40,0,4); */
  /* h_dRGenAK8J2 = new TH1D("dRGenAK8J2","deltaR(WGen,AK8J2)",40,0,4); */
  //  h_dPhiAK8J1J2 = new TH1D("dPhiAK8J1J2","dPhi(AK8J1,AK8J2)",40,0,4);

  /* h_AK8J1Pt = new TH1D("AK8Pt","Leading AK8 jets Pt",200,0,2000); */
  /* h_AK8J1Eta = new TH1D("AK8Eta","AK8 Eta",120,-6,6); */
  /* h_AK8J1Mass = new TH1D("AK8Mass","AK8 Mass",60,0,300); */
  /* h_AK8J1Tau21 = new TH1D("AK8Tau21","AK8 Tau21",100,0,1); */
  /* h_AK8J1wDisDC = new TH1D("AK8J1wDisDC","AK8 J1 w Discr.DeepDecorelated",100,0,1); */
  /* h_AK8J1zhDisDC = new TH1D("AK8J1zhDisDC","AK8 J1 zh Discr.DeepDecorelated",100,0,1); */
  
  //h_AK8J1wDis = new TH1D("AK8J1wDis","AK8 J1 w Discr.Deep Corelated",100,0,1);
  //h_AK8J1zhDis = new TH1D("AK8J1zhDis","AK8 J1 zh Discr.Deep Corelated",20,0,1);
  //h_AK8J1zDis = new TH1D("AK8J1zDis","AK8 J1 z Discr.Deep Corelated",100,0,1);


  /* h_AK8J2Pt = new TH1D("AK8J2Pt","2nd leading AK8 jets Pt",200,0,2000); */
  /* h_AK8J2Eta = new TH1D("AK8J2Eta","AK8J2 Eta",120,-6,6); */
  /* h_AK8J2Mass = new TH1D("AK8J2Mass","AK8J2 Mass",60,0,300); */
  /* h_AK8J2Tau21 = new TH1D("AK8J2Tau21","AK8J2 Tau21",100,0,1); */
  /* h_AK8J2wDisDC = new TH1D("h_AK8J2wDisDC","AK8 J2 w Discr.DeepDecorelated",100,0,1); */
  /* h_AK8J2zhDisDC = new TH1D("h_AK8J2zhDisDC","AK8 J2 zh Discr.DeepDecorelated",100,0,1); */
  
  // h_AK8J2wDis = new TH1D("AK8J2wDis","AK8 J2 w Discr.Deep Corelated",100,0,1);
  //h_AK8J2zhDis = new TH1D("AK8J2zhDis","AK8 J2 zh Discr.Deep Corelated",20,0,1);
  //h_AK8J2zDis = new TH1D("AK8J2zDis","AK8 J2 z Discr.Deep Corelated",100,0,1);
  //other mass histos
  /* h_AK8J1Mass2 = new TH1D("AK8J1Mass2","AK8j1 Mass2",60,0,300);   */
  /* h_AK8J2Mass2 = new TH1D("AK8J2Mass2","AK8j2 Mass2",60,0,300);   */
  /* h_AK8J1Mass3 = new TH1D("AK8J1Mass3","AK8j1 Mass3",60,0,300);   */
  /* h_AK8J2Mass3 = new TH1D("AK8J2Mass3","AK8j2 Mass3",60,0,300);   */
  /* h_AK8J1Mass4 = new TH1D("AK8J1Mass4","AK8j1 Mass4",60,0,300);   */
  /* h_AK8J2Mass4 = new TH1D("AK8J2Mass4","AK8j2 Mass4",60,0,300);   */

  // ---for wh----
  h_WHAK8J1Pt = new TH1D("WHAK8J1Pt","leading AK8 jets Pt",200,0,2000);
  h_WHAK8J1Eta = new TH1D("WHAK8J1Eta","AK8J1 Eta",120,-6,6);
  h_WHAK8J1Mass = new TH1D("WHAK8J1Mass","AK8J1 Mass",60,0,300);

  h_WHAK8J1Tau21 = new TH1D("WHAK8J1Tau21","AK8J1 Tau21",100,0,1);
  h_WHAK8J1wDis = new TH1D("WH_AK8J1wDis","AK8 J1 w Discr. corelated",100,0,1);
  h_WHMET = new TH1D("WHMET","MET",200,0,2000);
  h_WHMT = new TH1D("WHmT","mT(MET,AK8J)",200,0,2000);
  h_WHMT2J = new TH1D("WHmT2J","mT(MET,AK8J2)",200,0,2000);  
  h_WHMETa = new TH1D("WHMETa","MET",200,0,2000);
  h_WHMETc = new TH1D("WHMETc","MET",200,0,2000);

  h_WHAK8J2Pt = new TH1D("WHAK8J2Pt","2nd leading AK8 jets Pt",200,0,2000);
  h_WHAK8J2Eta = new TH1D("WHAK8J2Eta","AK8J2 Eta",120,-6,6);
  h_WHAK8J2Mass = new TH1D("WHAK8J2Mass","AK8J2 Mass",60,0,300);

  h_WHAK8J2Tau21 = new TH1D("WHAK8J2Tau21","AK8J2 Tau21",100,0,1);
  h_WHAK8J2wDis = new TH1D("WH_AK8J2wDis","AK8 J2 w Discr. corelated",100,0,1);
    
  // --- for hw
  h_HWAK8J1Pt = new TH1D("HWAK8J1Pt"," leading AK8 jets Pt",200,0,2000);
  h_HWAK8J1Eta = new TH1D("HWAK8J1Eta","AK8J1 Eta",120,-6,6);
  h_HWAK8J1Mass = new TH1D("HWAK8J1Mass","AK8J1 Mass",60,0,300);

  h_HWAK8J1Tau21 = new TH1D("HWAK8J1Tau21","AK8J1 Tau21",100,0,1);
  h_HWAK8J1wDis = new TH1D("HWAK8J1wDis","AK8 J2 w Discr.corelated",100,0,1);
  
  h_HWMET = new TH1D("HWMET","MET",200,0,2000);
  h_HWMT = new TH1D("HWmT","mT(MET,AK8J)",200,0,2000);
  h_HWMT2J = new TH1D("HWmT2J","mT(MET,AK8J2)",200,0,2000);  
  h_HWMETa = new TH1D("HWMETa","MET",200,0,2000);
  h_HWMETc = new TH1D("HWMETc","MET",200,0,2000);

  h_HWAK8J2Pt = new TH1D("HWAK8J2Pt","2nd leading AK8 jets Pt",200,0,2000);
  h_HWAK8J2Eta = new TH1D("HWAK8J2Eta","AK8J2 Eta",120,-6,6);
  h_HWAK8J2Mass = new TH1D("HWAK8J2Mass","AK8J2 Mass",60,0,300);
 
  h_HWAK8J2Tau21 = new TH1D("HWAK8J2Tau21","AK8J2 Tau21",100,0,1);
  h_HWAK8J2wDis = new TH1D("HWAK8J2wDis","AK8 J2 w Discr. corelated",100,0,1);
 
  // WW SR
  h_WWMET = new TH1D("WWMET","MET",200,0,2000);
  h_WWMT = new TH1D("WWmT","mT(MET,AK8J)",200,0,2000);
  h_WWMT2J = new TH1D("WWmT2J","mT(MET,AK8J2)",200,0,2000);  
  h_WWAK8J1Pt = new TH1D("WWAK8J1Pt"," leading AK8 jets Pt",200,0,2000);
  h_WWAK8J1Eta = new TH1D("WWAK8J1Eta","AK8J1 Eta",120,-6,6);
  h_WWAK8J1Mass = new TH1D("WWAK8J1Mass","AK8J1 Mass",60,0,300);


 // Bkgest MET
  h_WHMET_RegA = new TH1D("WHMET_RegA","WHMETA",200,0,2000);
  h_WHMET_RegB = new TH1D("WHMET_RegB","WHMETB",200,0,2000);
  h_WHMET_RegC = new TH1D("WHMET_RegC","WHMETC",200,0,2000);
  h_WHMET_RegD = new TH1D("WHMET_RegD","WHMETD",200,0,2000);

  h_WHAK8J2Mass_RegA = new TH1D("WHAK8J2Mass_RegA","WHAK8J2MassA",60,0,300);
  h_WHAK8J2Mass_RegB = new TH1D("WHAK8J2Mass_RegB","WHAK8J2MassB",60,0,300);
  h_WHAK8J2Mass_RegC = new TH1D("WHAK8J2Mass_RegC","WHAK8J2MassC",60,0,300);
  h_WHAK8J2Mass_RegD = new TH1D("WHAK8J2Mass_RegD","WHAK8J2MassD",60,0,300);

  h_HWMET_RegA = new TH1D("HWMET_RegA","HWMETA",200,0,2000);
  h_HWMET_RegB = new TH1D("HWMET_RegB","HWMETB",200,0,2000);
  h_HWMET_RegC = new TH1D("HWMET_RegC","HWMETC",200,0,2000);
  h_HWMET_RegD = new TH1D("HWMET_RegD","HWMETD",200,0,2000);

  h_HWAK8J1Mass_RegA = new TH1D("HWAK8J1Mass_RegA","HWAK8J1MassA",60,0,300);
  h_HWAK8J1Mass_RegB = new TH1D("HWAK8J1Mass_RegB","HWAK8J1MassB",60,0,300);
  h_HWAK8J1Mass_RegC = new TH1D("HWAK8J1Mass_RegC","HWAK8J1MassC",60,0,300);
  h_HWAK8J1Mass_RegD = new TH1D("HWAK8J1Mass_RegD","HWAK8J1MassD",60,0,300);
  //

  // Mass SB
  h_WHAK8J1MassSB = new TH1D("WHAK8J1MassSB","AK8J1 Mass",60,0,300);
  h_WHAK8J1MassNo2bTag = new TH1D("WHAK8J1MassNo2bTag","AK8J1 Mass",60,0,300);
  h_WHAK8J2MassSB = new TH1D("WHAK8J2MassSB","AK8J2 Mass",60,0,300);
  h_WHAK8J2MassNo2bTag = new TH1D("WHAK8J2MassNo2bTag","AK8J2 Mass",60,0,300);  
  h_HWAK8J1MassSB = new TH1D("HWAK8J1MassSB","AK8J1 Mass",60,0,300);
  h_HWAK8J1MassNo2bTag = new TH1D("HWAK8J1MassNo2bTag","AK8J1 Mass",60,0,300);
  h_HWAK8J2MassSB = new TH1D("HWAK8J2MassSB","AK8J2 Mass",60,0,300);
  h_HWAK8J2MassNo2bTag = new TH1D("HWAK8J2MassNo2bTag","AK8J1 Mass",60,0,300);

  // gen Reco match
  h_AK8J1doubleBDis = new TH1D("AK8J1doubleBDis","AK8 J1 doubleB disc.",100,0,1);
  h_AK8J2doubleBDis = new TH1D("AK8J2doubleBDis","AK8 J2 doubleB disc.",100,0,1);
  h_AK8J1deepdoubleBDis = new TH1D("AK8J1deepdoubleBDisH","AK8 J1 deepdoubleB disc.",100,0,1);
  h_AK8J2deepdoubleBDis = new TH1D("AK8J2deepdoubleBDisH","AK8 J2 deepdoubleB disc.",100,0,1);
  h_AK8J1Tau21 = new TH1D("AK8J1Tau21","AK8 J1 Tau21",100,0,1);
  h_AK8J2Tau21 = new TH1D("AK8J2Tau21","AK8 J2 Tau21",100,0,1);
  h_AK8J1zDis = new TH1D("AK8J1zDis","AK8 J1 zDisc.",100,0,1);
  h_AK8J2zDis = new TH1D("AK8J2zDis","AK8 J2 zDisc.",100,0,1);
  h_AK8J1deepdoubleBDisQ = new TH1D("AK8J1deepdoubleBDisQ","AK8 J1 deepdoubleB disc.",100,0,1);
  h_AK8J2deepdoubleBDisQ = new TH1D("AK8J2deepdoubleBDisQ","AK8 J2 deepdoubleB disc.",100,0,1);
  h_AK8J1zhDisMD = new TH1D("AK8J1zhDisMD","AK8 J1 zhDisc.",100,0,1);
  h_AK8J2zhDisMD = new TH1D("AK8J2zhDisMD","AK8 J2 zhDisc.",100,0,1);

  //

  //  h_mtbmin = new TH1D("mtbmin","mtbmin ",100,0,1000);
  // h_mct =  new TH1D("mcT","mcT ",100,0,1000);
  h_WHmtbmin = new TH1D("WHmtbmin","mtbmin ",100,0,1000);
  h_WHmct =  new TH1D("WHmcT","mcT ",100,0,1000);
  h2_WHmtbminMct = new TH2D("WHmtbminMct","x:mtbmin vs y:mct",100,0,1000,100,0,1000);
  h2_WHmtbminHmass = new TH2D("WHmtbminHmass","x:mtbmin vs y:mass",100,0,1000,60,0,300);

  h_HWmtbmin = new TH1D("HWmtbmin","mtbmin ",100,0,1000);
  h_HWmct =  new TH1D("HWmcT","mcT ",100,0,1000);
  h2_HWmtbminMct = new TH2D("HWmtbminMct","x:mtbmin vs y:mct",100,0,1000,100,0,1000);
  h2_HWmtbminHmass = new TH2D("HWmtbminHmass","x:mtbmin vs y:mass",100,0,1000,60,0,300);


  // ----------
  /* h2_dRAK8J1J2 = new TH2D("dRAK8J1J2","x:deltaR(GenW,AK8J2) vs y:deltaR(GenW,AK8J2)",50,0,5,50,0,5); */
  /* h2_AK8J1J2Tau21 = new TH2D("AK8J1J2Tau21","x:AK8J1 #tau21 vs y:AK8J2 #tau21",100,0,2,100,0,2); */
  /* h2_DisdRAK8J1 = new TH2D("disdRAK8J1","x:deltaR(GenW,AK8J1) vs y:deepDeCol w discr.",50,0,2,50,0,1); */
  /* h2_DisdRAK8J2 = new TH2D("disdRAK8J2","x:deltaR(GenW,AK8J2) vs y:deepDeCol w discr.",50,0,2,50,0,1);   */
  /* h2_Tau21dRAK8J1 = new TH2D("Tau21dRAK8J1","x:deltaR(GenW,AK8J1) vs y:AK8J1 #tau21",50,0,2,50,0,1); */
  /* h2_Tau21dRAK8J2 = new TH2D("Tau21dRAK8J2","x:deltaR(GenW,AK8J2) vs y:AK8J2 #tau21",50,0,2,50,0,1); */
  
  /* h2_AK8J1Mass_J1Tau21 = new TH2D("h2_AK8J1Mass_J1Tau21","x:SoftdropmassJ1 vs y:AK8J1 #tau21",50,0,200,100,0,1); */
  /* h2_AK8J2Mass_J2Tau21 = new TH2D("h2_AK8J2Mass_J2Tau21","x:SoftdropmassJ2 vs y:AK8J2 #tau21",50,0,200,100,0,1); */

  /* h_dPhi1 = new TH1D("DeltaPhi1","DeltaPhi1",40,0,4); */
  /* h_dPhi2 = new TH1D("DeltaPhi2","DeltaPhi2",40,0,4); */
  /* h_dPhi3 = new TH1D("DeltaPhi3","DeltaPhi3",40,0,4); */
  /* h_dPhi4 = new TH1D("DeltaPhi4","DeltaPhi4",40,0,4); */

}

SignalReg::SignalReg(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TString nameData2 = nameData;
  TChain *tree = new TChain("tree");
  //  if(nameData2.Contains("TChiWZ")) tree = new TChain("TreeMaker2/PreSelection");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t SignalReg::FillChain(TChain *chain, const TString &inputFileList) {
  int itr=0;
  TFile *filePointer;
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //    std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t SignalReg::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

SignalReg::~SignalReg() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

