#define SignalReg_cxx
#include "SignalReg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) 

using namespace std;
//

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

  TString s_data = data;
  TString s_runlist = inputFileList;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  bool isFastSim = false;
  // float xsec = 0.0, numEvents = 0.0;
  Long64_t nEvtSurv = 0;
  int ak8J1Idx = -1;

  h_cutflow->Fill("0",0);
  h_cutflow->Fill("Weighted",0);    
  h_cutflow->Fill("MET>200",0);
  h_cutflow->Fill("dPhiCuts",0);
  h_cutflow->Fill("photonVeto",0);
  h_cutflow->Fill("LVeto",0);
  h_cutflow->Fill("Filters",0);
  h_cutflow->Fill("HTRatioDphi",0);
  h_cutflow->Fill("PassedTrigger",0);
  h_cutflow->Fill("HEMaffected",0);
  h_cutflow->Fill("NEvtsNoWtLeft",0);

  h_filters->Fill("TotEvnts",0);
  h_filters->Fill("globalSuperTightHalo2016Filter",0);
  h_filters->Fill("HBHENoiseFilter",0);
  h_filters->Fill("HBHEIsoNoiseFilter",0);
  h_filters->Fill("eeBadScFilter",0);
  h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",0);
  h_filters->Fill("BadChargedCandidateFilter",0);
  h_filters->Fill("BadPFMuonFilter",0);
  h_filters->Fill("NVtx>0",0);
  h_filters->Fill("JetID",0);
  h_filters->Fill("(MET/CaloMET<5.)",0);

  int dataRun = 0;
  if(s_data.Contains("MC_2016")){ dataRun = -2016; lumiInfb = 35.815165; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("MC_2017")){ dataRun = -2017; lumiInfb = 41.486136; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("MC_2018")){ dataRun = -2018; lumiInfb = 59.546381; deepCSVvalue = 0.4184;}
  
  else if(s_data.Contains("2016")){ dataRun = 2016; isMC = false; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("2017")){ dataRun = 2017; isMC = false; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("2018")){ dataRun = 2018; isMC = false; deepCSVvalue = 0.4184;}

  if(s_data.Contains("TChi")){
    isFastSim = true;
    lumiInfb = 137.0;
    deepCSVvalue = 0.4184;
    deepAK8Wscore = 0.918;
    //   if(s_data.Contains("TChiWZ_1000")){ xsec = 1.34352e-3; numEvents = 28771;}
    //   else if(s_data.Contains("TChiWZ_800")){ xsec = 4.75843e-3; numEvents = 34036;}
    //   cout<<"Assigning xsec as: "<<xsec<<endl;
  }
  
  bool TTJets_nonHTcut = false;
  if (s_runlist.Contains("TTJets_DiLept") || s_runlist.Contains("TTJets_SingleLeptFromT") ){
    cout <<" *****  Applying madHT < 600 cut to add other HT samples > 600"<< endl;
    TTJets_nonHTcut = true;
  }

  lumiInfb = 137.0;
  cout<<"!!!! changing intLumi to 137/fb, although you should have used 2018 intLumi...."<<endl;

  if(dataRun>0) cout<<"Processing it as "<<dataRun<<" data"<<endl;
  else if(dataRun<0) cout<<"Processing it as "<<abs(dataRun)<<" MC"<<endl;
  else cout<<"No specific data/MC year"<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    print(jentry);  
    // if(dataRun==-2016 || dataRun==-2017) wt=Weight*1000.0*lumiInfb*NonPrefiringProb;
    // else if(dataRun <=0) wt=Weight*1000.0*lumiInfb;
    //    else wt = 1.0;
    // if(isFastSim && NumEvents==1 && CrossSection==1){
    //   CrossSection = xsec;
    //   NumEvents = numEvents;
    //   Weight = xsec/numEvents;
    // }
    wt=Weight*1000.0*lumiInfb;

    h_cutflow->Fill("0",1);
    h_cutflow->Fill("Weighted",wt);
    //--------------
    //    if(jentry>100) break;

    //#################### EWK Planned baseline cuts
    if(isFastSim) JetID = true;

    //   Calc. mtbmin, mt and mt2j
    Double_t mtbmin = -99999. , mct = -99999. ; //mCT is the contransverse mass variable 
    double mt = 0, mt2j = 0;
    //    mtbmin = find_bjets_mtbmin();           // loops over b-jets and finds mtmin

    // if(JetsAK8->size() > 0) mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
    if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2 ) mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
    if(JetsAK8->size()>=2) mt2j = sqrt(2*(*JetsAK8)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[1].Phi()))));
    
    if (TTJets_nonHTcut){
      if (madHT > 600) continue;  //madHT <= 600 for non HT
    }
    
    //for madHT stiching check in TTJets
    h_HT->Fill(HT,wt);  
    h_madHT->Fill(madHT,wt);      
    
    
    // Modified RA2b cuts
    if(NJets < 2                                                                //AK4 jets >=2
       || HT < 300 || MHT < 200 || MET < 250                                     // HT, MET, MHT cuts
       || NMuons!=0 || NElectrons!=0                                             // Veto e,muons
       || (MHT/HT > 1.0) || !JetID                                               // MHT<HT and Veto JET ID
       || !(DeltaPhi1 > 1.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3)  //Angle cuts
       || isoMuonTracks!=0 || isoElectronTracks!=0 || isoPionTracks!=0) continue; //Veto isolated tracks

    //Boosted AK8 jet cuts
    
    if (JetsAK8->size() < 2                                                     // require >=2 AK8 jets
	// || (*JetsAK8)[0].Pt() < 200 || (*JetsAK8)[1].Pt() < 200         // AK8 jets pT >200
	// || abs((*JetsAK8)[0].Eta()) > 2 || abs((*JetsAK8)[1].Eta()) > 2       // jets |Eta| < 2   
	) continue; 
    
    h_NJets->Fill(NJets,wt);
    h_MT->Fill(mt,wt);
    
    //+ some additional cuts
    if(NJets > 6                                                       //AK4 jets <= 6
       || mt < 500                                                      // mt > 500
       //||mtbmin < 200                                                   // mTbmin > 200 to reduce ttbar bkg    
       ) continue; 
          
    // ******* Choose  b veto
    if( BTagsDeepCSV!=0 ) continue;                                              //Veto B jets   old cut BTags!=0 

    //    if(MET < 200) continue;
    //h_cutflow->Fill("MET>200",wt);
    float dphi1=4, dphi2=4, dphi3=4, dphi4=4;
    //if(Jets->size() > 0 && (*Jets)[0].Pt() > 30 && abs((*Jets)[0].Eta()) < 6.0)
    //  dphi1 = (abs(DeltaPhi(METPhi,(*Jets)[0].Phi())));
    if(Jets->size() > 0){
      if ((*Jets)[0].Pt() > 30 && abs((*Jets)[0].Eta()) < 6.0)
	dphi1 = (abs(DeltaPhi(METPhi,(*Jets)[0].Phi())));
    }
    
    if(Jets->size() > 1 && (*Jets)[1].Pt() > 30 && abs((*Jets)[1].Eta()) < 6.0)
      dphi2 = (abs(DeltaPhi(METPhi,(*Jets)[1].Phi())));

    if(Jets->size() > 2 && (*Jets)[2].Pt() > 30 && abs((*Jets)[2].Eta()) < 6.0)
      dphi3 = (abs(DeltaPhi(METPhi,(*Jets)[2].Phi())));
    
    if(Jets->size() > 3 && (*Jets)[3].Pt() > 30 && abs((*Jets)[3].Eta()) < 6.0)
      dphi4 = (abs(DeltaPhi(METPhi,(*Jets)[3].Phi())));
    
    //----Photon veto
    int nPhotons=0;
    for(int i=0;i<Photons->size();i++){
      if((*Photons)[i].Pt() > 100 && (*Photons_fullID)[i] && (!(*Photons_hasPixelSeed)[i]) ){ nPhotons++; break;}
    }
    if(nPhotons>0) continue;
    //if(Photons->size()!=0) continue;
    h_cutflow->Fill("photonVeto",wt);

    h_filters->Fill("TotEvnts",1);
    h_filters->Fill("globalSuperTightHalo2016Filter",globalSuperTightHalo2016Filter);
    h_filters->Fill("HBHENoiseFilter",HBHENoiseFilter);
    h_filters->Fill("HBHEIsoNoiseFilter",HBHEIsoNoiseFilter);
    h_filters->Fill("eeBadScFilter",eeBadScFilter);
    h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",EcalDeadCellTriggerPrimitiveFilter);
    h_filters->Fill("BadChargedCandidateFilter",BadChargedCandidateFilter);
    h_filters->Fill("BadPFMuonFilter",BadPFMuonFilter);
    h_filters->Fill("NVtx>0",NVtx > 0);
    h_filters->Fill("JetID",JetID);
    h_filters->Fill("(MET/CaloMET<5.)",(MET/CaloMET < 5.));

    if(!isFastSim){
      if(!(globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0 && JetID && (MET/CaloMET < 5.))) continue;
      h_cutflow->Fill("Filters",wt);
      h_cutflow->Fill("HTRatioDphi",wt);
    }
    //--------------------------triggers
    if(!isMC){
      bool trgPass = false;
      TString trgName;
      for(int i=0;i<TriggerPass->size();i++){
	trgName = (*TriggerNames)[i];
	if(!(trgName.Contains("MET"))) continue;
	if((*TriggerPass)[i]==1 && (trgName.Contains("HLT_PFMET100_PFMHT100_IDTight_v") || trgName.Contains("HLT_PFMET110_PFMHT110_IDTight_v") ||
				    trgName.Contains("HLT_PFMET120_PFMHT120_IDTight_v") || trgName.Contains("HLT_PFMET130_PFMHT130_IDTight_v") ||
				    trgName.Contains("HLT_PFMET140_PFMHT140_IDTight_v") || trgName.Contains("HLT_PFMET90_PFMHT90_IDTight_v") || 
				    trgName.Contains("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") || trgName.Contains("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") ||
				    trgName.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") || trgName.Contains("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") || 
				    trgName.Contains("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") ||  trgName.Contains("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"))) trgPass = true;
      }
      if(trgPass) h_cutflow->Fill("PassedTrigger",wt);
      else continue;
    }
    bool HEMaffected = false;
    if(dataRun==2018 && RunNum >=319077){
      for(int i=0;i<Jets->size();i++){
	if((*Jets)[i].Pt() < 30) continue;
	if( (*Jets)[i].Eta() >= -3.20 && (*Jets)[i].Eta() <= -1.2 && 
	    (*Jets)[i].Phi() >= -1.77 && (*Jets)[i].Phi() <= -0.67 &&
	    (abs(DeltaPhi(METPhi,(*Jets)[i].Phi())) < 0.5) ){HEMaffected = true; break;}
      }
    }
    //--------------------------end of triggers

    if(HEMaffected){
      h_cutflow->Fill("HEMaffected",wt);
      continue;
    }
    
    // h_MET->Fill(MET,wt);
    // h_MHT->Fill(MHT,wt);
    // h_HT->Fill(HT,wt);
    // h_NJets->Fill(NJets,wt);
    // h_BTags->Fill(BTagsDeepCSV,wt);
    // h_dPhi1->Fill(DeltaPhi1,wt);
    // h_dPhi2->Fill(DeltaPhi2,wt);
    // h_dPhi3->Fill(DeltaPhi3,wt);
    // h_dPhi4->Fill(DeltaPhi4,wt);
    
    //cout<<"jets ak8 size "<<JetsAK8->size()<<endl;
    //cout<<"jets disc size "<<JetsAK8_zDiscriminatorDeep->size()<<endl;
    //cout<<"\n";
    
    bool Bkgest_MET = false;
    if (Bkgest_MET){
      // W H 
      if( 
	 ((*JetsAK8_softDropMass)[0] > massLowW) &&
	 ((*JetsAK8_softDropMass)[0] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[0] > deepAK8Wscore) )
	{
	  if(         //middle band
	     ((*JetsAK8_softDropMass)[1] > massLowH) &&
	     ((*JetsAK8_softDropMass)[1] < massHighH) ){
	    
	    if((*JetsAK8_doubleBDiscriminator)[1] < deepbbscore){ 
	      h_WHMET_RegC->Fill(MET,wt);
	      h_WHAK8J2Mass_RegC->Fill((*JetsAK8_softDropMass)[1],wt);
	    }
	    if( (*JetsAK8_doubleBDiscriminator)[1] > deepbbscore)  // Signal region **
	      {
		h_WHMET_RegA->Fill(MET,wt);
		h_WHAK8J2Mass_RegA->Fill((*JetsAK8_softDropMass)[1],wt);
	      }
	  }
	  
	  if( //side bands
	     (((*JetsAK8_softDropMass)[1] >50) && ((*JetsAK8_softDropMass)[1] < massLowH)) ||
	     (((*JetsAK8_softDropMass)[1] > 200) && ((*JetsAK8_softDropMass)[1] < 250)) ){

	    if( (*JetsAK8_doubleBDiscriminator)[1] < deepbbscore){
	      h_WHMET_RegD->Fill(MET,wt);
	      h_WHAK8J2Mass_RegD->Fill((*JetsAK8_softDropMass)[1],wt);
	    }
	    if( (*JetsAK8_doubleBDiscriminator)[1] > deepbbscore){
	      h_WHMET_RegB->Fill(MET,wt);
	      h_WHAK8J2Mass_RegB->Fill((*JetsAK8_softDropMass)[1],wt);
	    }
	  }
	} //end WH
    
      //     H W
      if( 
	 ((*JetsAK8_softDropMass)[1] > massLowW) &&
	 ((*JetsAK8_softDropMass)[1] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[1] > deepAK8Wscore) )
	{
	  if(         //middle band or signal band
	     ((*JetsAK8_softDropMass)[0] > massLowH) &&
	     ((*JetsAK8_softDropMass)[0] < massHighH) ){
	    
	    if( (*JetsAK8_doubleBDiscriminator)[0] < deepbbscore ){
	      h_HWMET_RegC->Fill(MET,wt);
	      h_HWAK8J1Mass_RegC->Fill((*JetsAK8_softDropMass)[0],wt);
	    }
	    if( (*JetsAK8_doubleBDiscriminator)[0] > deepbbscore ){
	      h_HWMET_RegA->Fill(MET,wt);
	      h_HWAK8J1Mass_RegA->Fill((*JetsAK8_softDropMass)[0],wt);
	    }
	  }
	  
	  if( //side bands
	     (((*JetsAK8_softDropMass)[0] > 50) && ((*JetsAK8_softDropMass)[0] < massLowH)) ||
	     (((*JetsAK8_softDropMass)[0] > 200) && ((*JetsAK8_softDropMass)[0] < 250)) ){
	    
	    if( (*JetsAK8_doubleBDiscriminator)[0] < deepbbscore ){
	      h_HWMET_RegD->Fill(MET,wt);
	      h_HWAK8J1Mass_RegD->Fill((*JetsAK8_softDropMass)[0],wt);
	    }	      
	    if( (*JetsAK8_doubleBDiscriminator)[0] > deepbbscore ){
	      h_HWMET_RegB->Fill(MET,wt);
	      h_HWAK8J1Mass_RegB->Fill((*JetsAK8_softDropMass)[0],wt);
	    }
	  }
	} //end HW
	
    } //end if 
    
    bool Mass_SB = false;
    if(Mass_SB){
      if( // W H
	 ((*JetsAK8_softDropMass)[0] > massLowW) &&
	 ((*JetsAK8_softDropMass)[0] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[0] > deepAK8Wscore) &&
	 ((*JetsAK8_softDropMass)[1] > 50) &&
	 ((*JetsAK8_softDropMass)[1] < 250) ){
	  
	//if( ((*JetsAK8_doubleBDiscriminator)[1] > bbscore) ) {
	if( ((*JetsAK8_deepDoubleBDiscriminatorH)[1] > deepbbscore) ){ // Sig reg mass range
	  h_WHAK8J1MassSB->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_WHAK8J2MassSB->Fill((*JetsAK8_softDropMass)[1],wt); 
	}
	  
	//if( ((*JetsAK8_doubleBDiscriminator)[1] < bbscore) ) {
	if( ((*JetsAK8_deepDoubleBDiscriminatorH)[1] < deepbbscore) ){ // 2b side band 
	  h_WHAK8J1MassNo2bTag->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_WHAK8J2MassNo2bTag->Fill((*JetsAK8_softDropMass)[1],wt); 
	}
      }
		
      if( // H W
	 ((*JetsAK8_softDropMass)[1] > massLowW) &&
	 ((*JetsAK8_softDropMass)[1] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[1] > deepAK8Wscore) &&
	 ((*JetsAK8_softDropMass)[0] > 50) &&
	 ((*JetsAK8_softDropMass)[0] < 250) ){
	  
	//if( ((*JetsAK8_doubleBDiscriminator)[0] > bbscore) ) {
	if( ((*JetsAK8_deepDoubleBDiscriminatorH)[0] > deepbbscore) ){
	  h_HWAK8J1MassSB->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_HWAK8J2MassSB->Fill((*JetsAK8_softDropMass)[1],wt); 
	}
	//if( ((*JetsAK8_doubleBDiscriminator)[0] < bbscore) ) {
	if( ((*JetsAK8_deepDoubleBDiscriminatorH)[0] < deepbbscore) ){
	  h_HWAK8J1MassNo2bTag->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_HWAK8J2MassNo2bTag->Fill((*JetsAK8_softDropMass)[1],wt); 
	}
      }
     
    } //end Mass SB
      
    bool WH_SR = false;
    if(WH_SR){
      //-----------    W H Signal region
      if(
	 ((*JetsAK8_softDropMass)[0] > massLowW) &&
	 ((*JetsAK8_softDropMass)[0] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[0] > deepAK8Wscore) &&
	 ((*JetsAK8_softDropMass)[1] > massLowH) &&
	 ((*JetsAK8_softDropMass)[1] < massHighH) &&
	 //	   ((*JetsAK8_doubleBDiscriminator)[1] >bbscore) ) {
	 ((*JetsAK8_deepDoubleBDiscriminatorH)[1] > deepbbscore) ) {      
	  
	h_WHMET->Fill(MET,wt);
	h_WHMT->Fill(mt,wt);
	h_WHMT2J->Fill(mt2j,wt);
      
	h_WHAK8J1Pt->Fill(((*JetsAK8)[0].Pt()),wt);
	h_WHAK8J1Eta->Fill(((*JetsAK8)[0].Eta()),wt);
	h_WHAK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt);
	h_WHAK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[0])/((*JetsAK8_NsubjettinessTau1)[0]),wt);
	h_WHAK8J1wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[0]),wt);
      
	h_WHAK8J2Pt->Fill(((*JetsAK8)[1].Pt()),wt);
	h_WHAK8J2Eta->Fill(((*JetsAK8)[1].Eta()),wt);
	h_WHAK8J2Mass->Fill((*JetsAK8_softDropMass)[1],wt);
	h_WHAK8J2Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[1])/((*JetsAK8_NsubjettinessTau1)[1]),wt);
	h_WHAK8J2wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[1]),wt);
	h_WHmtbmin->Fill(mtbmin,wt);
	h_WHmct->Fill(mct,wt);
	h2_WHmtbminMct->Fill(mtbmin,mct,wt);
	h2_WHmtbminHmass->Fill(mtbmin,(*JetsAK8_softDropMass)[1],wt);
      }
      //------------    H W  Signal region    
      if( 
	 ((*JetsAK8_softDropMass)[1] > massLowW) &&
	 ((*JetsAK8_softDropMass)[1] < massHighW) &&
	 ((*JetsAK8_wDiscriminatorDeep)[1] > deepAK8Wscore) &&
	 ((*JetsAK8_softDropMass)[0] > massLowH) &&
	 ((*JetsAK8_softDropMass)[0] < massHighH) &&
	 //	   ((*JetsAK8_doubleBDiscriminator)[0] > bbscore) ) {
	 ((*JetsAK8_deepDoubleBDiscriminatorH)[0] > deepbbscore) ){
	  
	h_HWMET->Fill(MET,wt);
	h_HWMT->Fill(mt,wt);
	h_HWMT2J->Fill(mt2j,wt);

	h_HWAK8J1Pt->Fill(((*JetsAK8)[0].Pt()),wt);
	h_HWAK8J1Eta->Fill(((*JetsAK8)[0].Eta()),wt);
	h_HWAK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt);
	h_HWAK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[0])/((*JetsAK8_NsubjettinessTau1)[0]),wt);
	h_HWAK8J1wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[0]),wt);

	h_HWAK8J2Pt->Fill(((*JetsAK8)[1].Pt()),wt);
	h_HWAK8J2Eta->Fill(((*JetsAK8)[1].Eta()),wt);
	h_HWAK8J2Mass->Fill((*JetsAK8_softDropMass)[1],wt);
	h_HWAK8J2Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[1])/((*JetsAK8_NsubjettinessTau1)[1]),wt);
	h_HWAK8J2wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[1]),wt);
	h_HWmtbmin->Fill(mtbmin,wt);
	h_HWmct->Fill(mct,wt);
	h2_HWmtbminMct->Fill(mtbmin,mct,wt);
	h2_HWmtbminHmass->Fill(mtbmin,(*JetsAK8_softDropMass)[0],wt);

      } //end H W
     
    } // end WH Signal region  
     
    bool new_WZWH_SR = false; // b is vetoed so this doesnt work
    if(new_WZWH_SR){
      vector<bool> jetsAK8hasb;
      for(int i=0;i<JetsAK8->size();i++){
	bool foundbInAK8 = 0;
	for(int b=0;b<bjets.size();b++){
	  if(bjets[b].DeltaR((*JetsAK8)[i]) < 0.8){
	    foundbInAK8 = 1;
	    break;
	  }
	}
	jetsAK8hasb.push_back(foundbInAK8);    
      }// AK8 loop for b content

      bool ZorHcand = false;
      bool Wcand = false;
      for(int i=0;i<JetsAK8->size();i++){
	if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
	  if(jetsAK8hasb[i]){ //FOR Z OR H CAND.
	
	    if(
	       ((*JetsAK8_softDropMass)[i] > massLowZ) &&
	       ((*JetsAK8_softDropMass)[i] < massHighZ) &&
	       ((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepbbscore) ) { 
	      ZorHcand = true;
	    }
	  }
	  else{ //W cand
	    
	    if(
	       ((*JetsAK8_softDropMass)[i] > massLowW) &&                                   	   
	       ((*JetsAK8_softDropMass)[i] < massHighW) &&                                                
	       ((*JetsAK8_wDiscriminatorDeep)[i] > deepAK8Wscore) ){
	      Wcand = true;
	    }
	  }
	}
      }
      if (ZorHcand && Wcand){
	h_whMET->Fill(MET,wt);
	h_whMT->Fill(mt,wt);
	h_whMT2J->Fill(mt2j,wt);
      }
    }

    bool WW_SR = false;
    if(WW_SR){
      //-----------    W W Signal region
      if(
	 ((*JetsAK8_softDropMass)[0] > massLowW) &&
	 ((*JetsAK8_softDropMass)[0] < massHighW) &&
	 //((*JetsAK8_NsubjettinessTau2)[0]/(*JetsAK8_NsubjettinessTau1)[0]) < 0.35 &&
	 ((*JetsAK8_wDiscriminatorDeep)[0] > deepAK8Wscore) &&
	 ((*JetsAK8_softDropMass)[1] > massLowW) &&
	 ((*JetsAK8_softDropMass)[1] < massHighW) &&
	 //((*JetsAK8_NsubjettinessTau2)[1]/(*JetsAK8_NsubjettinessTau1)[1]) < 0.35  ) { 
	 ((*JetsAK8_wDiscriminatorDeep)[1] > deepAK8Wscore) ) { 
	
	h_WWMET->Fill(MET,wt);
	h_WWMT->Fill(mt,wt);
	h_WWMT2J->Fill(mt2j,wt);
      
	h_WWAK8J1Pt->Fill(((*JetsAK8)[0].Pt()),wt);
	h_WWAK8J1Eta->Fill(((*JetsAK8)[0].Eta()),wt);
	h_WWAK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt);
      }
    } //end WW SR
    
    bool WZW_SR = true;
    if(WZW_SR){
      vector<double>  Wdis, WdisMD;
      double Wdismax,WdisMDmax = 0;
      Wdis.resize(0);
      WdisMD.resize(0);
      //-----------    WW/ WZ Signal region
      for(int w=0; w < JetsAK8->size(); w++){
	if(
	   ((*JetsAK8_softDropMass)[w] > 65) &&
	   ((*JetsAK8_softDropMass)[w] < 105) &&
	   ((*JetsAK8)[w].Pt()) > 200 &&                                                                      
	   (abs((*JetsAK8)[w].Eta()) < 2) ){ // W jets mass,eta,pt

	  if( (*JetsAK8_wDiscriminatorDeep)[w] > 0.918 ){
	    Wdis.push_back((*JetsAK8_wDiscriminatorDeep)[w]);
	    Wdismax = *max_element(Wdis.begin(), Wdis.end());  //WDis with highest value passing the cut
	  }
	  if((*JetsAK8_wDiscriminatorDeepDecorrel)[w] > 0.704 ){
	    WdisMD.push_back((*JetsAK8_wDiscriminatorDeepDecorrel)[w]);
	    WdisMDmax = *max_element(WdisMD.begin(), WdisMD.end()); 
	  } 
	}
      }	
      
      
      // if(Wdismax >0)  cout<<"W dis max "<<   Wdismax <<endl;	
      // if(Wdismax >0)  cout<<"W dis MD max "<<   WdisMDmax <<endl;	
      
      bool Wcand = false;
      bool worZcand = false;
      bool disandtau21 = false;
      bool nothighWdisc = false;
      bool exception = false;
      int widx = 99;
      for(int i=0; i < JetsAK8->size(); i++){ // W
      	if ( Wdismax > 0 &&
	     (*JetsAK8_wDiscriminatorDeep)[i] == Wdismax &&  // For max W-tag
	     //(*JetsAK8_wDiscriminatorDeepDecorrel)[i] == WdisMDmax &&
	     ((*JetsAK8_softDropMass)[i] > 65) &&     // ?? already passes by above cut
	     ((*JetsAK8_softDropMass)[i] < 105) &&
	     ((*JetsAK8)[i].Pt()) > 200 &&                                                                      
	     (abs((*JetsAK8)[i].Eta()) < 2) ){ // highest score W cand.
	  //      	  cout<<" i "<< i<<endl;
	  widx = i;
	  Wcand = true;
	  //cout<<"Wcand true"<<endl;
      	}
      }
      // cout<<"w idx  " <<widx<<endl;
      
      for(int j=0; j < JetsAK8->size(); j++){  //W/Z
	//if( (*JetsAK8_wDiscriminatorDeep)[j] > 0.918 || (*JetsAK8_zDiscriminatorDeep)[j] > 0.918 ){
	if( (j != widx) &&
	    // (*JetsAK8_wDiscriminatorDeepDecorrel)[j] > 0.704 &&
	    (((*JetsAK8_NsubjettinessTau2)[j]/(*JetsAK8_NsubjettinessTau1)[j]) < 0.35) &&
	    ((*JetsAK8_softDropMass)[j] > 65) &&
	    ((*JetsAK8_softDropMass)[j] < 105) &&
	    ((*JetsAK8)[j].Pt()) > 200 &&                                                                      
	    (abs((*JetsAK8)[j].Eta()) < 2) ){ // score W/Z cand.
	  //  cout<<" j "<< j<<endl;
	  worZcand=true;
	  //cout<<"worZcand true"<<endl;
	} // fails if the highest Wdisc. cand also passes tau21
      }

      for(int k=0; k < JetsAK8->size(); k++){  //W/Z
	if( (k == widx) &&
	    // (*JetsAK8_wDiscriminatorDeepDecorrel)[j] > 0.704 &&
	    (((*JetsAK8_NsubjettinessTau2)[k]/(*JetsAK8_NsubjettinessTau1)[k]) < 0.35) &&
	    ((*JetsAK8_softDropMass)[k] > 65) &&
	    ((*JetsAK8_softDropMass)[k] < 105) &&
	    ((*JetsAK8)[k].Pt()) > 200 &&                                                                      
	    (abs((*JetsAK8)[k].Eta()) < 2) ){ // score W/Z cand.
	  disandtau21 = true;
	}
      	if( (k != widx) &&
	    // (*JetsAK8_wDiscriminatorDeepDecorrel)[j] > 0.704 &&
	    (*JetsAK8_wDiscriminatorDeep)[k] > 0.918 &&
	    ((*JetsAK8_softDropMass)[k] > 65) &&
	    ((*JetsAK8_softDropMass)[k] < 105) &&
	    ((*JetsAK8)[k].Pt()) > 200 &&                                                                      
	    (abs((*JetsAK8)[k].Eta()) < 2) ){ // score W/Z cand.
	  nothighWdisc = true;
	}
	if (disandtau21 && nothighWdisc){
	  exception = true;
	}
      }
	
      if(Wcand && (worZcand || exception) ){
	//	cout<<"both true"<<endl;
	h_WZW1MET->Fill(MET,wt);
	h_METvBin->Fill(MET,wt);
	h_WZW1MT->Fill(mt,wt);
	h_WZW1MT2J->Fill(mt2j,wt);
      }
      //      cout<<"******"<<endl;
      /*  
	  if(
	 ((*JetsAK8_softDropMass)[0] > 65) &&
	 ((*JetsAK8_softDropMass)[0] < 90) &&
	 ((*JetsAK8_softDropMass)[1] > 65) &&
	 ((*JetsAK8_softDropMass)[1] < 90) ){
	
	if(
	   ( ((*JetsAK8_NsubjettinessTau2)[0]/(*JetsAK8_NsubjettinessTau1)[0]) < 0.35) &&	
	   ( ((*JetsAK8_NsubjettinessTau2)[1]/(*JetsAK8_NsubjettinessTau1)[1]) < 0.35) ){	
	  cout<<"tau21"<<endl;
	  h_WZW1MET->Fill(MET,wt);
	  h_WZW1MT->Fill(mt,wt);
	  h_WZW1MT2J->Fill(mt2j,wt);
      	}
	double W_disc1,W_disc2 =0;
	if(
	   ((*JetsAK8_wDiscriminatorDeep)[0] > 0.918) &&  // 1% mistag 
	   ((*JetsAK8_wDiscriminatorDeep)[1] > 0.918) ) { 
	  W_disc1 = MAX((*JetsAK8_wDiscriminatorDeep)[0],(*JetsAK8_wDiscriminatorDeep)[1]);
	  cout<<"Wdis1"<<endl;
	  h_WZW2MET->Fill(MET,wt);
	  h_WZW2MT->Fill(mt,wt);
	  h_WZW2MT2J->Fill(mt2j,wt);
      	}
	if(
	   ((*JetsAK8_wDiscriminatorDeep)[0] > 0.961) &&   // 0.5% mistag 
	   ((*JetsAK8_wDiscriminatorDeep)[1] > 0.961) ) { 
	  W_disc2 = MAX((*JetsAK8_wDiscriminatorDeep)[0],(*JetsAK8_wDiscriminatorDeep)[1]);
	   cout<<"Wdis2"<<endl;
	  h_WZW3MET->Fill(MET,wt);
	  h_WZW3MT->Fill(mt,wt);
	  h_WZW3MT2J->Fill(mt2j,wt);
      	}
	
	// highest of above 2 for one jet
	double high_W_disc = 0;
	high_W_disc = MAX(W_disc1, W_disc2);
	if(
	   ((*JetsAK8_wDiscriminatorDeep)[1] ==  high_W_disc) &&   // highest W in above 2 cases 
	   ((*JetsAK8_NsubjettinessTau2)[0]/(*JetsAK8_NsubjettinessTau1)[0]) < 0.35 ){
	  cout<<"W and tau21"<<endl;
	  h_WZW4MET->Fill(MET,wt);
	  h_WZW4MT->Fill(mt,wt);
	  h_WZW4MT2J->Fill(mt2j,wt);
      	}
      */
    } //end WZWWW SR
  
    

    bool WorZ_genrecomatch = false; // not b
    if (WorZ_genrecomatch){
    
      if( 
	 ((*JetsAK8_softDropMass)[1] > 65 && (*JetsAK8_softDropMass)[1] < 90) && 
	 ((*JetsAK8_softDropMass)[0] > 75 && (*JetsAK8_softDropMass)[0] < 100) ){ 
	
    	for (int gen=0; gen < GenParticles_PdgId->size(); gen++){
	  
    	  if( abs((*GenParticles_PdgId)[gen] == 2) || abs((*GenParticles_PdgId)[gen] == 1) || abs((*GenParticles_PdgId)[gen] == 3) || abs((*GenParticles_PdgId)[gen] == 4) ){
	    
	    if( abs((*GenParticles_ParentId)[gen]) == 24 ){ //W to u d s c
	      int GenWIdx = (*GenParticles_ParentIdx)[gen]; 
	      if(abs(DeltaR((*GenParticles)[GenWIdx].Eta(),(*GenParticles)[GenWIdx].Phi(),(*JetsAK8)[0].Eta(),(*JetsAK8)[0].Phi())) < 0.1){
		//      cout<<" ************  In gen reco match"<<endl;
	
		h_wAK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[1])/((*JetsAK8_NsubjettinessTau1)[1]),wt);
		h_wAK8J1zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[1]),wt);
		h_wAK8J1wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[1]),wt);
		h_wAK8J1wDisDC->Fill(((*JetsAK8_wDiscriminatorDeepDecorrel)[1]),wt);
	      }
	    }
	  
	    if( abs((*GenParticles_ParentId)[gen]) == 23 ){ //Z to u d s c
	      int GenZIdx = (*GenParticles_ParentIdx)[gen]; 
	      if(abs(DeltaR((*GenParticles)[GenZIdx].Eta(),(*GenParticles)[GenZIdx].Phi(),(*JetsAK8)[0].Eta(),(*JetsAK8)[0].Phi())) < 0.1){
		//cout<<" ************  In gen reco match  Z"<<endl;
		h_zAK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[0])/((*JetsAK8_NsubjettinessTau1)[0]),wt);
		h_zAK8J1zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[0]),wt);
		h_zAK8J1wDis->Fill(((*JetsAK8_wDiscriminatorDeep)[0]),wt);
		h_zAK8J1wDisDC->Fill(((*JetsAK8_wDiscriminatorDeepDecorrel)[0]),wt);
	      }
	    }
	  } //ucsd
	} //gen loop end
      }
    }
  
            
  

    // -----------------------     tweaking wdisc score
    bool wDisTweak = false;
    if(wDisTweak){
      if( (massLowW< (*JetsAK8_softDropMass)[0] && (*JetsAK8_softDropMass)[0] < massHighW) &&
	  (massHighW< (*JetsAK8_softDropMass)[1] && (*JetsAK8_softDropMass)[1] < massHighH) &&
	  (*JetsAK8_doubleBDiscriminator)[1] >bbscore  ){
	if( (*JetsAK8_wDiscriminatorDeep)[0] > 0.779) {
	  h_WHMETa->Fill(MET,wt);
	}
	if( (*JetsAK8_wDiscriminatorDeep)[0] > 0.981) {
	  h_WHMETc->Fill(MET,wt);
	}
      }
      if( (massLowW< (*JetsAK8_softDropMass)[1] && (*JetsAK8_softDropMass)[1] < massHighW) &&
	  (massHighW< (*JetsAK8_softDropMass)[0] && (*JetsAK8_softDropMass)[0] < massHighH) &&
	  (*JetsAK8_doubleBDiscriminator)[0] > bbscore)  {
	if( (*JetsAK8_wDiscriminatorDeep)[1] > 0.779) {
	  h_HWMETa->Fill(MET,wt);
	}
	if( (*JetsAK8_wDiscriminatorDeep)[1] > 0.981) {
	  h_HWMETc->Fill(MET,wt);
	}
      }
    } // end wDis tweak
    
    //-----------------     Gen Reco match for ROC Curves
    bool ZorHtobb_genrecomatch = false;
    if (ZorHtobb_genrecomatch){
      if(
	 ((*JetsAK8_softDropMass)[1] > massLowW && (*JetsAK8_softDropMass)[1] < massHighW) && 
	 ((*JetsAK8_softDropMass)[0] > massLowH && (*JetsAK8_softDropMass)[0] < massHighH) ){ 
	
 	for (int gen=0; gen < GenParticles_PdgId->size(); gen++){
	  
	  if((*GenParticles_PdgId)[gen] == 5 && abs((*GenParticles_ParentId)[gen]) == 25){ //Z->bb //removing double count by having only 1 b and its parrent Z
	    int GenZIdx = (*GenParticles_ParentIdx)[gen]; 
	    if(abs(DeltaR((*GenParticles)[GenZIdx].Eta(),(*GenParticles)[GenZIdx].Phi(),(*JetsAK8)[0].Eta(),(*JetsAK8)[0].Phi())) < 0.1){ //1st ak8 Zbb match
	      // cout<<" ************  In gen reco match"<<endl;
	      h_AK8J1doubleBDis->Fill(((*JetsAK8_doubleBDiscriminator)[0]),wt);
	      h_AK8J1deepdoubleBDis->Fill(((*JetsAK8_deepDoubleBDiscriminatorH)[0]),wt);
	      h_AK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[0])/((*JetsAK8_NsubjettinessTau1)[0]),wt);
	      h_AK8J1zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[0]),wt);
	      h_AK8J1deepdoubleBDisQ->Fill(((*JetsAK8_deepDoubleBDiscriminatorQ)[0]),wt);
	      h_AK8J1zhDisMD->Fill(((*JetsAK8_zhDiscriminatorDeepDecorrel)[0]),wt);
	    }
	    if(abs(DeltaR((*GenParticles)[GenZIdx].Eta(),(*GenParticles)[GenZIdx].Phi(),(*JetsAK8)[1].Eta(),(*JetsAK8)[1].Phi())) < 0.1){ //2nd ak8 Zbb match
	      h_AK8J2doubleBDis->Fill(((*JetsAK8_doubleBDiscriminator)[1]),wt);
	      h_AK8J2deepdoubleBDis->Fill(((*JetsAK8_deepDoubleBDiscriminatorH)[1]),wt);  
	      h_AK8J2Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[1])/((*JetsAK8_NsubjettinessTau1)[1]),wt);
	      h_AK8J2zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[1]),wt);
	      h_AK8J2deepdoubleBDisQ->Fill(((*JetsAK8_deepDoubleBDiscriminatorQ)[1]),wt);
	      h_AK8J2zhDisMD->Fill(((*JetsAK8_zhDiscriminatorDeepDecorrel)[1]),wt);
	      
	    }
	  } // end 2nd Zbb match
	} // end for
      } // end tune mass window
    }// end if
     //print(jentry);
  
    
    nEvtSurv++;
    h_cutflow->Fill("NEvtsNoWtLeft",1);
    //    cout<<" ____ End evt loop ----------   "<<endl;
  } // loop over entries
  cout<<"No. of entries survived: "<<nEvtSurv<<endl;
}

void SignalReg::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"MomMass:"<<SusyMotherMass<<" Kid Mass:"<<SusyLSPMass<<endl;
  for(int i=0;i<GenParticles->size();i++){  
    //    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx :"<<(*GenParticles)[i].Px()<<" Py :"<<(*GenParticles)[i].Py()<<" Pz :"<<(*GenParticles)[i].Pz()<<" E: "<<(*GenParticles)[i].Energy()<<" M:"<<(*GenParticles)[i].M()<<endl;
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" M:"<<(*GenParticles)[i].M()<<endl;
  }
  cout<<endl;
  for(int i=0;i<GenJets->size();i++){
    cout<< "Size of GenJets: "<< GenJets->size() << endl;
    cout<<"GenJetPt:"<<(*GenJets)[i].Pt()<<" GenJetEta:"<<(*GenJets)[i].Eta()<<" GenJetPhi:"<<(*GenJets)[i].Phi()<<endl;
  }
  cout<<endl;
  for(int i=0;i<Jets->size();i++){
    cout<< "Size of Jets: "<< Jets->size() << endl;
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  cout<<endl;
  cout<<"MET:"<<MET<<" METPhi:"<<METPhi<<" MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
  cout<<endl;
  for(int i=0;i<JetsAK8->size();i++){
    cout<< "Size of AK8 Jets: "<< JetsAK8->size() << endl;
    cout<<"AK8 pT:"<<(*JetsAK8)[i].Pt()<<" eta:"<<(*JetsAK8)[i].Eta()<<" phi:"<<(*JetsAK8)[i].Phi()<<" softDropM:"<<(*JetsAK8_softDropMass)[i]<<" tau21:"<<(*JetsAK8_NsubjettinessTau2)[i]/(*JetsAK8_NsubjettinessTau1)[i]<<" DeepDoubleBDiscr:"<<(*JetsAK8_deepDoubleBDiscriminatorH)[i]<<" WDiscr:"<<(*JetsAK8_wDiscriminatorDeep)[i]<<" WdiscrDecorr:"<<(*JetsAK8_wDiscriminatorDeepDecorrel)[i]<<endl;
  }
}
