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
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  bool isFastSim = false;
  // float xsec = 0.0, numEvents = 0.0;
  if(s_data.Contains("TChi")){
    isFastSim = true;
    //   if(s_data.Contains("TChiWZ_1000")){ xsec = 1.34352e-3; numEvents = 28771;}
    //   else if(s_data.Contains("TChiWZ_800")){ xsec = 4.75843e-3; numEvents = 34036;}
    //   cout<<"Assigning xsec as: "<<xsec<<endl;
  }
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
    if (JetsAK8->size() < 2) continue;                                          
    
    Double_t mtbmin = -99999. , mct = -99999. ; //mCT is the contransverse mass variable 
    mtbmin = find_bjets_mtbmin();
    /*
    bjets.resize(0);
    Double_t mtbmin = 0, mct = 0; //mCT is the contransverse mass variable 
    //vector<TLorentzVector> bjets;
    for(int i=0; i<Jets_bJetTagDeepCSVprobb->size(); i++){
      if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
      if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)
	bjets.push_back((*Jets)[i]);
    }
    sortTLorVec(&bjets);
    // cout<<"Bjets Pt: "<< (*Jets)[0].Pt() <<"\t"<< (*Jets)[1].Pt() <<"\t"<< (*Jets)[2].Pt() <<endl;
    
    if(bjets.size()==1){
      mtbmin = sqrt(2*(*Jets)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[0].Phi()))));
      //h_mtbmin->Fill(mtbmin,wt);
    }
    else if(bjets.size() >= 2){ 
      mtbmin = min( sqrt(2*(*Jets)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[0].Phi())))),
		    sqrt(2*(*Jets)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[1].Phi())))) );
      mct =  sqrt(2*(*Jets)[0].Pt()*(*Jets)[1].Pt()*(1+cos(DeltaPhi((*Jets)[0].Phi(),(*Jets)[1].Phi()))));
      // return mtbmin;
      //   h_mtbmin->Fill(mtbmin,wt);
      //h_mct->Fill(mct,wt);
    } 
    */
    
    // Double_t bjets_mTbmin = bjets_mtb();
    if(NJets < 2 || HT < 300 || MHT < 200 || MET <250                           // HT, MET, MHT cuts
       || (*JetsAK8)[0].Pt() < 200 || (*JetsAK8)[1].Pt() < 200                  // pT >200
       || abs((*JetsAK8)[0].Eta()) >2 || abs((*JetsAK8)[1].Eta()) >2            //  |Eta|<2
       || NMuons!=0 || NElectrons!=0                                            // Veto e,muons
       || (MHT/HT > 1.0) || !JetID                                              // MHT<HT and Veto JET ID
       || !(DeltaPhi1 > 1.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3)  //Angle cuts
       //|| BTags!=0                                                                //Veto B jets
       ||mtbmin < 200                                                     // mTbmin > 200 to reduce ttbar bkg
       || isoMuonTracks!=0 || isoElectronTracks!=0 || isoPionTracks!=0) continue; //Veto isolated tracks
    
    //    cout<<"after ra2b"<<endl;
    //----MET
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
    // print(jentry);
    // cout<<"METPhi:"<<METPhi<<" dphi1:"<<dphi1<<" dphi2:"<<dphi2<<" dphi3:"<<dphi3<<" dphi4:"<<dphi4<<endl;
    
    //    if(!(dphi1 > 0.5 && dphi2 > 0.5 && dphi3 > 0.5 && dphi4 > 0.5)) continue;
    //    h_cutflow->Fill("dPhiCuts",wt);

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
    //----MT
    double mt = 0, mt2j = 0;
    if(JetsAK8->size() > 0) mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
    if(JetsAK8->size()>=2) mt2j = sqrt(2*(*JetsAK8)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[1].Phi()))));

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
    
    // cout<<"jets ak8 size "<<JetsAK8->size()<<endl;
    //cout<<"jets disc size "<<JetsAK8_zDiscriminatorDeep->size()<<endl;
    //cout<<"\n";
   
    // SR =1 , Mass_SB =2 , MET_regions =3
    int Region = 2;
    switch (Region){
    case 3:
      {
	// W H 
	if( 
	   ((*JetsAK8_softDropMass)[0] >65) &&
	   ((*JetsAK8_softDropMass)[0] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[0] > 0.9) )
	  {
	    if(         //middle band
	       ((*JetsAK8_softDropMass)[1] > 85) &&
	       ((*JetsAK8_softDropMass)[1] < 135))
	      {
		h_WHMET_RegC->Fill(MET,wt);
		h_WHAK8J2Mass_RegC->Fill((*JetsAK8_softDropMass)[1],wt);
		//cout<<"Fill C WH"<<endl;
		if( (*JetsAK8_doubleBDiscriminator)[1] >0.3)  // Signal region **
		  {
		    h_WHMET_RegA->Fill(MET,wt);
		    h_WHAK8J2Mass_RegA->Fill((*JetsAK8_softDropMass)[1],wt);
		    cout<<"Fill A WH"<<endl;
		  }
	      }
	    if( //side bands
	       (((*JetsAK8_softDropMass)[1] >50) && ((*JetsAK8_softDropMass)[1] < 85)) ||
	       (((*JetsAK8_softDropMass)[1] > 135) && ((*JetsAK8_softDropMass)[1] < 250)) )
	      {
		h_WHMET_RegD->Fill(MET,wt);
		h_WHAK8J2Mass_RegD->Fill((*JetsAK8_softDropMass)[1],wt);
		// cout<<"Fill D WH"<<endl;
		if( (*JetsAK8_doubleBDiscriminator)[1] >0.3)
		  {
		    h_WHMET_RegB->Fill(MET,wt);
		    h_WHAK8J2Mass_RegB->Fill((*JetsAK8_softDropMass)[1],wt);
		    //cout<<"Fill B WH"<<endl;
		  }
	      }
	  }
	//     H W
	if( 
	   ((*JetsAK8_softDropMass)[1] > 65) &&
	   ((*JetsAK8_softDropMass)[1] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[1] > 0.9) )
	  {
	    if(         //middle band or signal band
	       ((*JetsAK8_softDropMass)[0] > 85) &&
	       ((*JetsAK8_softDropMass)[0] < 135))
	      {
		h_HWMET_RegC->Fill(MET,wt);
		h_HWAK8J1Mass_RegC->Fill((*JetsAK8_softDropMass)[0],wt);
	    	//cout<<"Fill C HW"<<endl;
		if( (*JetsAK8_doubleBDiscriminator)[0] >0.3)
		  {
		    h_HWMET_RegA->Fill(MET,wt);
		    h_HWAK8J1Mass_RegA->Fill((*JetsAK8_softDropMass)[0],wt);
		  }
	      }
	    if( //side bands
	       (((*JetsAK8_softDropMass)[0] >50) && ((*JetsAK8_softDropMass)[0] < 85)) ||
	       (((*JetsAK8_softDropMass)[0] >135) && ((*JetsAK8_softDropMass)[0] < 250)) )
	      {
		h_HWMET_RegD->Fill(MET,wt);
		h_HWAK8J1Mass_RegD->Fill((*JetsAK8_softDropMass)[0],wt);
		if( (*JetsAK8_doubleBDiscriminator)[0] >0.3)
		  {
		    h_HWMET_RegB->Fill(MET,wt);
		    h_HWAK8J1Mass_RegB->Fill((*JetsAK8_softDropMass)[0],wt);
		  }
	      }
	  }
       break;
      } // end Mass_SB case
 
    case 2:
      {
	if( // W H
	   ((*JetsAK8_softDropMass)[0] > 65) &&
	   ((*JetsAK8_softDropMass)[0] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[0] > 0.9) &&
	   ((*JetsAK8_softDropMass)[1] > 50) &&
	   ((*JetsAK8_softDropMass)[1] < 250) ){
	  h_WHAK8J1MassNo2bTag->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_WHAK8J2MassNo2bTag->Fill((*JetsAK8_softDropMass)[1],wt); 
	  if (((*JetsAK8_deepDoubleBDiscriminatorH)[1] > 0.45) ){ // Sig reg mass range
	    h_WHAK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt); 
	    h_WHAK8J2Mass->Fill((*JetsAK8_softDropMass)[1],wt); 
	  }
	}
	
	if( // H W
	   ((*JetsAK8_softDropMass)[1] > 65) &&
	   ((*JetsAK8_softDropMass)[1] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[1] > 0.9) &&
	   ((*JetsAK8_softDropMass)[0] > 50) &&
	   ((*JetsAK8_softDropMass)[0] < 250) ){
	  h_HWAK8J1MassNo2bTag->Fill((*JetsAK8_softDropMass)[0],wt); 
	  h_HWAK8J2MassNo2bTag->Fill((*JetsAK8_softDropMass)[1],wt); 
	  if ( ((*JetsAK8_deepDoubleBDiscriminatorH)[0] > 0.45) ){
	    h_HWAK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt); 
	    h_HWAK8J2Mass->Fill((*JetsAK8_softDropMass)[1],wt); 
	  } 
	}
	break;
      }
            
    case 1:
      {
    	//-----------    W H Signal region
	if(
	   ((*JetsAK8_softDropMass)[0] > 65) &&
	   ((*JetsAK8_softDropMass)[0] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[0] > 0.9) &&
	   ((*JetsAK8_softDropMass)[1] > 85) &&
	   ((*JetsAK8_softDropMass)[1] < 135) &&
	   //	   ((*JetsAK8_doubleBDiscriminator)[1] >0.3) ) {
	   ((*JetsAK8_deepDoubleBDiscriminatorH)[1] >0.45) ) {      
	  
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
	   ((*JetsAK8_softDropMass)[1] > 65) &&
	   ((*JetsAK8_softDropMass)[1] < 90) &&
	   ((*JetsAK8_wDiscriminatorDeep)[1] > 0.9) &&
	   ((*JetsAK8_softDropMass)[0] > 85) &&
	   ((*JetsAK8_softDropMass)[0] < 135) &&
	   //	   ((*JetsAK8_doubleBDiscriminator)[0] > 0.3) ) {
	   ((*JetsAK8_deepDoubleBDiscriminatorH)[0] > 0.45) ){
	  
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
	break;
      } // end Signal region case  
      
    } // end switch
    
    // -----------------------     tweaking wdisc score
    bool wDisTweak = false;
    if(wDisTweak){
      if( (65< (*JetsAK8_softDropMass)[0] && (*JetsAK8_softDropMass)[0] < 90) &&
	  (90< (*JetsAK8_softDropMass)[1] && (*JetsAK8_softDropMass)[1] < 135) &&
	  (*JetsAK8_doubleBDiscriminator)[1] >0.3  ){
	if( (*JetsAK8_wDiscriminatorDeep)[0] > 0.779) {
	  h_WHMETa->Fill(MET,wt);
	}
	if( (*JetsAK8_wDiscriminatorDeep)[0] > 0.981) {
	  h_WHMETc->Fill(MET,wt);
	}
      }
      if( (65< (*JetsAK8_softDropMass)[1] && (*JetsAK8_softDropMass)[1] < 90) &&
	  (90< (*JetsAK8_softDropMass)[0] && (*JetsAK8_softDropMass)[0] < 135) &&
	  (*JetsAK8_doubleBDiscriminator)[0] > 0.3)  {
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
	 ((*JetsAK8_softDropMass)[1] > 65 && (*JetsAK8_softDropMass)[1] < 100) && 
	 ((*JetsAK8_softDropMass)[0] > 65 && (*JetsAK8_softDropMass)[0] < 100) ){ 
	
	for (int gen=0; gen < GenParticles_PdgId->size(); gen++){
	  
	  if((*GenParticles_PdgId)[gen] == 5 && abs((*GenParticles_ParentId)[gen]) == 23){ //Z->bb //removing double count by having only 1 b and its parrent Z
	    int GenZIdx = (*GenParticles_ParentIdx)[gen]; 
	    if(abs(DeltaR((*GenParticles)[GenZIdx].Eta(),(*GenParticles)[GenZIdx].Phi(),(*JetsAK8)[0].Eta(),(*JetsAK8)[0].Phi())) < 0.1){ //1st ak8 Zbb match
	      cout<<" ************  In gen reco match"<<endl;
	      h_AK8J1doubleBDis->Fill(((*JetsAK8_doubleBDiscriminator)[0]),wt);
	      h_AK8J1deepdoubleBDis->Fill(((*JetsAK8_deepDoubleBDiscriminatorH)[0]),wt);
	      h_AK8J1Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[0])/((*JetsAK8_NsubjettinessTau1)[0]),wt);
	      h_AK8J1zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[0]),wt);
	    }
	    if(abs(DeltaR((*GenParticles)[GenZIdx].Eta(),(*GenParticles)[GenZIdx].Phi(),(*JetsAK8)[1].Eta(),(*JetsAK8)[1].Phi())) < 0.1){ //2nd ak8 Zbb match
	      h_AK8J2doubleBDis->Fill(((*JetsAK8_doubleBDiscriminator)[1]),wt);
	      h_AK8J2deepdoubleBDis->Fill(((*JetsAK8_deepDoubleBDiscriminatorH)[1]),wt);  
	      h_AK8J2Tau21->Fill(((*JetsAK8_NsubjettinessTau2)[1])/((*JetsAK8_NsubjettinessTau1)[1]),wt);
	      h_AK8J2zDis->Fill(((*JetsAK8_zDiscriminatorDeep)[1]),wt);
	    }
	  } // end 2nd Zbb match
	} // end for
      } // end tune mass window
    } // end if
    //print(jentry);


    nEvtSurv++;
    h_cutflow->Fill("NEvtsNoWtLeft",1);
    //    cout<<" ____ End evt loop ----------   "<<endl;
  } // loop over entries
  cout<<"No. of entries survived: "<<nEvtSurv<<endl;
}

Double_t SignalReg::find_bjets_mtbmin(){
  Double_t mtbmin = 0, mct = 0; //mCT is the contransverse mass variable 
  //vector<TLorentzVector> bjets;
  bjets.resize(0);
  for(int i=0; i < Jets_bJetTagDeepCSVprobb -> size(); i++){
    if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
    if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)
      bjets.push_back((*Jets)[i]);
  }
  sortTLorVec(&bjets); // Sort b jets by decreasing mass
  // cout<<"Bjets Pt: "<< (*Jets)[0].Pt() <<"\t"<< (*Jets)[1].Pt() <<"\t"<< (*Jets)[2].Pt() <<endl;
    
  if(bjets.size()==1){
    mtbmin = sqrt(2*(*Jets)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[0].Phi()))));
    return mtbmin;
    //h_mtbmin->Fill(mtbmin,wt);
  }
  else if(bjets.size() >= 2){ 
    mtbmin = min( sqrt(2*(*Jets)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[0].Phi())))),
		  sqrt(2*(*Jets)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Jets)[1].Phi())))) );
    mct =  sqrt(2*(*Jets)[0].Pt()*(*Jets)[1].Pt()*(1+cos(DeltaPhi((*Jets)[0].Phi(),(*Jets)[1].Phi()))));
    return mtbmin;
    //   h_mtbmin->Fill(mtbmin,wt);
    //h_mct->Fill(mct,wt);
  }
}


void SignalReg::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"MomMass:"<<SusyMotherMass<<" Kid Mass:"<<SusyLSPMass<<endl;
  for(int i=0;i<GenParticles->size();i++){  
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx :"<<(*GenParticles)[i].Px()<<" Py :"<<(*GenParticles)[i].Py()<<" Pz :"<<(*GenParticles)[i].Pz()<<" E: "<<(*GenParticles)[i].Energy()<<" M:"<<(*GenParticles)[i].M()<<endl;
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" E:"<<(*GenParticles)[i].Energy()<<endl;
    
  }
  for(int i=0;i<Photons->size();i++){
    double dR=0;//DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*Photons)[i].Eta(),(*Photons)[i].Phi() );
    //cout<<jentry<<" i:"<<i<<" phoSize:"<<Photons->size()<<" Pt:"<<bestPhoton.Pt()<<" eta:"<<bestPhoton.Eta()<<" phi:"<<bestPhoton.Phi()<<" otherP:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<" dR:"<<dR<<endl;
  }
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  cout<<"MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
}
