#!/bin/sh

#root -l -q 'findFailedJobs.C("TTJets_DiLept_MC_MC2018")'
#root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2018")'
#root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2018")'
#root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2018")'
#root -l -q 'findFailedJobs.C("QCD_HT_MC2018")'
#root -l -q 'findFailedJobs.C("ST__MC2018")'
root -l -q 'findFailedJobs.C("Rare_MC2018")'
#
root -l -q 'findFailedJobs.C("TChiWZ_1000_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_400_MC2018")'

root -l -q 'findFailedJobs.C("TChiWH_1000_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_400_MC2018")'

root -l -q 'findFailedJobs.C("QCD_HT_MC2016")'
root -l -q 'findFailedJobs.C("QCD_HT_MC2017")'
root -l -q 'findFailedJobs.C("QCD_HT_MC2018")'

root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2018")'

root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2016")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2017")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2018")'

root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2016")'
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2017")'
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC2018")'
root -l -q 'findFailedJobs.C("TTJets_HT_MC2016")'
root -l -q 'findFailedJobs.C("TTJets_HT_MC2017")'
root -l -q 'findFailedJobs.C("TTJets_HT_MC2018")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2016")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2017")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2018")'

hadd -f TTJets_MC2016.root TTJets_HT_MC2016.root TTJets_DiLept_MC2016.root TTJets_SingleLeptFromT_MC2016.root 
hadd -f TTJets_MC2017.root TTJets_HT_MC2017.root TTJets_DiLept_MC2017.root TTJets_SingleLeptFromT_MC2017.root
hadd -f TTJets_MC2018.root TTJets_HT_MC2018.root TTJets_DiLept_MC2018.root TTJets_SingleLeptFromT_MC2018.root

root -l -q 'findFailedJobs.C("ST__MC2016")'
root -l -q 'findFailedJobs.C("ST__MC2017")'
root -l -q 'findFailedJobs.C("ST__MC2018")'

hadd -f Rare_MCRun2.root Rare_MC2018.root
hadd -f QCD_HT_MCRun2.root QCD_HT_MC2016.root QCD_HT_MC2017.root QCD_HT_MC2018.root
hadd -f WJetsToLNu_HT_MCRun2.root WJetsToLNu_HT_MC2016.root WJetsToLNu_HT_MC2017.root WJetsToLNu_HT_MC2018.root
hadd -f ZJetsToNuNu_HT_MCRun2.root ZJetsToNuNu_HT_MC2016.root ZJetsToNuNu_HT_MC2017.root ZJetsToNuNu_HT_MC2018.root
hadd -f TTJets_MCRun2.root TTJets_MC2016.root TTJets_MC2017.root TTJets_MC2018.root
hadd -f ST__MCRun2.root ST__MC2016.root ST__MC2017.root ST__MC2018.root


#root -l -q 'findFailedJobs.C("MET_Run2016")'
#root -l -q 'findFailedJobs.C("MET_Run2017")'
#root -l -q 'findFailedJobs.C("MET_Run2018")'
