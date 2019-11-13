#!/bin/sh
for i in 1
do
    for j in 600
    do
	echo "root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/Autumn18Fast.TChiWZ_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root" > smallrunList.txt
	./signalReg smallrunList.txt Autumn18Fast.TChiWZ_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root TChiWZ_${j}_${i}
	mv Autumn18Fast.TChiWZ_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root /eos/uscms/store/user/vhegde/RA2bTrees_V17/
    done
done