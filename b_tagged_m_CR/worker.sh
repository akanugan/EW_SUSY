#!/bin/sh

executable=$1
inputFileTag=$2
outputFileTag=$3
#commitHash=$4
datasetName=$4
currDir=$(pwd)
######################################
# SETUP CMSSW STUFF...
######################################
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
pwd

######################################
# SETUP GRID STUFF...
######################################
export X509_USER_PROXY=/cms/data/$USER/.x509_user_proxy

######################################
# SETUP PRIVATE STUFF...
######################################
echo "ls"
echo "RUNNING ANALYSIS"
./$executable $inputFileTag $outputFileTag $datasetName
echo "processed. ls"

