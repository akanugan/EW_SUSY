# EW SUSY

## CMSSW and ROOT
Setup CMSSW_10_3_1 once. You need it for using ROOT only. Other CMSSW versions may not work because of different ROOT versions.

```
cmsrel CMSSW_10_3_1
cd CMSSW_10_3_1/src/
cmsenv
```

## Using the code
```
git clone https://github.com/akanugan/EW_SUSY
cd EW_SUSY
```

Main analysis code with cuts and sample specific options: SignalReg.cc

Declare and book histograms in SignalReg.h

Tree variables are defined in NtupleVariables.h. Some usefule functions are defined in NtupleVariables.cc.

Specify the file name on which you want to run the code in a text file.

If you modify the code, compile it and then run it.

```
make
./signalReg smallrunList.txt simpleTest.root MC_2018
```
These create a file named `simpleTest.root` with histograms in it. `smallrunList.txt` file contains name of the root file. Last arguement is needed to process the type of sample. For signal sample, it will be TCHiWZ_1000_100_MC_2018 etc.

## Running on condor

To submit jobs to condor for making histograms for all samples and signal:
```
voms-proxy-init --voms cms
./submitMany.sh
```
or submit jobs to pbs (recommended)
```
voms-proxy-init --voms cms
./submitManyPbs.sh
```

Once the jobs are finished,
```
./Check_FailedJobsMany.sh
```

Take a look at some `.stdout` and `.stderr` files and see if the jobs were successful. You may ignore erros like, "Unkown branch XXX..." if you are not using those branches.


## Making plots

Now you have histograms for BG and signal. Make stacked plots:

```
root -l -b -q plotKinStack.C
```

This creates a pdf file with stacked BG and signal. If you want to plot some other variables, open plotKinStack.C file and edit `varName`, `xLabel`, `rebin`.

## To clean the directory off many condor files

Use  ```./cleaUpBatchFiles.sh```



