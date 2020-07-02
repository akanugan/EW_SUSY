#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<stdlib.h>

using namespace std;
void splitRunList(string infile,int nfPerJob,int batch=1){
  //------------ needed for condor files --------------
  string exeCondor  = "worker2.sh";
  string exePbs     = "worker.sh";
  string exeAna     = "signalReg";
  string datasetAna = "";
  string filesToTransfer = "";
  //---------------------------------------------------
  TString fileName = infile;
  if(fileName.Contains("MC2016")) datasetAna = "MC_2016";
  if(fileName.Contains("MC2017")) datasetAna = "MC_2017";
  if(fileName.Contains("MC2018")) datasetAna = "MC_2018";

  if(fileName.Contains("MET_Run2016")) datasetAna = "2016";
  if(fileName.Contains("MET_Run2017")) datasetAna = "2017";
  if(fileName.Contains("MET_Run2018")) datasetAna = "2018";

  // if(fileName.Contains("TChiWZ_1000_1")) datasetAna = "TChiWZ_1000_1_"+datasetAna;
  // if(fileName.Contains("TChiWZ_800_1")) datasetAna = "TChiWZ_800_1_"+datasetAna;
  // if(fileName.Contains("TChiWZ_600_1")) datasetAna = "TChiWZ_600_1_"+datasetAna;
  // if(fileName.Contains("TChipmWW_1000_100")) datasetAna = "TChipmWW_1000_100_"+datasetAna;
  // if(fileName.Contains("TChipmWW_800_100")) datasetAna = "TChipmWW_800_100_"+datasetAna;
  // if(fileName.Contains("TChipmWW_600_100")) datasetAna = "TChipmWW_600_100_"+datasetAna;

  if(fileName.Contains("TChi")) datasetAna = "TChi_"+datasetAna;
  if(fileName.Contains("fast")) datasetAna = "TChi_"+datasetAna;
  //---------------------------------------------------
  if (batch==1)
    cout<<"executable at worker node : "<<exeCondor<<endl
	<<"Analysis executable : "<<exeAna<<endl
	<<"dataset name for analysis : "<<datasetAna<<endl;
  else if (batch==2)
    cout<<"executable at worker node : "<<exePbs<<endl
	<<"Analysis executable : "<<exeAna<<endl
	<<"dataset name for analysis : "<<datasetAna<<endl;
  else {
    cout << "batch should be either 1 (condor) or 2 (pbs)" << endl;
    return;
  }
  //----------------- split the input files into smaller ones ------------------
  ifstream file(infile);
  if(!file){cout<<"Couldn't Open File "<<infile<<endl;}
  string str,dataset=infile;
  dataset.pop_back();  dataset.pop_back();  dataset.pop_back();  dataset.pop_back();
  vector<string> fname;
  while (std::getline(file, str))
    {
      fname.push_back(str);
    }
  file.close();

  int jobid=0;
  char name[1000];
  ofstream outf;
  for(int i=0,j=0;i<fname.size();){
    sprintf(name,"FileList_%s_job%i.txt",dataset.c_str(),jobid);
    outf.open(name);
    for(j=0;j<nfPerJob && i<fname.size();j++){
      outf<<fname[i]<<endl;
      i++;
    }
    jobid++;
    outf.close();
  }

  //--------------------- make files for condor (make dummy .jdl for pbs) ------------------------------------
  char fileListName[200],logFile[200],jobName[100],command[100];
  for(int i=0;i<jobid;i++){
    sprintf(name,"%s_job%i.jdl",dataset.c_str(),i);
    if (batch==2){ // pbs
      sprintf(command,"touch %s",name); // creating dummy .jdl file for findFailedJobs.C
      system(command);
    }
    if (batch==1){ // condor    
      sprintf(fileListName,"FileList_%s_job%i.txt",dataset.c_str(),i);
      sprintf(logFile,"%s_job%i",dataset.c_str(),i);
      outf.open(name);
      outf<<"universe = vanilla"<<endl
	  <<"Executable = "<<exeCondor<<endl
	  <<"Should_Transfer_Files = YES"<<endl
	  <<"WhenToTransferOutput = ON_EXIT_OR_EVICT"<<endl
	  <<"Transfer_Input_Files = "<<filesToTransfer<<","<<exeAna<<","<<fileListName<<endl
	  <<"Output = "<<logFile<<".stdout"<<endl
	  <<"Error = "<<logFile<<".stderr"<<endl
	  <<"Log = "<<logFile<<".condor"<<endl
	  <<"Arguments = "<<exeAna<<" "<<fileListName<<" "<<logFile<<".root "<<datasetAna<<endl
	  <<"+LENGTH=\"SHORT\""<<endl
	  <<endl
	  <<"Queue 1";
      outf.close();
    }    
  }
  //------------------------ submit to condor or pbs --------------------------------------
  int t1=100;
  cout<<"Do you want to submit "<<jobid<<" jobs? If yes enter 100"<<endl;
  //  cin>>t1;
  for(int i=0;i<jobid && t1==100;i++){
    if (batch==1){
      sprintf(name,"condor_submit %s_job%i.jdl",dataset.c_str(),i);
      system(name);
    } else {
      sprintf(fileListName,"FileList_%s_job%i.txt",dataset.c_str(),i);
      sprintf(jobName,"%s_job%i",dataset.c_str(),i);
      sprintf(name,"qsub -N %s -o %s.stdout -e %s.stderr -v exeAna=%s,fileListName=%s,outputFile=%s.root,datasetAna=%s submit.pbs",
	      jobName,jobName,jobName,
	      exeAna.c_str(),fileListName,jobName,datasetAna.c_str());
      cout << name << endl;
      system(name);
    }
  }
  
}

void splitRunList(){
  cout<<"Please specify the input txt file to use and no. of files per job"<<endl;
  cout<<"splitRunList(string infile,int nfPerJob)"<<endl;
}
