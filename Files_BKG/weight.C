#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include"TFile.h"
#include"TH1.h"
#include "TTree.h"


using namespace std;
void weight(string infile){

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

  cout<<fname.size()<<endl;
  //  TFile *f[fname.size()];

  for(int i=0;i< fname.size();i++){
    cout<<fname[i]<<endl;
    //TTree *tr=(TTree*)fname[i]->Get("tree");

//    f[i] = new TFile(fname[i]);
    //TTree* tree;
  
    //tree->Draw("abs(Weight)>>h");

    //cout<<"Weight: "<< fname[i] << h->GetMean()<<endl;
  }
  
}
