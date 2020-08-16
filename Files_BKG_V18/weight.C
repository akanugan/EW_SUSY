#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include"TFile.h"
#include"TH1.h"
#include "TTree.h"


//using namespace std;
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


  char filename[200];
  for(int i=0;i< fname.size();i++){
    //    cout<<fname[i]<<endl;
    sprintf(filename,"%s",fname[i].c_str()); 
    TFile *_file0 = TFile::Open(filename);
    TTree *tree = (TTree*) _file0->Get("tree");
    tree->Draw("abs(Weight)>>hist");
    TH1F *hist = (TH1F*)gPad->GetPrimitive("hist");
    //std::cout << hist->GetMean() << "+/-" << hist->GetStdDev() << std::endl; 

    // Get shortened file name
    size_t pos = 0;
    std::string s = fname[i];
    std::string token;
    std::string delimiter = "/";
    while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      //std::cout << token << std::endl;
      //std::cout << s << std::endl;
      s.erase(0, pos + delimiter.length());
    }
    //std::cout << s << std::endl;

    printf("%-30s: %3.2e \u00B1% 7.2e\n",s.c_str(),hist->GetMean(),hist->GetStdDev());
    
    
//std::cout <<  << hist->GetStdDev() << shist-> +- GetMean() << histGetPrimitive(histDraw(treeGet(
    //TTree *tr=(TTree*)fname[i]->Get("tree");

//    f[i] = new TFile(fname[i]);
    //TTree* tree;
  
    //tree->Draw("abs(Weight)>>h");

    //cout<<"Weight: "<< fname[i] << h->GetMean()<<endl;
  }
  
}
