void checkWeights(string arg){

  char filename[200];
  sprintf(filename,"%s",arg.c_str());

  TFile *_file0 = TFile::Open(filename);
  TTree *tree = (TTree*) _file0->Get("tree");  
  tree->Draw("abs(Weight)>>hist");
  TH1F *hist = (TH1F*)gPad->GetPrimitive("hist");
  //std::cout << hist->GetMean() << " +- " << hist->GetStdDev() << std::endl; 

  // Get shortened file name
  size_t pos = 0;
  std::string s=arg;
  std::string token;
  std::string delimiter = "/";
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    //std::cout << s << std::endl;
    s.erase(0, pos + delimiter.length());
  }
  //std::cout << s << std::endl;

  printf("%-30s: %10.5e +- %10.5e\n",s.c_str(),hist->GetMean(),hist->GetStdDev());

}
