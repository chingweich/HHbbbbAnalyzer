#include "TFile.h"
#include "TTree.h"
#include "iostream"

void skimB(int w , string st){
  TFile* f;
  f =  TFile::Open(Form("%sNCUGlobalTuples_%d.root",st.data(),w));
  TDirectory * dir= (TDirectory*)f->Get(Form("%sNCUGlobalTuples_%d.root:/tree",st.data(),w));
  TTree* tree;
  dir->GetObject("treeMaker",tree);
  //TTree* newtree = tree->CopyTree("AK8PuppijetSDmass[0]>50","",1000000000,0);
  TFile* outputfile = new TFile(Form("%sskimmed/NCUGlobalTuples_%d.root",st.data(),w),"RECREATE");
  TTree* newtree = tree->CopyTree("FATjetPuppiSDmass[0]>43 && FATjetPuppiSDmass[1]>43 && FATjetP4[0].Pt()>200 && FATjetP4[1].Pt()>200 && fabs(FATjetP4[0].Eta())<2.4 && fabs(FATjetP4[1].Eta())<2.4 && fabs(FATjetP4[0].Eta()-FATjetP4[1].Eta())<1.3 ","",1000000000,0);
  outputfile->cd();
  newtree->Write();
  outputfile->Write();
    
}
void skimTree(int a,int b){

  string st1[40]={
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/B/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/C/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/D/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/E/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/F/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/G/",
    "/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/H/",
  };
  string  fileName[40]={
    "dataB",
    "dataC",
    "dataD",
    "dataE",
    "dataF",
    "dataG",
    "dataH",
  };
  int aa[40]={1045,348,584,497,361,854,24};
  cout<<st1[a]<<b<<".root"<<endl;
  skimB(b,st1[a]);
  
}
