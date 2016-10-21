#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>

TFile *f;
TTree *tree;

void dataPrinter(){
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	for (int massP=0;massP<1;massP++){
		
		//signal
		//f = TFile::Open(Form("/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root",masspoint[massP].data()));
		//if (!f || !f->IsOpen())continue;
		//TDirectory * dir = (TDirectory*)f->Get(Form("/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root:/tree",masspoint[massP].data()));
		//for (int w=0;w<1;w++){
			
		
		//QCD  1000-1500
		for (int w=0;w<1;w++){
		//f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084507/0000/NCUGlobalTuples_%d.root",w));
		//if (!f || !f->IsOpen())continue;
		//TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084507/0000/NCUGlobalTuples_%d.root:/tree",w));
		f=TFile::Open("/data7/chchen/8_0_20/B/B1200.root");
		TDirectory * dir = (TDirectory*)f->Get("/data7/chchen/8_0_20/B/B1200.root:/tree");
		//QCD  1500-2000
		//for (int w=0;w<103;w++){
		//f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084638/0000/NCUGlobalTuples_%d.root",w));
		//if (!f || !f->IsOpen())continue;
		//TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084638/0000/NCUGlobalTuples_%d.root:/tree",w));
		
		//QCD  2000-inf
		//for (int w=0;w<71;w++){
		//f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084730/0000/NCUGlobalTuples_%d.root",w));
		//if (!f || !f->IsOpen())continue;
		//TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084730/0000/NCUGlobalTuples_%d.root:/tree",w));
		
		
		
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		data.Print();
		int nPass=0;
		data.GetEntry(1);
		std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		 	const Int_t nsize = data.GetPtrStringSize();
		for(int it=0; it< nsize; it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				cout<<it<<"="<<thisTrig<<endl;
				// std::cout << thisTrig << " : " << results << std::endl;
				
			}
		cout<<"entries="<<data.GetEntriesFast()<<",nPass="<<nPass<<endl;
		
		
		}
		
	}
	
}