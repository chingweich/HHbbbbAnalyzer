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

void HHbbbbAnalyzer(){
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	for (int massP=0;massP<11;massP++){
		
		bool isSignal=0;
		
		//signal
		f = TFile::Open(Form("/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root",masspoint[massP].data()));
		if (!f || !f->IsOpen())continue;isSignal=1;
		TDirectory * dir = (TDirectory*)f->Get(Form("/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-%s_13TeV-madgraph.root:/tree",masspoint[massP].data()));
		for (int w=0;w<1;w++){
			
		
		//QCD  1000-1500
		//for (int w=0;w<155;w++){
		//f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084507/0000/NCUGlobalTuples_%d.root",w));
		//if (!f || !f->IsOpen())continue;
		//TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0830/150831_084507/0000/NCUGlobalTuples_%d.root:/tree",w));
		
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
		
		int nPass=0;
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				if(FATjetPRmass[i]>135||FATjetPRmass[i]<110)continue;
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<30)continue;
				if(thisHiggs->Eta()>2.5)continue;
				
				FatjetPreSelection.push_back(i);
			}
			
			float*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2"); 
			float* FATjetTau1=data.GetPtrFloat("FATjetTau1");
		    float* FATjetTau2=data.GetPtrFloat("FATjetTau2");
			Int_t* FATnSubSDJet=data.GetPtrInt("FATnSubSDJet");
			vector<float>   *FATsubjetSDCSV = data.GetPtrVectorFloat("FATsubjetSDCSV");
		    vector<float>   *FATsubjetSDPx =  data.GetPtrVectorFloat("FATsubjetSDPx");
		    vector<float>   *FATsubjetSDPy =  data.GetPtrVectorFloat("FATsubjetSDPy");
		    vector<float>   *FATsubjetSDPz =  data.GetPtrVectorFloat("FATsubjetSDPz");
			
			string originFATsubjetSDE="FATsubjetSDE";
			string signalFATsubjetSDE="FATsubjetSDCE";
			
			if (isSignal) vector<float>   *FATsubjetSDE =  data.GetPtrVectorFloat(signalFATsubjetSDE.data());
		    else vector<float>   *FATsubjetSDE =  data.GetPtrVectorFloat(originFATsubjetSDE.data());
			
			int thisHiggsNum=0,thatHiggsNum=0;
			
			
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				if(FATnSubSDJet[FatjetPreSelection[i]]<2)continue;
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
					TLorentzVector mjj=*thisHiggs+*thatHiggs;
					if(mjj.M()<890)continue;
					if(thisHiggs->DeltaPhi(*thatHiggs)<1.3)continue;
					if(FATnSubSDJet[FatjetPreSelection[j]]<2)continue;
					TLorentzVector  FATsubjet_1,FATsubjet_2;
					FATsubjet_1.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[i]][0],FATsubjetSDPy[FatjetPreSelection[i]][0],FATsubjetSDPz[FatjetPreSelection[i]][0],FATsubjetSDE[FatjetPreSelection[i]][0]);
					FATsubjet_2.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[i]][1],FATsubjetSDPy[FatjetPreSelection[i]][1],FATsubjetSDPz[FatjetPreSelection[i]][1],FATsubjetSDE[FatjetPreSelection[i]][1]);
					bool subjet1PassCSV=1,subjet2PassCSV=1;
					if(FATsubjet_1.DeltaR(FATsubjet_2)<0.3 && FATjetCISVV2[FatjetPreSelection[i]]<0.244)subjet1PassCSV=0;
					if(FATsubjet_1.DeltaR(FATsubjet_2)>0.3 && (FATsubjetSDCSV[FatjetPreSelection[i]][0]<0.244 || FATsubjetSDCSV[FatjetPreSelection[i]][1]<0.244))subjet1PassCSV=0;
					
					
					FATsubjet_1.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[j]][0],FATsubjetSDPy[FatjetPreSelection[j]][0],FATsubjetSDPz[FatjetPreSelection[j]][0],FATsubjetSDE[FatjetPreSelection[j]][0]);
					FATsubjet_2.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[j]][1],FATsubjetSDPy[FatjetPreSelection[j]][1],FATsubjetSDPz[FatjetPreSelection[j]][1],FATsubjetSDE[FatjetPreSelection[j]][1]);
					if(FATsubjet_1.DeltaR(FATsubjet_2)<0.3 && FATjetCISVV2[FatjetPreSelection[i]]<0.244)subjet2PassCSV=0;
					if(FATsubjet_1.DeltaR(FATsubjet_2)>0.3 && (FATsubjetSDCSV[FatjetPreSelection[i]][0]<0.244 || FATsubjetSDCSV[FatjetPreSelection[i]][1]<0.244))subjet2PassCSV=0;
					if(subjet1PassCSV==0 &&  subjet2PassCSV==0 )continue;
					
					thisHiggsNum=FatjetPreSelection[i];
					thatHiggsNum=FatjetPreSelection[j];
					break;
					
					
					
				}
				if(thisHiggsNum!=0)break;
			}
			
			if(thisHiggsNum ==0 || thatHiggsNum==0) continue;
			nPass++;
			
			
		}
		
		cout<<"entries="<<data.GetEntriesFast()<<",nPass="<<nPass<<endl;
		
		
		}
		
	}
	
}