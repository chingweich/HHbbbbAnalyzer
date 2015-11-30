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
#include "setNCUStyle.C"
#include <TH2.h>




TCanvas* c1;

void HHbbbbAnalyzerBase(int wM, string st,string st2,double Xsec=1,int SignalBkgNum=0){
	
	//for (int massP=0;massP<1;massP++){
		
		TFile *f;
		TTree *tree;
	   
		
		int nPass[10]={0};
		int total=0;
		
		/*
		TH2F* th2f[10];
		TH1D* th1[10];
		TH1D* th2[10];
		th1[0]=new TH1D("1000_FATjetPRmassL2L3Corr","1000_FATjetPRmassL2L3Corr",600,0,600);
		th1[1]=new TH1D("1200_FATjetPRmassL2L3Corr","1200_FATjetPRmassL2L3Corr",600,0,600);
		th1[2]=new TH1D("1400_FATjetPRmassL2L3Corr","1400_FATjetPRmassL2L3Corr",600,0,600);
		th1[3]=new TH1D("1600_FATjetPRmassL2L3Corr","1600_FATjetPRmassL2L3Corr",600,0,600);
		th1[4]=new TH1D("2000_FATjetPRmassL2L3Corr","2000_FATjetPRmassL2L3Corr",600,0,600);
		th1[5]=new TH1D("2500_FATjetPRmassL2L3Corr","2500_FATjetPRmassL2L3Corr",600,0,600);
		th1[6]=new TH1D("3000_FATjetPRmassL2L3Corr","3000_FATjetPRmassL2L3Corr",600,0,600);
		th1[7]=new TH1D("3500_FATjetPRmassL2L3Corr","3500_FATjetPRmassL2L3Corr",600,0,600);
		th1[8]=new TH1D("4000_FATjetPRmassL2L3Corr","4000_FATjetPRmassL2L3Corr",600,0,600);
		th1[9]=new TH1D("4500_FATjetPRmassL2L3Corr","4500_FATjetPRmassL2L3Corr",600,0,600);
		
		th2[0]=new TH1D("1000_FATjetPRmassL2L3Corr_2","1000_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[1]=new TH1D("1200_FATjetPRmassL2L3Corr_2","1200_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[2]=new TH1D("1400_FATjetPRmassL2L3Corr_2","1400_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[3]=new TH1D("1600_FATjetPRmassL2L3Corr_2","1600_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[4]=new TH1D("2000_FATjetPRmassL2L3Corr_2","2000_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[5]=new TH1D("2500_FATjetPRmassL2L3Corr_2","2500_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[6]=new TH1D("3000_FATjetPRmassL2L3Corr_2","3000_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[7]=new TH1D("3500_FATjetPRmassL2L3Corr_2","3500_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[8]=new TH1D("4000_FATjetPRmassL2L3Corr_2","4000_FATjetPRmassL2L3Corr_2",600,0,600);
		th2[9]=new TH1D("4500_FATjetPRmassL2L3Corr_2","4500_FATjetPRmassL2L3Corr_2",600,0,600);
		
		th2f[0]=new TH2F("1000_jetPRmassL2L3Corr","1000_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[1]=new TH2F("1200_jetPRmassL2L3Corr","1200_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[2]=new TH2F("1400_jetPRmassL2L3Corr","1400_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[3]=new TH2F("1600_jetPRmassL2L3Corr","1600_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[4]=new TH2F("2000_jetPRmassL2L3Corr","2000_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[5]=new TH2F("2500_jetPRmassL2L3Corr","2500_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[6]=new TH2F("3000_jetPRmassL2L3Corr","3000_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[7]=new TH2F("3500_jetPRmassL2L3Corr","3500_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[8]=new TH2F("4000_jetPRmassL2L3Corr","4000_jetPRmassL2L3Corr",600,0,600,600,0,600);
		th2f[9]=new TH2F("4500_jetPRmassL2L3Corr","4500_jetPRmassL2L3Corr",600,0,600,600,0,600);
		for(int i=0 ;i< 10;i++){
			th2f[i]->Sumw2();
			th1[i]->Sumw2();
			th2[i]->Sumw2();
		}
		*/
		
		TH1D* th1[10];
		th1[0]=new TH1D("Mass_dEta8","Mass_dEta8",1000,0,5000);
		th1[1]=new TH1D("Mass_dEta9","Mass_dEta9",1000,0,5000);
		th1[2]=new TH1D("Mass_dEta10","Mass_dEta10",1000,0,5000);
		th1[3]=new TH1D("Mass_dEta11","Mass_dEta11",1000,0,5000);
		th1[4]=new TH1D("Mass_dEta12","Mass_dEta12",1000,0,5000);
		th1[5]=new TH1D("Mass_dEta13","Mass_dEta13",1000,0,5000);
		th1[6]=new TH1D("Mass_dEta14","Mass_dEta14",1000,0,5000);
		th1[7]=new TH1D("Mass_dEta15","Mass_dEta15",1000,0,5000);
		th1[8]=new TH1D("dEta","dEta",100,0,4);
		//th1[8]=new TH1D("Mass_dEta8","4000_dEta",400,0,4);
		//th1[9]=new TH1D("Mass_dEta8","4500_dEta",400,0,4);
		for(int i=0 ;i< 9;i++)th1[i]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			f = TFile::Open(Form("%s%d.root",st.data(),w));
		if (!f || !f->IsOpen())continue;
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		
		cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	
		//cout<<data.GetEntriesFast()<<endl;
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2"); 
			float* FATjetTau1=data.GetPtrFloat("FATjetTau1");
		    float* FATjetTau2=data.GetPtrFloat("FATjetTau2");
			Int_t* FATnSubSDJet=data.GetPtrInt("FATnSubSDJet");
			vector<float>   *FATsubjetSDCSV = data.GetPtrVectorFloat("FATsubjetSDCSV");
		    vector<float>   *FATsubjetSDPx =  data.GetPtrVectorFloat("FATsubjetSDPx");
		    vector<float>   *FATsubjetSDPy =  data.GetPtrVectorFloat("FATsubjetSDPy");
		    vector<float>   *FATsubjetSDPz =  data.GetPtrVectorFloat("FATsubjetSDPz");
			Float_t*  jetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
		    vector<float>   *FATsubjetSDE =  data.GetPtrVectorFloat("FATsubjetSDE");
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<200)continue;
				if(fabs(thisHiggs->Eta())>2.4)continue;
				if(FATjetPRmass[i]<105 || FATjetPRmass[i]>135)continue;
				FatjetPreSelection.push_back(i);
			}
			if(FatjetPreSelection.size()<2)continue;
			//cout<<FatjetPreSelection.size()<<endl;
			/*
		    //bool isDeltaEta13=0;
			vector<int> FatjetPreSelection3i,FatjetPreSelection3j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
				    if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.3){
						FatjetPreSelection3i.push_back(FatjetPreSelection[i]);
						FatjetPreSelection3j.push_back(FatjetPreSelection[j]);
					}
				}
				//if(isDeltaEta13)break;
			}
			if(FatjetPreSelection3i.size()<1)continue;
			*/
			

			/*
			
			
			
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				if(FATnSubSDJet[FatjetPreSelection3i[i]]<2||FATnSubSDJet[FatjetPreSelection3j[i]]<2)continue;
				TLorentzVector  FATsubjet_1,FATsubjet_2;
				bool thisHiggsBTagTight=0,thisHiggsBTagLoose=0;
				FATsubjet_1.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection3i[i]][0],FATsubjetSDPy[FatjetPreSelection3i[i]][0],FATsubjetSDPz[FatjetPreSelection3i[i]][0],FATsubjetSDE[FatjetPreSelection3i[i]][0]);
				FATsubjet_2.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection3i[i]][1],FATsubjetSDPy[FatjetPreSelection3i[i]][1],FATsubjetSDPz[FatjetPreSelection3i[i]][1],FATsubjetSDE[FatjetPreSelection3i[i]][1]);
				if(FATsubjet_1.DeltaR(FATsubjet_2)<0.3 && FATjetCISVV2[FatjetPreSelection3i[i]]<0.244)continue;
				if(FATsubjet_1.DeltaR(FATsubjet_2)>0.3 && (FATsubjetSDCSV[FatjetPreSelection3i[i]][0]>0.244 || FATsubjetSDCSV[FatjetPreSelection3i[i]][1]>0.244))thisHiggsBTagLoose=1;
				if(FATsubjet_1.DeltaR(FATsubjet_2)>0.3 && (FATsubjetSDCSV[FatjetPreSelection3i[i]][0]>0.244 && FATsubjetSDCSV[FatjetPreSelection3i[i]][1]>0.244))thisHiggsBTagTight=1;
				if(FATsubjet_1.DeltaR(FATsubjet_2)<0.3){
					thisHiggsBTagTight=1;
					thisHiggsBTagLoose=1;
				}
				FATsubjet_1.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection3j[i]][0],FATsubjetSDPy[FatjetPreSelection3j[i]][0],FATsubjetSDPz[FatjetPreSelection3j[i]][0],FATsubjetSDE[FatjetPreSelection3j[i]][0]);
				FATsubjet_2.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection3j[i]][1],FATsubjetSDPy[FatjetPreSelection3j[i]][1],FATsubjetSDPz[FatjetPreSelection3j[i]][1],FATsubjetSDE[FatjetPreSelection3j[i]][1]);
				TLorentzVector  FATsubjet_3,FATsubjet_4;
				bool thatHiggsBTagTight=0,thatHiggsBTagLoose=0;
				if(FATsubjet_3.DeltaR(FATsubjet_4)<0.3 && FATjetCISVV2[FatjetPreSelection3j[i]]<0.244)continue;
				if(FATsubjet_3.DeltaR(FATsubjet_4)>0.3 && (FATsubjetSDCSV[FatjetPreSelection3j[i]][0]>0.244 || FATsubjetSDCSV[FatjetPreSelection3j[i]][1]>0.244))thatHiggsBTagLoose=1;
				if(FATsubjet_3.DeltaR(FATsubjet_4)>0.3 && (FATsubjetSDCSV[FatjetPreSelection3j[i]][0]>0.244 && FATsubjetSDCSV[FatjetPreSelection3j[i]][1]>0.244))thatHiggsBTagTight=1;
				if(FATsubjet_3.DeltaR(FATsubjet_4)<0.3){
					thatHiggsBTagTight=1;
					thatHiggsBTagLoose=1;
				}
				//if(FATsubjet_3.DeltaR(FATsubjet_4)>0.3)cout<<"Y"<<endl;
				if(thisHiggsBTagLoose==0 || thatHiggsBTagLoose==0 ){
					//cout<<thisHiggsBTagLoose<<","<<thatHiggsBTagLoose<<endl;
					//if(FATsubjet_3.DeltaR(FATsubjet_4)>0.3)cout<<FATsubjetSDCSV[FatjetPreSelection3j[i]][0]<<","<<FATsubjetSDCSV[FatjetPreSelection3j[i]][1]<<endl;
					continue;
					
				}
				if(thisHiggsBTagTight==0 && thatHiggsBTagTight==0) continue;
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				
				
			}
			//if(numHiggsBTagLoose<2||numHiggsBTagTight<1)continue;
            if(FatjetPreSelection4i.size()<1)continue;
			nPass[1]++;
			
			//int numHiggsHP=0,numHiggsLP=0;
			
			bool isPassingTau=0;
			for(unsigned int i=0;i<FatjetPreSelection4i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection4i[i]]/FATjetTau1[FatjetPreSelection4i[i]]<0.5)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection4i[i]]/FATjetTau1[FatjetPreSelection4i[i]]<0.75 && FATjetTau2[FatjetPreSelection4i[i]]/FATjetTau1[FatjetPreSelection4i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
			
				
				if(FATjetTau2[FatjetPreSelection4j[i]]/FATjetTau1[FatjetPreSelection4j[i]]<0.5)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection4j[i]]/FATjetTau1[FatjetPreSelection4j[i]]<0.75 && FATjetTau2[FatjetPreSelection4j[i]]/FATjetTau1[FatjetPreSelection4j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				isPassingTau=1;
			}
			if(isPassingTau==0)continue;
			
			*/
			
					
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
				
				TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
				TLorentzVector mjj=*thisHiggs+*thatHiggs;
				//cout<<mjj.M()<<endl;
				if(mjj.M()<1000)continue;
				th1[8]->Fill(fabs(thisHiggs->Eta()-thatHiggs->Eta()));
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<0.8)th1[0]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<0.9)th1[1]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.0)th1[2]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.1)th1[3]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.2)th1[4]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.3)th1[5]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.4)th1[6]->Fill(mjj.M());
				if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.5)th1[7]->Fill(mjj.M());
				}
			}
			
		}
		
		
		}
		cout<<"entries="<<total<<endl;
		//cout<<double(nPass[3])/total<<endl;
		for(int i=0 ;i< 9;i++){
			//th2f[i]->Scale(3000*Xsec/total);
			th1[i]->Scale(3000*Xsec/total);
			//th2[i]->Scale(3000*Xsec/total);
		}
		TFile* outFile = new TFile(Form("root_files_Op_DeltaEta_shape/%s.root",st2.data()),"recreate");
		for(int i=0 ;i< 9;i++){
			//th2f[i]->Write();
			th1[i]->Write();
			//th2[i]->Write();
		}
		outFile->Close();
		for(int i=0 ;i< 9;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
}


void HHbbbbOptimizer(int a){
	
	
	string st1[20]={
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215300/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215340/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215429/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215514/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215644/0000/NCUGlobalTuples_",
	/*6*/"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215734/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215815/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215903/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215944/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_220034/0000/NCUGlobalTuples_",
	/*11*/
		//"/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220208/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220428/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT500to7000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/0000/NCUGlobalTuples_",
	/*16*/
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_",
	
	};
	
	string  masspoint[10]={"1000","1200","1400","1600","2000","2500","3000","3500","4000","4500"};
	
	TH1D* th1[20];
	if(a==9)for(int i=0;i<10;i++)HHbbbbAnalyzerBase(7,st1[i],masspoint[i],1,1);
	if(a==10)HHbbbbAnalyzerBase(813,st1[10],"QCD100" ,27850000,0);
	if(a==11)HHbbbbAnalyzerBase(516,st1[11],"QCD200" ,1717000,0);
	if(a==12)HHbbbbAnalyzerBase(565,st1[12],"QCD300" ,351300,0);
	if(a==13)HHbbbbAnalyzerBase(496,st1[13],"QCD500" ,31630,0);
	if(a==14)HHbbbbAnalyzerBase(394,st1[14],"QCD700" ,6802,0);
	if(a==15)HHbbbbAnalyzerBase(155,st1[15],"QCD1000",1206,0);
	if(a==16)HHbbbbAnalyzerBase(103,st1[16],"QCD1500",120.4,0);
	if(a==17)HHbbbbAnalyzerBase(71 ,st1[17],"QCD2000",25.24,0);

}
/*
void HHbbbbOptimizer(){
	
	int width [nWidth]={20,25,30,35,40,45,50,55,60};
	int bmin[nBmin]={90,95,100,105,110};
	
	for(int i=0;i<nWidth;i++){
		for(int j=0;j<nBmin;j++){
			HHbbbbOptimizerBase(width[i]+bmin[j],bmin[j]);
		}
	}
}

*/