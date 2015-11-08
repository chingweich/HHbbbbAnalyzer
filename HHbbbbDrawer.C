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

void HHbbbbAnalyzerBase(int wM, string st,string st2,int SignalBkgNum=0){
	
	//for (int massP=0;massP<1;massP++){
		
		TFile *f;
		TTree *tree;
	    
		
		int nPass[10]={0};
		int total=0;
		
		
		TH1D* th1=new TH1D("deltaEta","deltaEta",30,0,3);
			TH1D* th2=new TH1D("Mjj","Mjj",2000,400,4000);
			TH1D* th3=new TH1D("SDMass1","SDMass1",200,0,200);
			TH1D* th4=new TH1D("SDMass2","SDMass2",200,0,200);
			TH1D* th5=new TH1D("PRMass1","PRMass1",200,0,200);
			TH1D* th6=new TH1D("PRMass2","PRMass2",200,0,200);
			TH1D* th7=new TH1D("deltaR1","deltaR1",20,0,1);
			TH1D* th8=new TH1D("deltaR2","deltaR2",20,0,1);
			TH1D* th9=new TH1D("CSV1","CSV1",20,0,1);
			TH1D* th10=new TH1D("CSV2","CSV2",20,0,1);
			TH1D* th11=new TH1D("CSVSum","CSVSum",20,0,2);
			TH1D* th12=new TH1D("CSVMax","CSVMax",20,0,1);
			TH1D* th13=new TH1D("CSVMin","CSVMin",20,0,1);
			TH1D* th14=new TH1D("tau21_1","tau21_1",20,0,1);
			TH1D* th15=new TH1D("tau21_2","tau21_2",20,0,1);
			TH1D* th16=new TH1D("jetPt1","jetPt1",1500,0,1500);
			TH1D* th17=new TH1D("jetPt2","jetPt2",1500,0,1500);
			
			TH1F* thf=new TH1F("FATjetPRmassL2L3Corr","FATjetPRmassL2L3Corr",600,0,600);
			TH1F* thf2=new TH1F("FATjetPRmassL2L3Corr_2","FATjetPRmassL2L3Corr_2",600,0,600);
			
			TH1D* cs1=new TH1D("case1","case1",2000,400,6000);
			TH1D* cs2=new TH1D("case2","case2",2000,400,6000);
			TH1D* cs3=new TH1D("case3","case3",2000,400,6000);
		
			TH2F* th2f=new TH2F("2d","2d",600,0,600,600,0,600);
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
		
		
		
			
			

		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*   FATjetSDmass= data.GetPtrFloat("FATjetSDmass");
			float*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2"); 
			float* FATjetTau1=data.GetPtrFloat("FATjetTau1");
		    float* FATjetTau2=data.GetPtrFloat("FATjetTau2");
			Int_t* FATnSubSDJet=data.GetPtrInt("FATnSubSDJet");
			vector<float>   *FATsubjetSDCSV = data.GetPtrVectorFloat("FATsubjetSDCSV");
		    vector<float>   *FATsubjetSDPx =  data.GetPtrVectorFloat("FATsubjetSDPx");
		    vector<float>   *FATsubjetSDPy =  data.GetPtrVectorFloat("FATsubjetSDPy");
		    vector<float>   *FATsubjetSDPz =  data.GetPtrVectorFloat("FATsubjetSDPz");
			Float_t*  jetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			
			string fixName ="";if(SignalBkgNum==1)fixName="";
			
		    vector<float>   *FATsubjetSDE =  data.GetPtrVectorFloat(Form("FATsubjetSD%sE",fixName.data()));
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<200)continue;
				if(fabs(thisHiggs->Eta())>2.5)continue;
				FatjetPreSelection.push_back(i);
			}
			if(FatjetPreSelection.size()<2)continue;
			
			//TH1D* th1=new TH1D("","",,0,);
			
			TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[0]);
			TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[1]);
			TLorentzVector mjj=*thisHiggs+*thatHiggs;
			th1->Fill(fabs(fabs(thisHiggs->Eta())-fabs(thatHiggs->Eta())));
			th2->Fill(mjj.M());
			th3->Fill(FATjetSDmass[FatjetPreSelection[0]]);
			th4->Fill(FATjetSDmass[FatjetPreSelection[1]]);
			th5->Fill(FATjetPRmass[FatjetPreSelection[0]]);
			th6->Fill(FATjetPRmass[FatjetPreSelection[1]]);
			
			TLorentzVector  FATsubjet_1,FATsubjet_2,FATsubjet_3,FATsubjet_4;;
			FATsubjet_1.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[0]][0],FATsubjetSDPy[FatjetPreSelection[0]][0],FATsubjetSDPz[FatjetPreSelection[0]][0],FATsubjetSDE[FatjetPreSelection[0]][0]);
			FATsubjet_2.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[0]][1],FATsubjetSDPy[FatjetPreSelection[0]][1],FATsubjetSDPz[FatjetPreSelection[0]][1],FATsubjetSDE[FatjetPreSelection[0]][1]);
			FATsubjet_3.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[1]][0],FATsubjetSDPy[FatjetPreSelection[1]][0],FATsubjetSDPz[FatjetPreSelection[1]][0],FATsubjetSDE[FatjetPreSelection[1]][0]);
			FATsubjet_4.SetPxPyPzE(FATsubjetSDPx[FatjetPreSelection[1]][1],FATsubjetSDPy[FatjetPreSelection[1]][1],FATsubjetSDPz[FatjetPreSelection[1]][1],FATsubjetSDE[FatjetPreSelection[1]][1]);
				
				
			th7->Fill(FATsubjet_1.DeltaR(FATsubjet_2));
			th8->Fill(FATsubjet_3.DeltaR(FATsubjet_4));
			
			
			th9->Fill(FATjetCISVV2[FatjetPreSelection[0]]);
			th10->Fill(FATjetCISVV2[FatjetPreSelection[1]]);
			
			th11->Fill(FATjetCISVV2[FatjetPreSelection[0]]+FATjetCISVV2[FatjetPreSelection[1]]);
			double temp=FATjetCISVV2[FatjetPreSelection[0]]>FATjetCISVV2[FatjetPreSelection[1]]?FATjetCISVV2[FatjetPreSelection[0]]:FATjetCISVV2[FatjetPreSelection[1]];
			th12->Fill(temp);
			temp=FATjetCISVV2[FatjetPreSelection[0]]<FATjetCISVV2[FatjetPreSelection[1]]?FATjetCISVV2[FatjetPreSelection[0]]:FATjetCISVV2[FatjetPreSelection[1]];
			th13->Fill(temp);
			
			
			
			th14->Fill(FATjetTau2[FatjetPreSelection[0]]/FATjetTau1[FatjetPreSelection[0]]);
			th15->Fill(FATjetTau2[FatjetPreSelection[1]]/FATjetTau1[FatjetPreSelection[1]]);
			
			th16->Fill(thisHiggs->Pt());
			//cout<<thisHiggs->Pt()<<endl;
			th17->Fill(thatHiggs->Pt());
			

			vector<int> FatjetPreSelection1i,FatjetPreSelection1j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
				    if(fabs(fabs(thisHiggs->Eta())-fabs(thatHiggs->Eta()))<1.3){
						FatjetPreSelection1i.push_back(FatjetPreSelection[i]);
						FatjetPreSelection1j.push_back(FatjetPreSelection[j]);
					}
				}
				//if(isDeltaEta13)break;
			}
			if(FatjetPreSelection1i.size()<1)continue;
	
			nPass[0]++;
			
			vector<int> FatjetPreSelection2i,FatjetPreSelection2j;
			for(unsigned int i=0;i<FatjetPreSelection1i.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection1i[i]);
				TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection1j[i]);
				TLorentzVector mjj=*thisHiggs+*thatHiggs;
				//cout<<mjj.M()<<endl;
				if(mjj.M()>1000){
					//cout<<"y"<<endl;
					FatjetPreSelection2i.push_back(FatjetPreSelection1i[i]);
				    FatjetPreSelection2j.push_back(FatjetPreSelection1j[i]);
					
				}
				
				
			}
			if(FatjetPreSelection2i.size()<1)continue;
			nPass[1]++;
			
			vector<int> FatjetPreSelection3i,FatjetPreSelection3j;
			for(unsigned int i=0;i<FatjetPreSelection2i.size();i++){
				//if(jetPRmassL2L3Corr[FatjetPreSelection2i[i]]>135||jetPRmassL2L3Corr[FatjetPreSelection2i[i]]<105)continue;
				//if(FATjetPRmass[FatjetPreSelection2j[i]]>130||FATjetPRmass[FatjetPreSelection2j[i]]<95)continue;
				//if(FATjetPRmass[FatjetPreSelection2i[i]]>150||FATjetPRmass[FatjetPreSelection2i[i]]<90)continue;
				//if(FATjetPRmass[FatjetPreSelection2j[i]]>150||FATjetPRmass[FatjetPreSelection2j[i]]<90)continue;
				thf->Fill(jetPRmassL2L3Corr[FatjetPreSelection2i[i]]);
				thf2->Fill(jetPRmassL2L3Corr[FatjetPreSelection2j[i]]);
				th2f->Fill(jetPRmassL2L3Corr[FatjetPreSelection2i[i]],jetPRmassL2L3Corr[FatjetPreSelection2j[i]]);
				
				
				
				FatjetPreSelection3i.push_back(FatjetPreSelection2i[i]);
				FatjetPreSelection3j.push_back(FatjetPreSelection2j[i]);
			}
			if(FatjetPreSelection3i.size()<1)continue;
			nPass[2]++;
			
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
			nPass[3]++;
			
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
			
  
			nPass[4]++;
			
			TLorentzVector* thisHiggsCs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4i[0]);
			TLorentzVector* thatHiggsCs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4j[0]);
			TLorentzVector mjjCs=*thisHiggsCs+*thatHiggsCs;
			cs1->Fill(mjjCs.M());
			
			
			thatHiggsCs ->SetPxPyPzE(thatHiggsCs->Px()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thatHiggsCs->Py()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thatHiggsCs->Pz()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thatHiggsCs->E()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]]);
			
			mjjCs=*thisHiggsCs+*thatHiggsCs;
			cs2->Fill(mjjCs.M());
			
			thisHiggsCs ->SetPxPyPzE(thisHiggsCs->Px()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thisHiggsCs->Py()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thisHiggsCs->Pz()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]],thisHiggsCs->E()*125/jetPRmassL2L3Corr[FatjetPreSelection4j[0]]);
			
			mjjCs=*thisHiggsCs+*thatHiggsCs;
			cs3->Fill(mjjCs.M());
		}
		
		
		}
		
		double Xsec=1;
		if(SignalBkgNum==2)Xsec=1064;
		else if (SignalBkgNum==3)Xsec=121.5;
		else if (SignalBkgNum==4)Xsec=25.42;
		else if (SignalBkgNum==5)Xsec=6524;
		thf->Scale(1000*Xsec/total);
		thf2->Scale(1000*Xsec/total);
		th1->Scale(1000*Xsec/total);
	th2->Scale(1000*Xsec/total);
	th3->Scale(1000*Xsec/total);
	th4->Scale(1000*Xsec/total);
	th5->Scale(1000*Xsec/total);
	th6->Scale(1000*Xsec/total);
	th7->Scale(1000*Xsec/total);
	th8->Scale(1000*Xsec/total);
	th9->Scale(1000*Xsec/total);
	th10->Scale(1000*Xsec/total);
	th11->Scale(1000*Xsec/total);
	th12->Scale(1000*Xsec/total);
	th13->Scale(1000*Xsec/total);
	th14->Scale(1000*Xsec/total);
	th15->Scale(1000*Xsec/total);
	th16->Scale(1000*Xsec/total);
	th17->Scale(1000*Xsec/total);
	
	cs1->Scale(1000*Xsec/total);
	cs2->Scale(1000*Xsec/total);
	cs3->Scale(1000*Xsec/total);
		
	TFile* outFile = new TFile(Form("root_files/%s.root",st2.data()),"recreate");   
	th1->Write();
	th2->Write();
	th3->Write();
	th4->Write();
	th5->Write();
	th6->Write();
	th7->Write();
	th8->Write();
	th9->Write();
	th10->Write();
	th11->Write();
	th12->Write();
	th13->Write();
	th14->Write();
	th15->Write();
	th16->Write();
	th17->Write();
	cs1->Write();
	cs2->Write();
	cs3->Write();
	
	thf->Write();
	thf2->Write();
	th2f->Write();
	cout<<th2f->GetCorrelationFactor();
	outFile->Close();
	
	
	th1->Clear();
	th2->Clear();
	th3->Clear();
	th4->Clear();
	th5->Clear();
	th6->Clear();
	th7->Clear();
	th8->Clear();
	th9->Clear();
	th10->Clear();
	th11->Clear();
	th12->Clear();
	th13->Clear();
	th14->Clear();
	th15->Clear();
	th16->Clear();
	th17->Clear();
	cs1->Clear();
	cs2->Clear();
	cs3->Clear();
	
	thf->Clear();
	thf2->Clear();
	th2f->Clear();
}


void HHbbbbDrawer(){
	setNCUStyle();
	
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
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_",
		};
	
	string  masspoint[10]={"1000","1200","1400","1600","2000","2500","3000","3500","4000","4500"};
	
	TH1D* th1[20];
	for(int i=0;i<10;i++)HHbbbbAnalyzerBase(7,st1[i],masspoint[i],1);
	HHbbbbAnalyzerBase(155,st1[11],"QCD1000",2);
	HHbbbbAnalyzerBase(103,st1[12],"QCD1500",3);
	HHbbbbAnalyzerBase(71,st1[13],"QCD2000",4);
	HHbbbbAnalyzerBase(394,st1[10],"QCD700",5);
	
	
	
	//c1 = new TCanvas("c1","",1360,768);
	
	//setNCUStyle(true);
	
	
}