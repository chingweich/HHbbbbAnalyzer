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

TCanvas* c1;

TH1D* HHbbbbAnalyzerBase(int wM, string st,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		
		TFile *f;
		TTree *tree;
	    TH1D* th1 =new TH1D("cut flow","",5,0,5);
		
		int nPass[10]={0};
		int total=0;
		
		
		
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
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<200)continue;
				if(fabs(thisHiggs->Eta())>2.5)continue;
				FatjetPreSelection.push_back(i);
			}
		
		    //bool isDeltaEta13=0;
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
	
			
			float*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2"); 
			float* FATjetTau1=data.GetPtrFloat("FATjetTau1");
		    float* FATjetTau2=data.GetPtrFloat("FATjetTau2");
			Int_t* FATnSubSDJet=data.GetPtrInt("FATnSubSDJet");
			vector<float>   *FATsubjetSDCSV = data.GetPtrVectorFloat("FATsubjetSDCSV");
		    vector<float>   *FATsubjetSDPx =  data.GetPtrVectorFloat("FATsubjetSDPx");
		    vector<float>   *FATsubjetSDPy =  data.GetPtrVectorFloat("FATsubjetSDPy");
		    vector<float>   *FATsubjetSDPz =  data.GetPtrVectorFloat("FATsubjetSDPz");
			Float_t*  jetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			
			string fixName ="";if(!isSignal)fixName="";
			
		    vector<float>   *FATsubjetSDE =  data.GetPtrVectorFloat(Form("FATsubjetSD%sE",fixName.data()));
			
			
			
			//bool isMass1000=0;
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
			//cout<<jEntry<<"1"<<endl;
			
			vector<int> FatjetPreSelection3i,FatjetPreSelection3j;
			for(unsigned int i=0;i<FatjetPreSelection2i.size();i++){
				if(jetPRmassL2L3Corr[FatjetPreSelection2i[i]]>135||jetPRmassL2L3Corr[FatjetPreSelection2i[i]]<105)continue;
				if(jetPRmassL2L3Corr[FatjetPreSelection2j[i]]>135||jetPRmassL2L3Corr[FatjetPreSelection2j[i]]<105)continue;
				//if(FATjetPRmass[FatjetPreSelection2i[i]]>150||FATjetPRmass[FatjetPreSelection2i[i]]<90)continue;
				//if(FATjetPRmass[FatjetPreSelection2j[i]]>150||FATjetPRmass[FatjetPreSelection2j[i]]<90)continue;
				
				
				FatjetPreSelection3i.push_back(FatjetPreSelection2i[i]);
				FatjetPreSelection3j.push_back(FatjetPreSelection2j[i]);
			}
			if(FatjetPreSelection3i.size()<1)continue;
			nPass[2]++;
			//cout<<jEntry<<"2"<<endl;
			
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
			
			
			
			
		}
		
		
		}
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th1->SetBinContent(i+1,double(nPass[i])/total);
			th1->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	return th1;
	
}


void HHbbbbAnalyzer(){
	//setNCUStyle();
	
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
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_",
		"",
	};
	string  masspoint[10]={"1000","1200","1400","1600","2000","2500","3000","3500","4000","4500"};
	double eff1[20],eff2[20],eff3[20],eff4[20],eff5[20];
	TH1D* th1[20];
	for(int i=0;i<10;i++)th1[i]=HHbbbbAnalyzerBase(2,st1[i],eff1[i],eff2[i],eff3[i],eff4[i],eff5[i],1);
	th1[10]=HHbbbbAnalyzerBase(155,st1[10],eff1[10],eff2[10],eff3[10],eff4[10],eff5[10]);
	th1[11]=HHbbbbAnalyzerBase(103,st1[11],eff1[11],eff2[11],eff3[11],eff4[11],eff5[11]);
	th1[12]=HHbbbbAnalyzerBase(71,st1[12],eff1[12],eff2[12],eff3[12],eff4[12],eff5[12]);
	
	/*
	
	c1 = new TCanvas("c1","",1360,768);
	
	setNCUStyle(true);
	TLegend* leg ;
	leg=new TLegend(0.791452,0.562447,0.890645,0.783966);
	//leg->SetFillColor(18);
	//leg->SetFillStyle(0);
	//leg->SetTextSize(0.02);
	//leg->SetBorderSize(2);
	
	
	
	TH1D* th2[5];
	for(int i=0;i<5;i++){
		th2[i]=new TH1D("","",13,0,13);
		for(int j=1;j<14;j++){
			if(i==0)th2[i]->SetBinContent(j,eff1[j-1]);
			if(i==1)th2[i]->SetBinContent(j,eff2[j-1]);
			if(i==2)th2[i]->SetBinContent(j,eff3[j-1]);
			if(i==3)th2[i]->SetBinContent(j,eff4[j-1]);
			if(i==4)th2[i]->SetBinContent(j,eff5[j-1]);
			th2[i]->GetXaxis()->SetBinLabel(j,Form("%s",masspoint[j-1].data()));
		}
	}
	//TLegend* leg ;
	//leg=new TLegend(0.691452,0.662447,0.890645,0.883966);
	
	th2[0]->GetXaxis()->SetBinLabel(11,"QCD1000-1500");
	th2[0]->GetXaxis()->SetBinLabel(12,"QCD1500-2000");
	th2[0]->GetXaxis()->SetBinLabel(13,"QCD2000-Inf.");
	
	th2[0]->SetMinimum(0);
	th2[0]->SetLineColor(1);
	th2[1]->SetLineColor(2);
	th2[2]->SetLineColor(3);
	th2[3]->SetLineColor(4);
	th2[4]->SetLineColor(5);
	
	leg->AddEntry(th2[0],"cut1");
	leg->AddEntry(th2[1],"cut2");
	leg->AddEntry(th2[2],"cut3");
	leg->AddEntry(th2[3],"cut4");
	leg->AddEntry(th2[4],"cut5");
	
	th2[0]->Draw();
	c1->Print("pdf/cutM5.png");
	
	th2[0]->Draw();
	//c1->Print("pdf/cutM1.png");
	th2[1]->Draw("same");
	//c1->Print("pdf/cutM2.png");
	th2[2]->Draw("same");
	//c1->Print("pdf/cutM3.png");
	th2[3]->Draw("same");
	//c1->Print("pdf/cutM4.png");
	th2[4]->Draw("same");
	leg->Draw("same");
	c1->Print("pdf/cutM6.png");
    //c1->SetLogy(1);
	
	//c1->Print("pdf/cutM7.png");
	
	*/
	
	ofstream myfile;
	myfile.open ("1105.txt");
	
	for(int i=0;i<14;i++)eff5[i]*=1000;
	eff5[11]*=1064;
	eff5[12]*=121.5;
	eff5[13]*=25.42 ;
	
	
	myfile<<"DATA 0"<<endl;
	myfile<<"QCD "<<eff5[11]+eff5[12]+eff5[13]<<endl;
	for(int i=0;i<11;i++)myfile<<"M"<<masspoint[i]<<" "<<eff5[i]<<endl;
	myfile.close();
}