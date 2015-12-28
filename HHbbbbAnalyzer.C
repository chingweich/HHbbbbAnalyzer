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
#include "jetEnergyScale.h"

TCanvas* c1;

TH1D* HHbbbbAnalyzerBase(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		int nPass[10]={0};
		int total=0;
		
		TH1D * th1[30];
		th1[0]=new TH1D("cat0","cat0",4000,1000,5000);
		th1[1]=new TH1D("cat1","cat1",4000,1000,5000);
		th1[2]=new TH1D("cat2","cat2",4000,1000,5000);
		th1[3]=new TH1D("cat3","cat3",4000,1000,5000);
		th1[4]=new TH1D("cat4","cat4",4000,1000,5000);
		th1[5]=new TH1D("dRcat0","dRcat0",4000,1000,5000);
		th1[6]=new TH1D("dRcat1","dRcat1",4000,1000,5000);
		th1[7]=new TH1D("dRcat2","dRcat0",4000,1000,5000);
		th1[8]=new TH1D("dRcat3","dRcat1",4000,1000,5000);
		th1[9]=new TH1D("dRcat4","dRcat1",4000,1000,5000);
		th1[10]=new TH1D("dR_com_cat0","dR_com_cat0",4000,1000,5000);
		th1[11]=new TH1D("dR_com_cat1","dR_com_cat1",4000,1000,5000);
		th1[12]=new TH1D("dR_com_cat2","dR_com_cat0",4000,1000,5000);
		th1[13]=new TH1D("dR_com_cat3","dR_com_cat1",4000,1000,5000);
		th1[14]=new TH1D("dR_com_cat4","dR_com_cat1",4000,1000,5000);
		
		for(int ii=0;ii<15;ii++)th1[ii]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			if (isSignal==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
			else {
				//cout<<st<<endl;
				f = TFile::Open(st.data());
			}
		if (!f || !f->IsOpen()){
			//cout<<"error"<<w<<endl;
			continue;
		}
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
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
				if(fabs(thisHiggs->Eta())>2.4)continue;
				FatjetPreSelection.push_back(i);
			}
		
		    //bool isDeltaEta13=0;
			vector<int> FatjetPreSelection1i,FatjetPreSelection1j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
				    if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.3){
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
			*/
			
			nPass[3]++;
			
			//int numHiggsHP=0,numHiggsLP=0;
			
			//bool isPassingTau=0;
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.5)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.75 && FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.5)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.75 && FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				//isPassingTau=1;
			}
			//if(isPassingTau==0)continue;
			if(FatjetPreSelection4i.size()<1)continue;
			nPass[3]++;
			
			
			//for(unsigned int i=0;i<FatjetPreSelection4i.size();i++){
			int numSub=0,numdRSub=0,numdR_com_sub;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][0]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][1]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][0]>0.605)numSub++;	
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][1]>0.605)numSub++;			
			//
			
			TLorentzVector* thisHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4i[0]);
			TLorentzVector* thatHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4j[0]);
			TLorentzVector mjjF=*thisHiggsF+*thatHiggsF;
			if(numSub==0)th1[0]->Fill(mjjF.M());
			if(numSub==1)th1[1]->Fill(mjjF.M());
			if(numSub==2)th1[2]->Fill(mjjF.M());
			if(numSub==3)th1[3]->Fill(mjjF.M());
			if(numSub==4)th1[4]->Fill(mjjF.M());
			
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	
	for(int i=0 ;i< 5;i++){
			//th2f[i]->Scale(3000*Xsec/total);
			if (isSignal==0)th1[i]->Scale(2700*Xsec/total);
			//th2[i]->Scale(3000*Xsec/total);
		}
		TFile* outFile = new TFile(Form("root_files_cat_newNtuple/%s.root",st2.data()),"recreate");
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Write();
			th1[i]->Write();
			//th2[i]->Write();
		}
		outFile->Close();
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}

TH1D* HHbbbbAnalyzerCompare(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		int nPass[10]={0};
		int total=0;
		
		TH1D * th1[30];
		th1[0]=new TH1D("cat0","cat0",4000,1000,5000);
		th1[1]=new TH1D("cat1","cat1",4000,1000,5000);
		th1[2]=new TH1D("cat2","cat2",4000,1000,5000);
		th1[3]=new TH1D("cat3","cat3",4000,1000,5000);
		th1[4]=new TH1D("cat4","cat4",4000,1000,5000);
		th1[5]=new TH1D("dRcat0","dRcat0",4000,1000,5000);
		th1[6]=new TH1D("dRcat1","dRcat1",4000,1000,5000);
		th1[7]=new TH1D("dRcat2","dRcat0",4000,1000,5000);
		th1[8]=new TH1D("dRcat3","dRcat1",4000,1000,5000);
		th1[9]=new TH1D("dRcat4","dRcat1",4000,1000,5000);
		th1[10]=new TH1D("dR_com_cat0","dR_com_cat0",4000,1000,5000);
		th1[11]=new TH1D("dR_com_cat1","dR_com_cat1",4000,1000,5000);
		th1[12]=new TH1D("dR_com_cat2","dR_com_cat0",4000,1000,5000);
		th1[13]=new TH1D("dR_com_cat3","dR_com_cat1",4000,1000,5000);
		th1[14]=new TH1D("dR_com_cat4","dR_com_cat1",4000,1000,5000);
		
		for(int ii=0;ii<15;ii++)th1[ii]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			if (isSignal==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
			else {
				//cout<<st<<endl;
				f = TFile::Open(st.data());
			}
		if (!f || !f->IsOpen()){
			//cout<<"error"<<w<<endl;
			continue;
		}
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
		cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	

		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			cout<<jEntry<<endl;
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*  FATjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp"); 
			float*  FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown"); 
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				// to compare with Henry's mode +1
				TLorentzVector new_vector_up = (*thisHiggs)*(1+FATjetCorrUncUp[i] );
				TLorentzVector Henry_vector_up = jetEnergyScale(1,thisHiggs->Pt(),thisHiggs->Eta(),*thisHiggs);
				
				if(new_vector_up.Px()!=Henry_vector_up.Px())cout<<new_vector_up.Px()<<","<<Henry_vector_up.Px()<<endl;
				if(new_vector_up.Py()!=Henry_vector_up.Py())cout<<new_vector_up.Py()<<","<<Henry_vector_up.Py()<<endl;
				if(new_vector_up.Pz()!=Henry_vector_up.Pz())cout<<new_vector_up.Pz()<<","<<Henry_vector_up.Pz()<<endl;
				if(new_vector_up.E()!=Henry_vector_up.E())cout<<new_vector_up.E()<<","<<Henry_vector_up.E()<<endl;
				if(new_vector_up.M()!=Henry_vector_up.M())cout<<new_vector_up.M()<<","<<Henry_vector_up.M()<<endl;
				// to compare with Henry's mode -1
				TLorentzVector new_vector_down = (*thisHiggs)*(1-FATjetCorrUncDown[i] );
				TLorentzVector Henry_vector_down =  jetEnergyScale(-1,thisHiggs->Pt(),thisHiggs->Eta(),*thisHiggs);
				
				if(new_vector_down.Px()!=Henry_vector_down.Px())cout<<new_vector_down.Px()<<","<<Henry_vector_down.Px()<<endl;
				if(new_vector_down.Py()!=Henry_vector_down.Py())cout<<new_vector_down.Py()<<","<<Henry_vector_down.Py()<<endl;
				if(new_vector_down.Pz()!=Henry_vector_down.Pz())cout<<new_vector_down.Pz()<<","<<Henry_vector_down.Pz()<<endl;
				if(new_vector_down.E()!=Henry_vector_down.E())cout<<new_vector_down.E()<<","<<Henry_vector_down.E()<<endl;
				if(new_vector_down.M()!=Henry_vector_down.M())cout<<new_vector_down.M()<<","<<Henry_vector_down.M()<<endl;
			}
		
		   
			
	
			
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	
	for(int i=0 ;i< 5;i++){
			//th2f[i]->Scale(3000*Xsec/total);
			if (isSignal==0)th1[i]->Scale(2700*Xsec/total);
			//th2[i]->Scale(3000*Xsec/total);
		}
		
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}

TH1D* HHbbbbAnalyzerUp(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		int nPass[10]={0};
		int total=0;
		
		TH1D * th1[30];
		th1[0]=new TH1D("cat0","cat0",4000,1000,5000);
		th1[1]=new TH1D("cat1","cat1",4000,1000,5000);
		th1[2]=new TH1D("cat2","cat2",4000,1000,5000);
		th1[3]=new TH1D("cat3","cat3",4000,1000,5000);
		th1[4]=new TH1D("cat4","cat4",4000,1000,5000);
		th1[5]=new TH1D("dRcat0","dRcat0",4000,1000,5000);
		th1[6]=new TH1D("dRcat1","dRcat1",4000,1000,5000);
		th1[7]=new TH1D("dRcat2","dRcat0",4000,1000,5000);
		th1[8]=new TH1D("dRcat3","dRcat1",4000,1000,5000);
		th1[9]=new TH1D("dRcat4","dRcat1",4000,1000,5000);
		th1[10]=new TH1D("dR_com_cat0","dR_com_cat0",4000,1000,5000);
		th1[11]=new TH1D("dR_com_cat1","dR_com_cat1",4000,1000,5000);
		th1[12]=new TH1D("dR_com_cat2","dR_com_cat0",4000,1000,5000);
		th1[13]=new TH1D("dR_com_cat3","dR_com_cat1",4000,1000,5000);
		th1[14]=new TH1D("dR_com_cat4","dR_com_cat1",4000,1000,5000);
		
		for(int ii=0;ii<15;ii++)th1[ii]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			if (isSignal==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
			else {
				//cout<<st<<endl;
				f = TFile::Open(st.data());
			}
		if (!f || !f->IsOpen()){
			//cout<<"error"<<w<<endl;
			continue;
		}
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
		cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	

		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*  FATjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp"); 
			
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				TLorentzVector new_vector_up = (*thisHiggs)*(1+FATjetCorrUncUp[i] );
				if(new_vector_up.Pt()<200)continue;
				if(fabs(new_vector_up.Eta())>2.4)continue;
				FatjetPreSelection.push_back(i);
			}
		
		    //bool isDeltaEta13=0;
			vector<int> FatjetPreSelection1i,FatjetPreSelection1j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				TLorentzVector thisHiggs_up = (*thisHiggs)*(1+FATjetCorrUncUp[i] );
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
					TLorentzVector thatHiggs_up = (*thatHiggs)*(1+FATjetCorrUncUp[j] );
				    if(fabs(thisHiggs_up.Eta()-thatHiggs_up.Eta())<1.3){
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
				TLorentzVector thisHiggs_up = (*thisHiggs)*(1+FATjetCorrUncUp[FatjetPreSelection1i[i]] );
				TLorentzVector thatHiggs_up = (*thatHiggs)*(1+FATjetCorrUncUp[FatjetPreSelection1j[i]] );
				TLorentzVector mjj=thisHiggs_up+thatHiggs_up;
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
			*/
			
			nPass[3]++;
			
			//int numHiggsHP=0,numHiggsLP=0;
			
			//bool isPassingTau=0;
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.5)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.75 && FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.5)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.75 && FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				//isPassingTau=1;
			}
			//if(isPassingTau==0)continue;
			if(FatjetPreSelection4i.size()<1)continue;
			nPass[3]++;
			
			
			//for(unsigned int i=0;i<FatjetPreSelection4i.size();i++){
			int numSub=0,numdRSub=0,numdR_com_sub;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][0]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][1]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][0]>0.605)numSub++;	
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][1]>0.605)numSub++;			
			//
			
			TLorentzVector* thisHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4i[0]);
			TLorentzVector* thatHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4j[0]);
			TLorentzVector thisHiggsF_up = (*thisHiggsF)*(1+FATjetCorrUncUp[FatjetPreSelection4i[0]] );
			TLorentzVector thatHiggsF_up = (*thatHiggsF)*(1+FATjetCorrUncUp[FatjetPreSelection4j[0]] );
			TLorentzVector mjjF=thisHiggsF_up+thatHiggsF_up;
			if(numSub==0)th1[0]->Fill(mjjF.M());
			if(numSub==1)th1[1]->Fill(mjjF.M());
			if(numSub==2)th1[2]->Fill(mjjF.M());
			if(numSub==3)th1[3]->Fill(mjjF.M());
			if(numSub==4)th1[4]->Fill(mjjF.M());
			
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	
	for(int i=0 ;i< 5;i++){
			//th2f[i]->Scale(3000*Xsec/total);
			if (isSignal==0)th1[i]->Scale(2700*Xsec/total);
			//th2[i]->Scale(3000*Xsec/total);
		}
		TFile* outFile = new TFile(Form("root_files_cat_UpNtuple/%s.root",st2.data()),"recreate");
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Write();
			th1[i]->Write();
			//th2[i]->Write();
		}
		outFile->Close();
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}

TH1D* HHbbbbAnalyzerDown(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		int nPass[10]={0};
		int total=0;
		
		TH1D * th1[30];
		th1[0]=new TH1D("cat0","cat0",4000,1000,5000);
		th1[1]=new TH1D("cat1","cat1",4000,1000,5000);
		th1[2]=new TH1D("cat2","cat2",4000,1000,5000);
		th1[3]=new TH1D("cat3","cat3",4000,1000,5000);
		th1[4]=new TH1D("cat4","cat4",4000,1000,5000);
		th1[5]=new TH1D("dRcat0","dRcat0",4000,1000,5000);
		th1[6]=new TH1D("dRcat1","dRcat1",4000,1000,5000);
		th1[7]=new TH1D("dRcat2","dRcat0",4000,1000,5000);
		th1[8]=new TH1D("dRcat3","dRcat1",4000,1000,5000);
		th1[9]=new TH1D("dRcat4","dRcat1",4000,1000,5000);
		th1[10]=new TH1D("dR_com_cat0","dR_com_cat0",4000,1000,5000);
		th1[11]=new TH1D("dR_com_cat1","dR_com_cat1",4000,1000,5000);
		th1[12]=new TH1D("dR_com_cat2","dR_com_cat0",4000,1000,5000);
		th1[13]=new TH1D("dR_com_cat3","dR_com_cat1",4000,1000,5000);
		th1[14]=new TH1D("dR_com_cat4","dR_com_cat1",4000,1000,5000);
		
		for(int ii=0;ii<15;ii++)th1[ii]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			if (isSignal==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
			else {
				//cout<<st<<endl;
				f = TFile::Open(st.data());
			}
		if (!f || !f->IsOpen()){
			//cout<<"error"<<w<<endl;
			continue;
		}
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
		cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	

		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*  FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown"); 
			
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				TLorentzVector new_vector_up = (*thisHiggs)*(1-FATjetCorrUncDown[i] );
				if(new_vector_up.Pt()<200)continue;
				if(fabs(new_vector_up.Eta())>2.4)continue;
				FatjetPreSelection.push_back(i);
			}
		
		    //bool isDeltaEta13=0;
			vector<int> FatjetPreSelection1i,FatjetPreSelection1j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				TLorentzVector thisHiggs_up = (*thisHiggs)*(1-FATjetCorrUncDown[i] );
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
					TLorentzVector thatHiggs_up = (*thatHiggs)*(1-FATjetCorrUncDown[j] );
				    if(fabs(thisHiggs_up.Eta()-thatHiggs_up.Eta())<1.3){
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
				TLorentzVector thisHiggs_up = (*thisHiggs)*(1-FATjetCorrUncDown[FatjetPreSelection1i[i]] );
				TLorentzVector thatHiggs_up = (*thatHiggs)*(1-FATjetCorrUncDown[FatjetPreSelection1j[i]] );
				TLorentzVector mjj=thisHiggs_up+thatHiggs_up;
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
			*/
			
			nPass[3]++;
			
			//int numHiggsHP=0,numHiggsLP=0;
			
			//bool isPassingTau=0;
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.5)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.75 && FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.5)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.75 && FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				//isPassingTau=1;
			}
			//if(isPassingTau==0)continue;
			if(FatjetPreSelection4i.size()<1)continue;
			nPass[3]++;
			
			
			//for(unsigned int i=0;i<FatjetPreSelection4i.size();i++){
			int numSub=0,numdRSub=0,numdR_com_sub;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][0]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][1]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][0]>0.605)numSub++;	
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][1]>0.605)numSub++;			
			//
			
			TLorentzVector* thisHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4i[0]);
			TLorentzVector* thatHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4j[0]);
			TLorentzVector thisHiggsF_up = (*thisHiggsF)*(1-FATjetCorrUncDown[FatjetPreSelection4i[0]] );
			TLorentzVector thatHiggsF_up = (*thatHiggsF)*(1-FATjetCorrUncDown[FatjetPreSelection4j[0]] );
			TLorentzVector mjjF=thisHiggsF_up+thatHiggsF_up;
			if(numSub==0)th1[0]->Fill(mjjF.M());
			if(numSub==1)th1[1]->Fill(mjjF.M());
			if(numSub==2)th1[2]->Fill(mjjF.M());
			if(numSub==3)th1[3]->Fill(mjjF.M());
			if(numSub==4)th1[4]->Fill(mjjF.M());
			
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	
	for(int i=0 ;i< 5;i++){
			//th2f[i]->Scale(3000*Xsec/total);
			if (isSignal==0)th1[i]->Scale(2700*Xsec/total);
			//th2[i]->Scale(3000*Xsec/total);
		}
		TFile* outFile = new TFile(Form("root_files_cat_DownNtuple/%s.root",st2.data()),"recreate");
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Write();
			th1[i]->Write();
			//th2[i]->Write();
		}
		outFile->Close();
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}

TH1D* HHbbbbAnalyzerCSVEff(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		db1=0;db2=0;db3=0;db4=0;db5=0;
		
		int nPass[10]={0};
		int total=0;
		
		TH1D * th1[30];
		th1[0]=new TH1D("cat0","cat0",4000,1000,5000);
		th1[1]=new TH1D("cat1","cat1",4000,1000,5000);
		th1[2]=new TH1D("cat2","cat2",4000,1000,5000);
		th1[3]=new TH1D("cat3","cat3",4000,1000,5000);
		th1[4]=new TH1D("cat4","cat4",4000,1000,5000);
		th1[5]=new TH1D("dRcat0","dRcat0",4000,1000,5000);
		th1[6]=new TH1D("dRcat1","dRcat1",4000,1000,5000);
		th1[7]=new TH1D("dRcat2","dRcat0",4000,1000,5000);
		th1[8]=new TH1D("dRcat3","dRcat1",4000,1000,5000);
		th1[9]=new TH1D("dRcat4","dRcat1",4000,1000,5000);
		th1[10]=new TH1D("dR_com_cat0","dR_com_cat0",4000,1000,5000);
		th1[11]=new TH1D("dR_com_cat1","dR_com_cat1",4000,1000,5000);
		th1[12]=new TH1D("dR_com_cat2","dR_com_cat0",4000,1000,5000);
		th1[13]=new TH1D("dR_com_cat3","dR_com_cat1",4000,1000,5000);
		th1[14]=new TH1D("dR_com_cat4","dR_com_cat1",4000,1000,5000);
		
		for(int ii=0;ii<15;ii++)th1[ii]->Sumw2();
		
		for (int w=1;w<wM;w++){
			//cout<<Form("%s%d.root",st.data(),w)<<endl;
		//if(SignalBkgNum==1)f = TFile::Open(Form("%s",st.data()));
		//else 
			if (isSignal==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
			else {
				//cout<<st<<endl;
				f = TFile::Open(st.data());
			}
		if (!f || !f->IsOpen()){
			//cout<<"error"<<w<<endl;
			continue;
		}
		TDirectory * dir;
		//if(SignalBkgNum==1)dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else 
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
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
				if(fabs(thisHiggs->Eta())>2.4)continue;
				FatjetPreSelection.push_back(i);
			}
		
		    //bool isDeltaEta13=0;
			vector<int> FatjetPreSelection1i,FatjetPreSelection1j;
			for(unsigned int i=0;i<FatjetPreSelection.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[i]);
				for(unsigned int j=i+1;j<FatjetPreSelection.size();j++){
					TLorentzVector* thatHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection[j]);
				    if(fabs(thisHiggs->Eta()-thatHiggs->Eta())<1.3){
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
			*/
			
			nPass[3]++;
			
			//int numHiggsHP=0,numHiggsLP=0;
			
			//bool isPassingTau=0;
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.5)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.75 && FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.5)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.75 && FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				//isPassingTau=1;
			}
			//if(isPassingTau==0)continue;
			if(FatjetPreSelection4i.size()<1)continue;
			nPass[3]++;
			
			
			//for(unsigned int i=0;i<FatjetPreSelection4i.size();i++){
			int numSub=0,numdRSub=0,numdR_com_sub;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][0]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4i[0]][1]>0.605)numSub++;
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][0]>0.605)numSub++;	
			if(FATsubjetSDCSV[FatjetPreSelection4j[0]][1]>0.605)numSub++;			
			//
			
			TLorentzVector* thisHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4i[0]);
			TLorentzVector* thatHiggsF = (TLorentzVector*)FATjetP4->At(FatjetPreSelection4j[0]);
			TLorentzVector mjjF=*thisHiggsF+*thatHiggsF;
			if(numSub==0)db1++;
			if(numSub==1)db2++;
			if(numSub==2)db3++;
			if(numSub==3)db4++;
			if(numSub==4)db5++;
			
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	
	cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	
	
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}

void HHbbbbAnalyzer(int a){
	//setNCUStyle();
	
	string st1[20]={
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215300/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215340/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215429/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215514/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215644/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
	
	/*6*///"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215734/0000/NCUGlobalTuples_",
	    //"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215815/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215903/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_215944/0000/NCUGlobalTuples_",
		//"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/BulkGravToHHBBBBSignal_eleIDjet_CMSSW7412_20151006/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph/crab_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraphMC25ns_eleIDjet_CMSSW7412_20151006/151007_220034/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
	/*11*/
		//"/data2/syu/13TeV/BulkGravTohhTohbbhbb/softdrop_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220208/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220428/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT500to7000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		
	/*16*/	
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_",
	
	};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	string  masspointN[19]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500","QCD100","QCD200","QCD300","QCD500","QCD700","QCD1000","QCD1500","QCD2000"};
	double eff1,eff2,eff3,eff4,eff5;
	TH1D* th1[20];
	int aa[20]={
		2,2,2,2,2,
		2,2,2,2,2,
		2,812,516,565,496,
		394,155,103,71
	};
	double Xsec[20]={
		1,1,1,1,1,
		1,1,1,1,1,
		1,27850000,1717000,351300,31630,
		6802,1206,120.4,25.24
	};
	bool sig=0;
	if (a<11)sig=1;
	//th1[a]=HHbbbbAnalyzerBase(aa[a],st1[a],masspointN[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerCompare(aa[a],st1[a],masspointN[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	th1[a]=HHbbbbAnalyzerCSVEff(aa[a],st1[a],masspointN[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	
	TFile * tf ;
	TH1D* th2[5];
	
	if(a==0){
		tf=new TFile("Analyzer1228.root","recreate");
		th2[0]=new TH1D("th1","",19,0,19);
		th2[1]=new TH1D("th2","",19,0,19);
		th2[2]=new TH1D("th3","",19,0,19);
		th2[3]=new TH1D("th4","",19,0,19);
		th2[4]=new TH1D("th5","",19,0,19);
		
		th2[0]->SetBinContent(a+1,eff1);
		th2[1]->SetBinContent(a+1,eff2);
		th2[2]->SetBinContent(a+1,eff3);
		th2[3]->SetBinContent(a+1,eff4);
		th2[4]->SetBinContent(a+1,eff5);
		
		for(int j=1;j<12;j++)th2[0]->GetXaxis()->SetBinLabel(j,Form("%s",masspoint[j-1].data()));
		th2[0]->GetXaxis()->SetBinLabel(12,"QCDHT100-200");
	th2[0]->GetXaxis()->SetBinLabel(13,"QCDHT200-300");
	th2[0]->GetXaxis()->SetBinLabel(14,"QCDHT300-400");
	th2[0]->GetXaxis()->SetBinLabel(15,"QCDHT500-700");
	th2[0]->GetXaxis()->SetBinLabel(16,"QCDHT700-1000");
	th2[0]->GetXaxis()->SetBinLabel(17,"QCDHT1000-1500");
	th2[0]->GetXaxis()->SetBinLabel(18,"QCDHT1500-2000");
	th2[0]->GetXaxis()->SetBinLabel(19,"QCDHT2000-Inf.");
		
		for (int i=0;i<5;i++)th2[i]->Write();
		tf->Close();
		
	}
    else {
		
		tf=new TFile("Analyzer1228.root","update");
		
		th2[0]=(TH1D*) tf->FindObjectAny("th1");
		th2[1]=(TH1D*) tf->FindObjectAny("th2");
		th2[2]=(TH1D*) tf->FindObjectAny("th3");
		th2[3]=(TH1D*) tf->FindObjectAny("th4");
		th2[4]=(TH1D*) tf->FindObjectAny("th5");
		
		th2[0]->SetBinContent(a+1,eff1);
		th2[1]->SetBinContent(a+1,eff2);
		th2[2]->SetBinContent(a+1,eff3);
		th2[3]->SetBinContent(a+1,eff4);
		th2[4]->SetBinContent(a+1,eff5);
		
		
		
		
		for (int i=0;i<5;i++)th2[i]->Write();
		tf->Close();
	}
	
	
	/*
	
	
	c1 = new TCanvas("c1","",1360,768);
	
	setNCUStyle(true);
	TLegend* leg ;
	leg=new TLegend(0.591452,0.562447,0.690645,0.783966);
	//leg->SetFillColor(18);
	//leg->SetFillStyle(0);
	//leg->SetTextSize(0.02);
	//leg->SetBorderSize(2);
	
	
	
	
	for(int i=0;i<5;i++){
		th2[i]=new TH1D("","",18,0,18);
		for(int j=1;j<19;j++){
			if(i==0)th2[i]->SetBinContent(j,eff1[j-1]);
			if(i==1)th2[i]->SetBinContent(j,eff2[j-1]);
			if(i==2)th2[i]->SetBinContent(j,eff3[j-1]);
			if(i==3)th2[i]->SetBinContent(j,eff4[j-1]);
			if(i==4)th2[i]->SetBinContent(j,eff5[j-1]);
			if(i<11)th2[i]->GetXaxis()->SetBinLabel(j,Form("%s",masspoint[j-1].data()));
		}
		
	}
	//TLegend* leg ;
	//leg=new TLegend(0.691452,0.662447,0.890645,0.883966);
	
	th2[0]->GetXaxis()->SetBinLabel(11,"QCDHT100-200");
	th2[0]->GetXaxis()->SetBinLabel(12,"QCDHT200-300");
	th2[0]->GetXaxis()->SetBinLabel(13,"QCDHT300-400");
	th2[0]->GetXaxis()->SetBinLabel(14,"QCDHT500-700");
	th2[0]->GetXaxis()->SetBinLabel(15,"QCDHT700-1000");
	th2[0]->GetXaxis()->SetBinLabel(16,"QCDHT1000-1500");
	th2[0]->GetXaxis()->SetBinLabel(17,"QCDHT1500-2000");
	th2[0]->GetXaxis()->SetBinLabel(18,"QCDHT2000-Inf.");
	
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
	c1->Print("pdf/cutFlow1110.png");
    //c1->SetLogy(1);
	
	//c1->Print("pdf/cutM7.png");
	
	
	
	
	
	ofstream myfile;
	myfile.open ("1110cout.txt");
	
	for(int i=0;i<18;i++)eff5[i]*=1000;
	eff5[10]*=27850000;
	eff5[11]*=1717000;
	eff5[12]*=351300;
	eff5[13]*=31630 ;
	eff5[14]*=6802;
	eff5[15]*=1206;
	eff5[16]*=120.4;
	eff5[17]*=25.24 ;
	
	
	
	myfile<<"DATA 0"<<endl;
	double bkgtemp=0;
	for(int i=10;i<18;i++)bkgtemp+= eff5[i];
	myfile<<"QCD "<<bkgtemp<<endl;
	for(int i=0;i<10;i++)myfile<<"M"<<masspoint[i]<<" "<<eff5[i]<<endl;
	myfile.close();
	
	*/
}