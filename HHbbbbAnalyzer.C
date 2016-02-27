#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
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
#include <algorithm>
//#include "BTagCalibrationStandalone.h"

TCanvas* c1;

//old cut flow
/*
TH1D* HHbbbbAnalyzerBase(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		int nPass[20]={0};
		int total=0;
		
		TH1D * th1[30],*thd;
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
		thd=new TH1D("dEta","dEta",40,-2,2);
		
		TH1D* thht=new TH1D("HT","HT",1000,0,3000);
		TH1D* thht2=new TH1D("HTcut","HT",1000,0,3000);
		
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
		
		if(w%10==0)cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	
		//data.Print();
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			Float_t HT=data.GetFloat("HT"); 
			
			thht->Fill(HT);
			
			if(nJet<1)continue;
			
			vector<int> FatjetPreSelection_00;
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<200)continue;
				FatjetPreSelection_00.push_back(i);
			}
			if(FatjetPreSelection_00.size()<1)continue;
			nPass[0]++;
			vector<int> FatjetPreSelection;
			for(unsigned int i=0;i<FatjetPreSelection_00.size();i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(FatjetPreSelection_00[i]);
				if(fabs(thisHiggs->Eta())>2.4)continue;
				FatjetPreSelection.push_back(i);
			}
			if(FatjetPreSelection.size()<1)continue;
			nPass[1]++;
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
					thd->Fill(thisHiggs->Eta()-thatHiggs->Eta());
				}
				//if(isDeltaEta13)break;
			}
			if(FatjetPreSelection1i.size()<1)continue;
	
			nPass[2]++;
	
			
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
			nPass[3]++;
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
			nPass[4]++;
			//cout<<jEntry<<"2"<<endl;
			
			
			
			
			//bool isPassingTau=0;
			vector<int> FatjetPreSelection4i,FatjetPreSelection4j;
			vector<int> FatjetPreSelection5i,FatjetPreSelection5j;
			for(unsigned int i=0;i<FatjetPreSelection3i.size();i++){
				bool thisHiggsHP=0,thisHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.6)thisHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]<0.75 && FATjetTau2[FatjetPreSelection3i[i]]/FATjetTau1[FatjetPreSelection3i[i]]>0.5)thisHiggsLP=1;
				if(thisHiggsHP)thisHiggsLP=1;
				bool thatHiggsHP=0,thatHiggsLP=0;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.6)thatHiggsHP=1;
				if(FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]<0.75 && FATjetTau2[FatjetPreSelection3j[i]]/FATjetTau1[FatjetPreSelection3j[i]]>0.5)thatHiggsLP=1;
				if(thatHiggsHP)thatHiggsLP=1;
				
				if(thisHiggsLP==0 || thatHiggsLP==0 )continue;
				if(thisHiggsHP==0 && thatHiggsHP==0) continue;
				
				FatjetPreSelection4i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection4j.push_back(FatjetPreSelection3j[i]);
				//isPassingTau=1;
				if(thisHiggsHP==0 || thatHiggsHP==0) continue;
				FatjetPreSelection5i.push_back(FatjetPreSelection3i[i]);
				FatjetPreSelection5j.push_back(FatjetPreSelection3j[i]);
			}
			//if(isPassingTau==0)continue;
			if(FatjetPreSelection4i.size()<1)continue;
			nPass[5]++;
			
			thht2->Fill(HT);
			
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
			if(numSub==0){
				th1[0]->Fill(mjjF.M());
				nPass[6]++;
			}
			if(numSub==1){
				th1[1]->Fill(mjjF.M());
				nPass[7]++;
			}
			if(numSub==2){
				th1[2]->Fill(mjjF.M());
				nPass[8]++;
			}
			if(numSub==3){
				th1[3]->Fill(mjjF.M());
				nPass[9]++;
			}
			if(numSub==4){
				th1[4]->Fill(mjjF.M());
				nPass[10]++;
			}
			
			
			if(FatjetPreSelection5i.size()<1)continue;
			int numSubHPHP=0;
			if(FATsubjetSDCSV[FatjetPreSelection5i[0]][0]>0.605)numSubHPHP++;
			if(FATsubjetSDCSV[FatjetPreSelection5i[0]][1]>0.605)numSubHPHP++;
			if(FATsubjetSDCSV[FatjetPreSelection5j[0]][0]>0.605)numSubHPHP++;	
			if(FATsubjetSDCSV[FatjetPreSelection5j[0]][1]>0.605)numSubHPHP++;	
			if(numSubHPHP==3)nPass[11]++;
			
			}
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			//cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			th2->SetBinContent(i+1,double(nPass[i])/total);
			th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	//return  double(nPass[4])/total;
	db1=double(nPass[0])/total;
	db2=double(nPass[1])/total;
	db3=double(nPass[2])/total;
	db4=double(nPass[3])/total;
	db5=double(nPass[4])/total; 
	//cout<<db1<<","<<db2<<","<<db3<<","<<db4<<","<<db5<<",p="<<nPass[4]<<endl;
	
	cout<<"nPass[0]="<<nPass[0]<<endl;
	cout<<"nPass[1]="<<nPass[1]<<endl;
	cout<<"nPass[2]="<<nPass[2]<<endl;
	cout<<"nPass[3]="<<nPass[3]<<endl;
	cout<<"nPass[4]="<<nPass[4]<<endl;
	cout<<"nPass[5]="<<nPass[5]<<endl;
	cout<<"nPass[6]="<<nPass[6]<<endl;
	cout<<"nPass[7]="<<nPass[7]<<endl;
	cout<<"nPass[8]="<<nPass[8]<<endl;
	cout<<"nPass[9]="<<nPass[9]<<endl;
	cout<<"nPass[10]="<<nPass[10]<<endl;
	cout<<"nPass[11]="<<nPass[11]<<endl;
	
	TH1D* cutflow=new TH1D("cutflow","",12,0,12);
	TH1D* cutflow2=new TH1D("cutflow2","",12,0,12);
	
	for(int i=0;i<12;i++){
		
		cutflow->SetBinContent(i+1,nPass[i]*2197.7*Xsec/total);
		cutflow2->SetBinContent(i+1,nPass[i]);
	}
	
	
	
	for(int i=0 ;i< 5;i++){
			//th2f[i]->Scale(3000*Xsec/total);
		if (isSignal==0){
				//cout<<Xsec<<","<<total<<endl;
				th1[i]->Scale(2197.7*Xsec/total);
				
		}
			//th2[i]->Scale(3000*Xsec/total);
	}
		
		if (isSignal==0){
			thd->Scale(2197.7*Xsec/total);
			thht->Scale(2197.7*Xsec/total);
			thht2->Scale(2197.7*Xsec/total);	
		}
		
		
		
		TFile* outFile = new TFile(Form("root_files_cat_newNtuple/%s.root",st2.data()),"recreate");
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Write();
			th1[i]->Write();
			//th2[i]->Write();
		}
		thd->Write();
		cutflow->Write();
		cutflow2->Write();
		outFile->Close();
		
		TFile* outFile2 = new TFile(Form("root_files_bkg_HT/%s.root",st2.data()),"recreate");
		thht->Write();
		thht2->Write();
		outFile2->Close();
		
		for(int i=0 ;i< 5;i++){
			//th2f[i]->Clear();
			th1[i]->Clear();
			thd->Clear();
			thht->Clear();
			thht2->Clear();
			//th2[i]->Clear();
		}
	
	
	return th2;
	
}
*/


void HHbbbbAnalyzerBaseC(int wMs,int wM, string st,string st2,double Xsec,bool nameRoot=0,bool isData=0){	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	
	TH1D * th1[50];
	
	th1[0]=new TH1D("Pt_j0_sj0_0b","Pt_j0_sj0_0b",200,0,2000);
	th1[1]=new TH1D("Pt_j0_sj1_0b","Pt_j0_sj1_0b",200,0,2000);
	th1[2]=new TH1D("Pt_j1_sj0_0b","Pt_j1_sj0_0b",200,0,2000);
	th1[3]=new TH1D("Pt_j1_sj1_0b","Pt_j1_sj1_0b",200,0,2000);
	
	th1[4]=new TH1D("Pt_j0_sj0_1b","Pt_j0_sj0_1b",200,0,2000);
	th1[5]=new TH1D("Pt_j0_sj1_1b","Pt_j0_sj1_1b",200,0,2000);
	th1[6]=new TH1D("Pt_j1_sj0_1b","Pt_j1_sj0_1b",200,0,2000);
	th1[7]=new TH1D("Pt_j1_sj1_1b","Pt_j1_sj1_1b",200,0,2000);
	
	th1[8]=new TH1D("Pt_j0_sj0_2b","Pt_j0_sj0_2b",200,0,2000);
	th1[9]=new TH1D("Pt_j0_sj1_2b","Pt_j0_sj1_2b",200,0,2000);
	th1[10]=new TH1D("Pt_j1_sj0_2b","Pt_j1_sj0_2b",200,0,2000);
	th1[11]=new TH1D("Pt_j1_sj1_2b","Pt_j1_sj1_2b",200,0,2000);
	
	th1[12]=new TH1D("deltaR0_0b","deltaR0_0b",20,0,1);
	th1[13]=new TH1D("deltaR1_0b","deltaR1_0b",20,0,1);
	th1[14]=new TH1D("deltaR0_1b","deltaR0_1b",20,0,1);
	th1[15]=new TH1D("deltaR1_1b","deltaR1_1b",20,0,1);
	th1[16]=new TH1D("deltaR0_2b","deltaR0_2b",20,0,1);
	th1[17]=new TH1D("deltaR1_2b","deltaR1_2b",20,0,1);
	
	th1[18]=new TH1D("pt0","pt0",480,200,1400);
	th1[19]=new TH1D("pt1","pt1",480,200,1400);
	
	th1[20]=new TH1D("DeltaEta","DeltaEta",16,-0.1,1.5);
	th1[21]=new TH1D("eta_j0","eta_j0",60,-3,3);
	th1[22]=new TH1D("eta_j1","eta_j1",480,200,1400);
	th1[23]=new TH1D("HT12","HT12",80,600,2600);
	th1[24]=new TH1D("prMass_j0","prMass_j0",12,90,150);
	th1[25]=new TH1D("prMassCut_j0","prMassCut_j0",12,90,150);
	th1[26]=new TH1D("prMass_j1","prMass_j1",12,90,150);
	th1[27]=new TH1D("prMassCut_j1","prMassCut_j1",12,90,150);
	th1[28]=new TH1D("tau21_j0","tau21_j0",10,0,1);
	th1[29]=new TH1D("tau21_j1","tau21_j1",10,0,1);
	
	th1[30]=new TH1D("totalMass_0b","totalMass_0b",75,1000,4000);
	th1[31]=new TH1D("totalMass_1b","totalMass_1b",75,1000,4000);
	th1[32]=new TH1D("totalMass_2b","totalMass_2b",75,1000,4000);
	th1[33]=new TH1D("totalMass_3b","totalMass_3b",75,1000,4000);
	th1[34]=new TH1D("totalMass_4b","totalMass_4b",75,1000,4000);
	
	
	std::vector<TString> eventlist;                                                                                                                                            
	if(isData==1){
		std::vector<TString> eventlist;                                                                                                                                            
		std::ifstream fin("somfilename.txt");
		eventlist.clear();
		std::string line;
		if (fin.is_open())  {
			while (fin.good())  {
				getline(fin,line);
				eventlist.push_back(line);
		}
		fin.close();
		}
	}
		
	for (int w=wMs;w<wM;w++){
		
		if (nameRoot==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		
		if (nameRoot==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		
		if(w%20==0)cout<<w<<endl;
		
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			
			//Int_t event        = data.GetInt("eventId");
			//0. has a good vertex
			Int_t nVtx        = data.GetInt("nVtx");
			if(nVtx<1)continue;
			nPass[0]++;
	
			int nFATJet         = data.GetInt("FATnJet");
			const int nJets=nFATJet;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			//Float_t*  fatjetTau4 = data.GetPtrFloat("FATjetTau4");
			Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			//  vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			//vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			if(nJets<2)continue;
			//-----------------------------------------
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(0);
			TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(1);
			//2. Pt and tightId-----------------------------------------
			if(thisJet->Pt()<200)continue;
			if(thatJet->Pt()<200)continue;
			if(FATjetPassIDTight[0]==0)continue;
			if(FATjetPassIDTight[1]==0)continue;
			nPass[1]++;
			//3. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[2]++;
			//4. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[3]++;
			//5. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[4]++;
			//6. fatjetPRmassL2L3Corr-----------------------------------------
			th1[24]->Fill(fatjetPRmassL2L3Corr[0]);
			th1[26]->Fill(fatjetPRmassL2L3Corr[1]);
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[5]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		           tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			th1[28]->Fill(tau21_1);	   
			th1[29]->Fill(tau21_2);	   
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[6]++;
			
			if(isData==1){
				//data:csc2015
				TString thisEvent;                                                   
				Int_t  runId = data.GetInt("runId");
				Int_t  lumiSection = data.GetInt("lumiSection");
				Int_t  eventId = data.GetInt("eventId");
				thisEvent.Form("%d:%d:%d",runId,lumiSection,eventId);                                                                                                                 
				if ( std::find(eventlist.begin(), eventlist.end(), thisEvent) != eventlist.end()){
					std::cout<<" match found to skip event "<<std::endl;
					continue;
				}                                                                              
				dataPassingcsc++;
			}
			
			th1[18]->Fill(thisJet->Pt());
			th1[19]->Fill(thatJet->Pt());
			th1[23]->Fill(thisJet->Pt()+thatJet->Pt());
			th1[20]->Fill(dEta);
			th1[21]->Fill(thisJet->Eta());
			th1[22]->Fill(thatJet->Eta());
			th1[25]->Fill(fatjetPRmassL2L3Corr[0]);
			th1[27]->Fill(fatjetPRmassL2L3Corr[1]);
			
			//8.btag
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			
	
			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			double dr1=thisSub1->DeltaR(*thisSub2),dr2=thatSub1->DeltaR(*thatSub2);
			
			if(nbtag==0){
				nPass[7]++;
				th1[0]->Fill(thisSub1->Pt());
				th1[1]->Fill(thisSub2->Pt());
				th1[2]->Fill(thatSub1->Pt());
				th1[3]->Fill(thatSub2->Pt());
				th1[12]->Fill(dr1);
				th1[13]->Fill(dr2);
				th1[30]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==1){
				nPass[8]++;
				th1[4]->Fill(thisSub1->Pt());
				th1[5]->Fill(thisSub2->Pt());
				th1[6]->Fill(thatSub1->Pt());
				th1[7]->Fill(thatSub2->Pt());
				th1[14]->Fill(dr1);
				th1[15]->Fill(dr2);
				th1[31]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==2){
				nPass[9]++;
				th1[8]->Fill(thisSub1->Pt());
				th1[9]->Fill(thisSub2->Pt());
				th1[10]->Fill(thatSub1->Pt());
				th1[11]->Fill(thatSub2->Pt());
				th1[16]->Fill(dr1);
				th1[17]->Fill(dr2);
				th1[32]->Fill((*thisJet+*thatJet).M());
				
			}
			if(nbtag==3){
				nPass[10]++;
				th1[33]->Fill((*thisJet+*thatJet).M());
			}
			if(nbtag==4){
				nPass[11]++;
				th1[34]->Fill((*thisJet+*thatJet).M());
			}
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(isData==1)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<12;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TH1D * th2=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	th2->SetBinContent(1,nPass[7]);
	th2->SetBinContent(2,nPass[8]);
	th2->SetBinContent(3,nPass[9]);
	th2->SetBinContent(4,nPass[10]);
	th2->SetBinContent(5,nPass[11]);
	
	TH1D * th2s=new TH1D("NbtagjetS","NbtagjetS",5,-0.5,4.5);
	th2s->SetBinContent(1,nPass[7]);
	th2s->SetBinContent(2,nPass[8]);
	th2s->SetBinContent(3,nPass[9]);
	th2s->SetBinContent(4,nPass[10]);
	th2s->SetBinContent(5,nPass[11]);
	
	
	TFile* outFile = new TFile(Form("root_files_btaggedScaleFactor/%s.root",st2.data()),"recreate");
	if(isData==0){
		for(int i=0;i<35;i++){
			th1[i]->Sumw2();
			th1[i]->Scale(2245.87*Xsec/total);
			
		}
	}
	for(int i=0;i<35;i++)th1[i]->Write();
	th2->Write();
	th2s->Scale(2245.87*Xsec/total);
	th2s->Write();
	outFile->Close();
	
}

/*
void HHbbbbAnalyzerBaseB(int wM, string st,string st2,double Xsec,bool nameRoot=0,bool isData=0){	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	
	TH1D * th1[30];
	
	th1[0]=new TH1D("subjet11_0b","subjet11_0b",200,0,2000);
	th1[1]=new TH1D("subjet12_0b","subjet12_0b",200,0,2000);
	th1[2]=new TH1D("subjet21_0b","subjet21_0b",200,0,2000);
	th1[3]=new TH1D("subjet22_0b","subjet22_0b",200,0,2000);
	
	th1[4]=new TH1D("subjet11_1b","subjet11_1b",200,0,2000);
	th1[5]=new TH1D("subjet12_1b","subjet12_1b",200,0,2000);
	th1[6]=new TH1D("subjet21_1b","subjet21_1b",200,0,2000);
	th1[7]=new TH1D("subjet22_1b","subjet22_1b",200,0,2000);
	
	th1[8]=new TH1D("subjet11_2b","subjet11_2b",200,0,2000);
	th1[9]=new TH1D("subjet12_2b","subjet12_2b",200,0,2000);
	th1[10]=new TH1D("subjet21_2b","subjet21_2b",200,0,2000);
	th1[11]=new TH1D("subjet22_2b","subjet22_2b",200,0,2000);
	
	th1[12]=new TH1D("deltaR1_0b","deltaR1_0b",50,0,4);
	th1[13]=new TH1D("deltaR2_0b","deltaR2_0b",50,0,4);
	th1[14]=new TH1D("deltaR1_1b","deltaR1_1b",50,0,4);
	th1[15]=new TH1D("deltaR2_1b","deltaR2_1b",50,0,4);
	th1[16]=new TH1D("deltaR1_2b","deltaR1_2b",50,0,4);
	th1[17]=new TH1D("deltaR2_2b","deltaR2_2b",50,0,4);
	
	th1[18]=new TH1D("nsubJet1_0b","nsubJet1_0b",10,0,10);
	th1[19]=new TH1D("nsubJet2_0b","nsubJet2_0b",10,0,10);
	th1[20]=new TH1D("nsubJet1_1b","nsubJet1_1b",10,0,10);
	th1[21]=new TH1D("nsubJet2_1b","nsubJet2_1b",10,0,10);
	th1[22]=new TH1D("nsubJet1_2b","nsubJet1_2b",10,0,10);
	th1[23]=new TH1D("nsubJet2_2b","nsubJet2_2b",10,0,10);
	
	BTagCalibration calib("csvv2", "CSVv2.csv");
	BTagCalibrationReader reader(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "comb",               // measurement type
                             "central");           // systematics type
	BTagCalibrationReader reader_up(&calib, BTagEntry::OP_LOOSE, "comb", "up");  // sys up
	BTagCalibrationReader reader_do(&calib, BTagEntry::OP_LOOSE, "comb", "down");  // sys down

	
		
	for (int w=1;w<wM;w++){
		
		if (nameRoot==0)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		
		if (nameRoot==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		
		if(w%20==0)cout<<w<<endl;
		
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			
			//Int_t event        = data.GetInt("eventId");
			//0. has a good vertex
			Int_t nVtx        = data.GetInt("nVtx");
			if(nVtx<1)continue;
			nPass[0]++;
	
			int nFATJet         = data.GetInt("FATnJet");
			const int nJets=nFATJet;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			//Float_t*  fatjetTau4 = data.GetPtrFloat("FATjetTau4");
			Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			//  vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			//vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			if(nJets<2)continue;
			//-----------------------------------------
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(0);
			TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(1);
			//2. Pt and tightId-----------------------------------------
			if(thisJet->Pt()<200)continue;
			if(thatJet->Pt()<200)continue;
			if(FATjetPassIDTight[0]==0)continue;
			if(FATjetPassIDTight[1]==0)continue;
			nPass[1]++;
			//3. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[2]++;
			//4. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[3]++;
			//5. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[4]++;
			//6. fatjetPRmassL2L3Corr-----------------------------------------
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[5]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		   tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[6]++;
			
		
			//8.btag
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			
	
			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			double dr1=thisSub1->DeltaR(*thisSub2),dr2=thatSub1->DeltaR(*thatSub2);
			
			double f1 = reader.eval(BTagEntry::FLAV_B, thisSub1->Eta(), thisSub1->Pt()); 
			double f2 = reader.eval(BTagEntry::FLAV_B, thisSub2->Eta(), thisSub2->Pt()); 
			double f3 = reader.eval(BTagEntry::FLAV_B, thatSub1->Eta(), thatSub1->Pt()); 
			double f4 = reader.eval(BTagEntry::FLAV_B, thatSub2->Eta(), thatSub2->Pt()); 
			
			thisSub1->SetPxPyPzE(subjetSDPx[0][0]*f1,subjetSDPy[0][0]*f1,subjetSDPz[0][0]*f1,subjetSDE[0][0]*f1);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1]*f2,subjetSDPy[0][1]*f2,subjetSDPz[0][1]*f2,subjetSDE[0][1]*f2);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0]*f3,subjetSDPy[1][0]*f3,subjetSDPz[1][0]*f3,subjetSDE[1][0]*f3);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1]*f4,subjetSDPy[1][1]*f4,subjetSDPz[1][1]*f4,subjetSDE[1][1]*f4);
			
			if(nbtag==0){
				nPass[7]++;
				th1[0]->Fill(thisSub1->Pt());
				th1[1]->Fill(thisSub2->Pt());
				th1[2]->Fill(thatSub1->Pt());
				th1[3]->Fill(thatSub2->Pt());
				th1[12]->Fill(dr1);
				th1[13]->Fill(dr2);
				th1[18]->Fill(nSubSoftDropJet[0]);
				th1[19]->Fill(nSubSoftDropJet[1]);
			}
			if(nbtag==1){
				nPass[8]++;
				th1[4]->Fill(thisSub1->Pt());
				th1[5]->Fill(thisSub2->Pt());
				th1[6]->Fill(thatSub1->Pt());
				th1[7]->Fill(thatSub2->Pt());
				th1[14]->Fill(dr1);
				th1[15]->Fill(dr2);
				th1[20]->Fill(nSubSoftDropJet[0]);
				th1[21]->Fill(nSubSoftDropJet[1]);
			}
			if(nbtag==2){
				nPass[9]++;
				th1[8]->Fill(thisSub1->Pt());
				th1[9]->Fill(thisSub2->Pt());
				th1[10]->Fill(thatSub1->Pt());
				th1[11]->Fill(thatSub2->Pt());
				th1[16]->Fill(dr1);
				th1[17]->Fill(dr2);
				th1[22]->Fill(nSubSoftDropJet[0]);
				th1[23]->Fill(nSubSoftDropJet[1]);
			}
			if(nbtag==3)nPass[10]++;
			if(nbtag==4)nPass[11]++;
			
			
			
			
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(isData==1)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<12;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	
	
	TFile* outFile = new TFile(Form("root_files_btaggedScaleFactorB/%s.root",st2.data()),"recreate");
	if(isData==0){
		for(int i=0;i<24;i++){
			th1[i]->Sumw2();
			th1[i]->Scale(2197.7*Xsec/total);
			th1[i]->Write();
		}
	}
	outFile->Close();
	
}
*/

void dataPrinter(int wM, string st,string st2,double Xsec,bool isSignal=0){
	
	//for (int massP=0;massP<1;massP++){
		//cout<<isSignal<<endl;
		TFile *f;
		TTree *tree;
		
		int nPass[20]={0};
		int total=0;
		
				
		for (int w=1;w<2;w++){
		
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
		
			if (isSignal==0)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
			else  {
				//cout<<Form("%s:/tree",st.data())<<endl;
				dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
			}
		
		if(w%10==0)cout<<w<<endl;
		
		
		dir->GetObject("treeMaker",tree);
		
		TreeReader data(tree);
		
		total+=data.GetEntriesFast();
	
		data.Print();
		
		
		for(Long64_t jEntry=0; jEntry<1 ;jEntry++){

  //  if (jEntry % 50000 == 0)
     // fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    //nTotal++;
    Int_t event        = data.GetInt("eventId");
    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++;



    //2. pass electron or muon trigger
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];
	cout<<it<<"="<<thisTrig<<endl;
	// std::cout << thisTrig << " : " << results << std::endl;
	
	if( (thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1)
	    )
	  {
	    passTrigger=true;
	    break;
	  }


      }

    if(!passTrigger)continue;
		nPass[1]++;
		
		
		
		}
		
		
		
		
		
			
		}
		
		
		
		
		
		
		cout<<"entries="<<total<<endl;
		//th1->SetBinContent(1,total);
		//string label[5]={"1","2","3","4","5"};
		for(int i=0;i<5;i++){
			//cout<<Form("nPass[%d]=",i)<<double(nPass[i])/total<<endl;
			//th2->SetBinContent(i+1,double(nPass[i])/total);
			//th2->GetXaxis()->SetBinLabel(i+1,Form("%d",i+1));
		}
		
	
	
	
}


//jet energy up
/*
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
*/

//jet energy down
/*
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
*/

//jetEff 
/*
TH1D* HHbbbbAnalyzerJetEff(int wM, string st,string st2,double Xsec,double & db1,double & db2,double & db3,double & db4,double & db5,bool isSignal=0){
	
	
		TFile *f;
		TTree *tree;
	    TH1D* th2 =new TH1D("cut flow","",5,0,5);
		
		db1=0;db2=0;db3=0;db4=0;db5=0;
		
		int nPass[10]={0};
		int total=0;
		
		TH1D* thpt1=new TH1D("fatPt","",50,200,4000);
		TH1D* thpt2=new TH1D("fatPtD","",50,200,4000);
		double bins[10]={0,30,50,70,100,140,200,300,670,2000};
		TH1D* pt11=new TH1D("sub1Pt","",9,bins);
		TH1D* pt12=new TH1D("sub1PtD","",9,bins);
		TH1D* pt21=new TH1D("sub2Pt","",9,bins);
		TH1D* pt22=new TH1D("sub2PtD","",9,bins);
		TH1D* pt31=new TH1D("sub1And2Pt","",9,bins);
		TH1D* pt32=new TH1D("sub1And2PtD","",9,bins);
		
		TH2D* ep11=new TH2D("FatEta_vs_fatPt","",9,bins,50,-2.4,2.4);
		TH2D* ep12=new TH2D("FatEta_vs_fatPtD","",9,bins,50,-2.4,2.4);
		
		TH2D* ep21=new TH2D("FatEta_vs_subPt","",9,bins,50,-2.4,2.4);
		TH2D* ep22=new TH2D("FatEta_vs_subPtD","",9,bins,50,-2.4,2.4);
		
		
		TH1D* shift1=new TH1D("leading","",50,0,0.5);
		TH1D* shift2=new TH1D("subleading","",50,0,0.5);
		
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
		
		//data.Print();
		total+=data.GetEntriesFast();

		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			Int_t nJet         = data.GetInt("FATnJet");
			TClonesArray* FATjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*   FATjetPRmass= data.GetPtrFloat("FATjetPRmass");  
			float*   FATjetPRmassL2L3Corr= data.GetPtrFloat("FATjetPRmassL2L3Corr");  
			vector<bool> &FATjetPassIDLoose=*((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
			vector<float>   *FATsubjetSDCSV = data.GetPtrVectorFloat("FATsubjetSDCSV");
			vector<float>   *FATsubjetSDPx =  data.GetPtrVectorFloat("FATsubjetSDPx");
		    vector<float>   *FATsubjetSDPy =  data.GetPtrVectorFloat("FATsubjetSDPy");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			
			if(nJet<1)continue;
			
			bool b1=0,b2=0;
			
			for(int i=0;i<nJet;i++){
				TLorentzVector* thisHiggs = (TLorentzVector*)FATjetP4->At(i);
				if(thisHiggs->Pt()<200)continue;
				if(fabs(thisHiggs->Eta())>2.4)continue;
				if(FATjetPassIDLoose[i]==0)continue;
				thpt1->Fill(thisHiggs->Pt());
				
				//if(FATsubjetSDHadronFlavor[i][0]!=0||FATsubjetSDHadronFlavor[i][1]!=0)continue;
				
				double temp= sqrt( pow(FATsubjetSDPx[i][0],2)+pow(FATsubjetSDPy[i][0],2));
				double temp2=sqrt( pow(FATsubjetSDPx[i][1],2)+pow(FATsubjetSDPy[i][1],2));
				//FatjetPreSelection.push_back(i);
				pt11->Fill(temp);
				pt21->Fill(temp2);
				pt31->Fill(temp);
				pt31->Fill(temp2);
				ep11->Fill(thisHiggs->Pt(),thisHiggs->Eta());
				ep21->Fill(temp,thisHiggs->Eta());
				ep21->Fill(temp2,thisHiggs->Eta());
				
				
				if(b1==1 && b2==0){
					shift2->Fill((FATjetPRmassL2L3Corr[i]-FATjetPRmass[i])/FATjetPRmass[i]);
					b2=1;
				}
				
				if(b1==0){
					shift1->Fill((FATjetPRmassL2L3Corr[i]-FATjetPRmass[i])/FATjetPRmass[i]);
					b1=1;
				}
				
				
				
				
				if(FATsubjetSDCSV[i][0]>0.605){
					pt12->Fill(temp);
					pt32->Fill(temp);
					ep22->Fill(temp,thisHiggs->Eta());
				}
				
				if(FATsubjetSDCSV[i][1]>0.605){
					pt22->Fill(temp2);
					pt32->Fill(temp2);
					ep22->Fill(temp2,thisHiggs->Eta());
				}
				
				
				if(FATsubjetSDCSV[i][0]<0.605||FATsubjetSDCSV[i][1]<0.605)continue;
				thpt2->Fill(thisHiggs->Pt());
				ep12->Fill(thisHiggs->Pt(),thisHiggs->Eta());
				
				
			}
		
			
		
			
			
		
			
			}
			
		}
		
		//c1 = new TCanvas("c1","",1360,768);
		//tg1->Draw("APL");
		
		//cout<<Form("0114/%s.png",st2.data())<<endl;
		//c1->Print(Form("0114/%s.png",st2.data()));
	
	TFile* outFile = new TFile(Form("0129_light/%s.root",st2.data()),"recreate");
	thpt1->Write();
	thpt2->Write();
	pt11->Write();
	pt21->Write();
	pt31->Write();
	ep11->Write();
	ep21->Write();
	
	pt12->Write();
	pt22->Write();
	pt32->Write();
	ep12->Write();
	ep22->Write();
	//tg2->Write();
	outFile->Close();
	
	TFile* outFile2 = new TFile(Form("0202_G/%s.root",st2.data()),"recreate");
	shift1->Write();
	shift2->Write();
	
	return th2;
	
}
*/

void HHbbbbAnalyzer(int a){
	//setNCUStyle();
	
	string st1[30]={
		/*0-11*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
	    "/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*12-21*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/4c7722e/RadionTohhTohbbhbb/RadionTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
		/*22-29*/
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220208/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220428/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT500to7000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_",
	    "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_",
		"/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_",
	
	};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	string  fileName[30]={"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD100","QCD200","QCD300","QCD500","QCD700","QCD1000","QCD1500","QCD2000"};
	double eff1,eff2,eff3,eff4,eff5;
	TH1D* th1[30];
	int aa[30]={
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		812,516,565,496,
		394,155,103,71
	};
	double Xsec[30]={
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		27990000,1712000,347700,32100,
		6831,1207,119.9,25.24
	};
	bool sigName=0;
	if (a<22)sigName=1;
	
	if(a==100){
	string data1="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-05Oct2015-v1/160223_142842/0000/NCUGlobalTuples_";
	string data2="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-PromptReco-v4/160224_140926/0000/NCUGlobalTuples_";
	string data2_2="/data7/syu/NCUGlobalTuples/Run2015D/eec7461/JetHT/crab_JetHT-Run2015D-PromptReco-v4/160224_140926/0001/NCUGlobalTuples_";
	string dataName1="JetHT-Run2015D-05Oct2015-v1";
	string dataName2="JetHT-Run2015D-PromptReco-v4";
	HHbbbbAnalyzerBaseC(1,412,data1,dataName1,1,0,1);	
	//HHbbbbAnalyzerBaseC(1,501,data2,Form("%s1",dataName2.data()),1,0,1);	
	//HHbbbbAnalyzerBaseC(501,1000,data2,Form("%s2",dataName2.data()),1,0,1);	
	HHbbbbAnalyzerBaseC(1000,1094,data2_2,Form("%s3",dataName2.data()),1,0,1);	
	//dataPrinter(412,data1,dataName1,1,0);	
	//dataPrinter(856,data2,dataName2,1,0);	
	}
	
	bool sig=0;
	if (a<22)sig=1;
	
	//th1[a]=HHbbbbAnalyzerCompare(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerCSVEff(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerJetEff(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerBaseDCut(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	//th1[a]=HHbbbbAnalyzerBaseHPHP(aa[a],st1[a],fileName[a],Xsec[a],eff1,eff2,eff3,eff4,eff5,sig);
	if(!(a==22||a==23||a==24||a==25||a==100))HHbbbbAnalyzerBaseC(1,aa[a],st1[a],fileName[a],Xsec[a],sigName,sig);
	//dataPrinter(aa[a],st1[a],fileName[a],Xsec[a],sig);
	
}