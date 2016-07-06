//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//root feature
#include <TLegend.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <TSystem.h>

//math 
#include <cmath>
#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../untuplizer.h"

void ZPrimeMassAnalyzerBase(int wMs,int wM, string st,string st2,string option=""){	
	//0=signal ,1=QCD ,2=data-----------------------------------------------------------
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	//tuple tree and cutflow variables------------------------------------------------------------------------------------
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0;
	double nPassB[6]={0};
	
	TH1D* th5[4];
	th5[0]=new TH1D("FATjetPRmass","FATjetPRmass",90,0,180);
	th5[1]=new TH1D("FATjetPRmassL2L3Corr","FATjetPRmassL2L3Corr",90,0,180);
	th5[2]=new TH1D("FATjetPuppiSDmass","FATjetPuppiSDmass",90,0,180);
	th5[3]=new TH1D("FATjetPuppiSDmassL2L3Corr","FATjetPuppiSDmassL2L3Corr",90,0,180);
	//NCUtuple loop----------------------------------------------------------------------------------------
	for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		//Get ntuple----------------------------------------------------------------------------------------
		cout<<st.data()<<endl;
		f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		TDirectory * dir;
		dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize(); it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//if(trigResult[it]==1)cout<<it<<","<<thisTrig<<","<<trigResult[it]<<endl;
				if( ((thisTrig.find("PFMET170_NoiseCleaned_v")!= std::string::npos||
						thisTrig.find("PFMET90_PFMHT90_IDTight_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			//if(!passTrigger)continue;
			nPass[0]++;
			Int_t nVtx        = data.GetInt("nVtx");
			if(nVtx<1)continue;
			nPass[1]++;
			TClonesArray* vertexP3 = (TClonesArray*) data.GetPtrTObject("vertexP3");
			TLorentzVector* thisvertex;
			thisvertex= (TLorentzVector*)vertexP3->At(0);
			double x=thisvertex->X(),y=thisvertex->Y(),z=thisvertex->Z();
			if(sqrt(x*x+y*y+z*z)>24)continue;
			if(sqrt(z*z)>2)continue;
			nPass[2]++;
			
			int nFATJet         = data.GetInt("FATnJet");
			if(nFATJet<1)continue;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			TLorentzVector* thisJet;
			thisJet= (TLorentzVector*)fatjetP4->At(0);
			if(thisJet->Pt()<200)continue;
			if(fabs(thisJet->Eta())>2.4)continue;
			nPass[3]++;
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			if(subjetSDCSV[0][0]<0.46)continue;
			if(subjetSDCSV[0][1]<0.46)continue;
			nPass[4]++;
			
			Float_t pfMetCorrPt= data.GetFloat("pfMetCorrPt");
			Float_t pfMetCorrPhi= data.GetFloat("pfMetCorrPhi");
			if(pfMetCorrPt<200)continue;
			nPass[5]++;
			
			TLorentzVector *METLorentz=new TLorentzVector(0,0,0,0);
			METLorentz->SetPtEtaPhiE (pfMetCorrPt,0,pfMetCorrPhi,10);
			
			int THINnJet         = data.GetInt("THINnJet");
			TClonesArray* THINjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
			int* THINjetHadronFlavor=data.GetPtrInt("THINjetHadronFlavor");
			bool AK4JetCheck=0;
			for(int i=0;i<THINnJet;i++){
				TLorentzVector* thatJet;
				thatJet= (TLorentzVector*)THINjetP4->At(0);
				if(thatJet->DeltaR(*thisJet)<0.4){
					if( THINjetHadronFlavor[i]==5&&thatJet->Pt()>30 && fabs(thatJet->DeltaPhi(*METLorentz))<0.4){
						//cout<<thatJet->Phi()<<","<<METLorentz->Phi()<<endl;
						AK4JetCheck=1;
						break;
					}
				}
				else {
					if(THINjetHadronFlavor[i]==5){
						AK4JetCheck=1;
						break;
					}
				}
			}
			if(AK4JetCheck)continue;
			nPass[6]++;
			
			int nFATJetPass=0;
			for(int i=1;i<nFATJet;i++){
				TLorentzVector* thisJeta;
				thisJeta= (TLorentzVector*)fatjetP4->At(i);
				if(thisJeta->Pt()>30 && fabs(thisJeta->Eta())<4.5)nFATJetPass++;
			}
			if(nFATJetPass>1)continue;
			nPass[7]++;
			
			int nMu         = data.GetInt("nMu");
			TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
			int nEle         = data.GetInt("nEle");
			TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
			int HPSTau_n         = data.GetInt("HPSTau_n");
			TClonesArray* HPSTau_4Momentum = (TClonesArray*) data.GetPtrTObject("HPSTau_4Momentum");
			
			bool isoLep=0;
			for(int i=0;i<nMu;i++){
				TLorentzVector* thisJeta;
				thisJeta= (TLorentzVector*)muP4->At(i);
				if(thisJeta->Pt()>10 && fabs(thisJeta->Eta())<2.4 && thisJeta->DeltaR(*thisJet)>0.8)isoLep=1;
				
			}
			if(isoLep)continue;
			nPass[8]++;
			for(int i=0;i<nEle;i++){
				TLorentzVector* thisJeta;
				thisJeta= (TLorentzVector*)eleP4->At(i);
				if(thisJeta->Pt()>10 && fabs(thisJeta->Eta())<2.5 && thisJeta->DeltaR(*thisJet)>0.8)isoLep=1;
				
			}
			if(isoLep)continue;
			nPass[9]++;
			for(int i=0;i<HPSTau_n;i++){
				TLorentzVector* thisJeta;
				thisJeta= (TLorentzVector*)HPSTau_4Momentum->At(i);
				if(thisJeta->Pt()>20 && fabs(thisJeta->Eta())<2.5  && thisJeta->DeltaR(*thisJet)>0.8 )isoLep=1;
				
			}
			if(isoLep)continue;
			nPass[10]++;
			
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			th5[0]->Fill(FATjetPRmass[0]);
			th5[1]->Fill(FATjetPRmassL2L3Corr[0]);
			th5[2]->Fill(FATjetPuppiSDmass[0]);
			th5[3]->Fill(FATjetPuppiSDmassL2L3Corr[0]);
			
			
		}//end event loop----------------------------------------------------------------------------------------
	}	//end ntuple loop----------------------------------------------------------------------------------------
	
	cout<<"total="<<total<<endl;
	for(int i=0;i<16;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TFile* outFile = new TFile(Form("output/%s.root",st2.data()),"recreate");
	th5[0]->Write();
	th5[1]->Write();
	th5[2]->Write();
	th5[3]->Write();
	outFile->Close();
}