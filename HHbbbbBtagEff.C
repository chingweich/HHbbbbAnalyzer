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
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

TCanvas* c1;
/*
void HHbbbbAnalyzerBaseC(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	
	TH2D* th1[6];
	
	th1[0]=new TH2D("effD_b","effD_b",200,0,2000,60,-3,3);
	th1[1]=new TH2D("effN_b","effN_b",200,0,2000,60,-3,3);
	
	th1[2]=new TH2D("effD_c","effD_c",200,0,2000,60,-3,3);
	th1[3]=new TH2D("effN_c","effN_c",200,0,2000,60,-3,3);
	
	th1[4]=new TH2D("effD_l","effD_l",200,0,2000,60,-3,3);
	th1[5]=new TH2D("effN_l","effN_l",200,0,2000,60,-3,3);
	
	
	std::vector<TString> eventlist;                                                                                                                                            
	if(nameRoot==2){
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
		
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
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
	
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		 	const Int_t nsize = data.GetPtrStringSize();

			bool passTrigger=false;
			for(int it=0; it< nsize; it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//cout<<it<<"="<<thisTrig<<endl;
				// std::cout << thisTrig << " : " << results << std::endl;
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			if(!passTrigger && (nameRoot!=0))continue;
			nPass[1]++;
		
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
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
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
			nPass[2]++;
			//3. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[3]++;
			//4. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[4]++;
			//5. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[5]++;
			//6. fatjetPRmassL2L3Corr-----------------------------------------
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[6]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		           tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[7]++;
			
			bool isHPHP=0;
			if(tau21_1<0.6 && tau21_2<0.6 )isHPHP=1;
			
			if(nameRoot==2){
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
			
			
			//8.btag
			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			
			if(FATsubjetSDHadronFlavor[0][0]==5)th1[0]->Fill(thisSub1->Pt(),thisSub1->Eta());
			else if(FATsubjetSDHadronFlavor[0][0]==4)th1[2]->Fill(thisSub1->Pt(),thisSub1->Eta());
			else th1[4]->Fill(thisSub1->Pt(),thisSub1->Eta());
			
			if(FATsubjetSDHadronFlavor[0][1]==5)th1[0]->Fill(thisSub2->Pt(),thisSub2->Eta());
			else if(FATsubjetSDHadronFlavor[0][1]==4)th1[2]->Fill(thisSub2->Pt(),thisSub2->Eta());
			else th1[4]->Fill(thisSub2->Pt(),thisSub2->Eta());
			
			if(FATsubjetSDHadronFlavor[1][0]==5)th1[0]->Fill(thatSub1->Pt(),thatSub1->Eta());
			else if(FATsubjetSDHadronFlavor[1][0]==4)th1[2]->Fill(thatSub1->Pt(),thatSub1->Eta());
			else th1[4]->Fill(thatSub1->Pt(),thatSub1->Eta());
			
			if(FATsubjetSDHadronFlavor[1][1]==5)th1[0]->Fill(thatSub2->Pt(),thatSub2->Eta());
			else if(FATsubjetSDHadronFlavor[1][1]==4)th1[2]->Fill(thatSub2->Pt(),thatSub2->Eta());
			else th1[4]->Fill(thatSub2->Pt(),thatSub2->Eta());
		
			//if(FATsubjetSDHadronFlavor[0][0]==5||FATsubjetSDHadronFlavor[0][0]==4)cout<<"Yes"<<endl;
			cout<<FATsubjetSDHadronFlavor[0][0]<<",";
			
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605){
				if(FATsubjetSDHadronFlavor[0][0]==5)th1[1]->Fill(thisSub1->Pt(),thisSub1->Eta());
				else if(FATsubjetSDHadronFlavor[0][0]==4)th1[3]->Fill(thisSub1->Pt(),thisSub1->Eta());
				else th1[5]->Fill(thisSub1->Pt(),thisSub1->Eta());
			}
			if(subjetSDCSV[0][1]>0.605){
				if(FATsubjetSDHadronFlavor[0][1]==5)th1[1]->Fill(thisSub2->Pt(),thisSub2->Eta());
				else if(FATsubjetSDHadronFlavor[0][1]==4)th1[3]->Fill(thisSub2->Pt(),thisSub2->Eta());
				else th1[5]->Fill(thisSub2->Pt(),thisSub2->Eta());
			}
			if(subjetSDCSV[1][0]>0.605){
				if(FATsubjetSDHadronFlavor[1][0]==5)th1[1]->Fill(thatSub1->Pt(),thatSub1->Eta());
				else if(FATsubjetSDHadronFlavor[1][0]==4)th1[3]->Fill(thatSub1->Pt(),thatSub1->Eta());
				else th1[5]->Fill(thatSub1->Pt(),thatSub1->Eta());
			}
			if(subjetSDCSV[1][1]>0.605){
				if(FATsubjetSDHadronFlavor[1][1]==5)th1[1]->Fill(thatSub2->Pt(),thatSub2->Eta());
				else if(FATsubjetSDHadronFlavor[1][1]==4)th1[3]->Fill(thatSub2->Pt(),thatSub2->Eta());
				else th1[5]->Fill(thatSub2->Pt(),thatSub2->Eta());
			}
			
			
	
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(nameRoot==2)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<14;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TFile* outFile = new TFile(Form("root_files_btaggedEff/%s.root",st2.data()),"recreate");
	
	for (int i=0;i<6;i++)th1[i]->Write();
	outFile->Close();
	
}
*/


void HHbbbbAnalyzerBaseC(int wMs,int wM, string st,string st2,double Xsec,int nameRoot=0){	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("root_files_btaggedEff/data.root");
	else f1=TFile::Open(Form("root_files_btaggedEff/%s.root",st2.data()));
	TH2D* th1[6];
	th1[0]=(TH2D*)f1->FindObjectAny("effD_b");
	th1[1]=(TH2D*)f1->FindObjectAny("effN_b");
	th1[2]=(TH2D*)f1->FindObjectAny("effD_c");
	th1[3]=(TH2D*)f1->FindObjectAny("effN_c");
	th1[4]=(TH2D*)f1->FindObjectAny("effD_l");
	th1[5]=(TH2D*)f1->FindObjectAny("effN_l");
	
	th1[1]->Divide(th1[0]);
	th1[3]->Divide(th1[2]);
	th1[5]->Divide(th1[4]);
	
	TFile *f2;
	f2=TFile::Open("root_files_btaggedEff/data.root");
	TH2D* th2[6];
	th2[0]=(TH2D*)f1->FindObjectAny("effD_b");
	th2[1]=(TH2D*)f1->FindObjectAny("effN_b");
	th2[2]=(TH2D*)f1->FindObjectAny("effD_c");
	th2[3]=(TH2D*)f1->FindObjectAny("effN_c");
	th2[4]=(TH2D*)f1->FindObjectAny("effD_l");
	th2[5]=(TH2D*)f1->FindObjectAny("effN_l");
	
	th2[1]->Divide(th2[0]);
	th2[3]->Divide(th2[2]);
	th2[5]->Divide(th2[4]);
	
	
	
	BTagCalibration calib("CSVv2L", "CSVv2.csv");
	BTagCalibrationReader LF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "comb",               // measurement type
                             "central");           // systematics type
	
	BTagCalibrationReader HF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "mujets",               // measurement type
                             "central");           // systematics type
	
	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};int total=0;
	int dataPassingcsc=0;
	double nPassB[5]={0};
	
	for (int w=wMs;w<wM;w++){
		
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
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
	
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		 	const Int_t nsize = data.GetPtrStringSize();

			bool passTrigger=false;
			for(int it=0; it< nsize; it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//cout<<it<<"="<<thisTrig<<endl;
				// std::cout << thisTrig << " : " << results << std::endl;
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			if(!passTrigger )continue;
			nPass[1]++;
		
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
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
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
			nPass[2]++;
			//3. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[3]++;
			//4. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[4]++;
			//5. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			if(mjj<1000)continue;
			nPass[5]++;
			//6. fatjetPRmassL2L3Corr-----------------------------------------
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[6]++;
			//7.-----------------------------------------
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),
		           tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			if(tau21_1>0.75||tau21_2>0.75)continue;
			if(tau21_1>0.6 &&tau21_2>0.6) continue;
			nPass[7]++;
			
			bool isHPHP=0;
			if(tau21_1<0.6 && tau21_2<0.6 )isHPHP=1;
			
			
			//8.btag

			TLorentzVector* thisSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thisSub2=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub1=new TLorentzVector(0,0,0,0);
			TLorentzVector* thatSub2=new TLorentzVector(0,0,0,0);
			thisSub1->SetPxPyPzE(subjetSDPx[0][0],subjetSDPy[0][0],subjetSDPz[0][0],subjetSDE[0][0]);
			thisSub2->SetPxPyPzE(subjetSDPx[0][1],subjetSDPy[0][1],subjetSDPz[0][1],subjetSDE[0][1]);
			thatSub1->SetPxPyPzE(subjetSDPx[1][0],subjetSDPy[1][0],subjetSDPz[1][0],subjetSDE[1][0]);
			thatSub2->SetPxPyPzE(subjetSDPx[1][1],subjetSDPy[1][1],subjetSDPz[1][1],subjetSDE[1][1]);
			
			
			double sf1=1,sf2=1,sf3=1,sf4=1;
			
			
			float MaxBJetPt = 670., MaxLJetPt = 1000.;
			
			double pt1=thisSub1->Pt(),pt2=thisSub2->Pt(),pt3=thatSub1->Pt(),pt4=thatSub2->Pt();
			if(FATsubjetSDHadronFlavor[0][0]!=0 && pt1>MaxBJetPt )pt1=MaxBJetPt;
			if(FATsubjetSDHadronFlavor[0][0]==0 && pt1>MaxLJetPt )pt1=MaxLJetPt;
			
			if(FATsubjetSDHadronFlavor[0][0]!=0 && pt2>MaxBJetPt )pt2=MaxBJetPt;
			if(FATsubjetSDHadronFlavor[0][0]==0 && pt2>MaxLJetPt )pt2=MaxLJetPt;
			
			if(FATsubjetSDHadronFlavor[0][0]!=0 && pt3>MaxBJetPt )pt3=MaxBJetPt;
			if(FATsubjetSDHadronFlavor[0][0]==0 && pt3>MaxLJetPt )pt3=MaxLJetPt;
			
			if(FATsubjetSDHadronFlavor[0][0]!=0 && pt4>MaxBJetPt )pt4=MaxBJetPt;
			if(FATsubjetSDHadronFlavor[0][0]==0 && pt4>MaxLJetPt )pt4=MaxLJetPt;
			
			//cout<<pt1<<","<<pt2<<","<<pt3<<","<<pt4<<endl;
			//cout<<thisSub1->Eta()<<","<<thisSub2->Eta()<<","<<thatSub1->Eta()<<","<<thatSub2->Eta()<<endl;
			
			if(FATsubjetSDHadronFlavor[0][0]==5)sf1=HF.eval(BTagEntry::FLAV_B,thisSub1->Eta(),pt1); 
			else if(FATsubjetSDHadronFlavor[0][0]==4)sf1=HF.eval(BTagEntry::FLAV_C,thisSub1->Eta(),pt1); 
			else sf1=LF.eval(BTagEntry::FLAV_UDSG,thisSub1->Eta(),pt1); 
			
			//cout<<"1"<<endl;
			
			if(FATsubjetSDHadronFlavor[0][1]==5)sf2=HF.eval(BTagEntry::FLAV_B,thisSub2->Eta(),pt2); 
			else if(FATsubjetSDHadronFlavor[0][1]==4)sf2=HF.eval(BTagEntry::FLAV_C,thisSub2->Eta(),pt2); 
			else sf2=LF.eval(BTagEntry::FLAV_UDSG,thisSub2->Eta(),pt2); 
			
			//cout<<"2"<<endl;
			
			if(FATsubjetSDHadronFlavor[1][0]==5)sf3=HF.eval(BTagEntry::FLAV_B,thatSub1->Eta(),pt3); 
			else if(FATsubjetSDHadronFlavor[1][0]==4)sf3=HF.eval(BTagEntry::FLAV_C,thatSub1->Eta(),pt3); 
			else sf3=LF.eval(BTagEntry::FLAV_UDSG,thatSub1->Eta(),pt3); 
			
			//cout<<"3"<<endl;
			
			if(FATsubjetSDHadronFlavor[1][1]==5)sf4=HF.eval(BTagEntry::FLAV_B,thatSub2->Eta(),pt4); 
			else if(FATsubjetSDHadronFlavor[1][1]==4)sf4=HF.eval(BTagEntry::FLAV_C,thatSub2->Eta(),pt4); 
			else sf4=LF.eval(BTagEntry::FLAV_UDSG,thatSub2->Eta(),pt4); 
			
			//cout<<"4"<<endl;
			
			double eff1,eff2,eff3,eff4;

			if(FATsubjetSDHadronFlavor[0][0]==5)eff1=th1[1]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[0][0]==4)eff1=th1[3]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			else eff1=th1[5]->GetBinContent(ceil(thisSub1->Pt()/10),ceil(thisSub1->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[0][1]==5)eff2=th1[1]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[0][1]==4)eff2=th1[3]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			else eff2=th1[5]->GetBinContent(ceil(thisSub2->Pt()/10),ceil(thisSub2->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[1][0]==5)eff3=th1[1]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[1][0]==4)eff3=th1[3]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			else eff3=th1[5]->GetBinContent(ceil(thatSub1->Pt()/10),ceil(thatSub1->Eta()/0.1)+30);
			
			if(FATsubjetSDHadronFlavor[1][1]==5)eff4=th1[1]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			else if(FATsubjetSDHadronFlavor[1][1]==4)eff4=th1[3]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			else eff4=th1[5]->GetBinContent(ceil(thatSub2->Pt()/10),ceil(thatSub2->Eta()/0.1)+30);
			
			double scaleFactor=1;
			if(subjetSDCSV[0][0]>0.605)scaleFactor*=sf1;
			else scaleFactor*=((1-eff1*sf1)/(1-eff1));
			if(subjetSDCSV[0][1]>0.605)scaleFactor*=sf2;
			else scaleFactor*=((1-eff2*sf2)/(1-eff2));
			if(subjetSDCSV[1][0]>0.605)scaleFactor*=sf3;
			else scaleFactor*=((1-eff3*sf3)/(1-eff3));
			if(subjetSDCSV[1][1]>0.605)scaleFactor*=sf4;
			else scaleFactor*=((1-eff4*sf4)/(1-eff4));
			
			/*
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[0][0]>0.605)cout<<"1Y="<<sf1<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[0][0]<0.605)cout<<"1N="<<eff1<<"=="<<((1-eff1*sf1)/(1-eff1))<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[0][1]>0.605)cout<<"2Y="<<sf2<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[0][1]<0.605)cout<<"2N="<<eff2<<"=="<<((1-eff2*sf2)/(1-eff2))<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[1][0]>0.605)cout<<"3Y="<<sf3<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[1][0]<0.605)cout<<"3N="<<eff3<<"=="<<((1-eff3*sf3)/(1-eff3))<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[1][1]>0.605)cout<<"4Y="<<sf4<<",";
			if((scaleFactor<0.1||scaleFactor>2) && subjetSDCSV[1][1]<0.605)cout<<"4N="<<eff4<<"=="<<((1-eff4*sf4)/(1-eff4))<<",";
			if(scaleFactor<0.1||scaleFactor>2)cout<<scaleFactor<<endl;
			*/
			
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.605)nbtag++;
			if(subjetSDCSV[0][1]>0.605)nbtag++;
			if(subjetSDCSV[1][0]>0.605)nbtag++;
			if(subjetSDCSV[1][1]>0.605)nbtag++;
			
			if(nbtag==0)nPassB[0]+=scaleFactor;
			if(nbtag==1)nPassB[1]+=scaleFactor;
			if(nbtag==2)nPassB[2]+=scaleFactor;
			if(nbtag==3)nPassB[3]+=scaleFactor;
			if(nbtag==4)nPassB[4]+=scaleFactor;
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	if(nameRoot==2)cout<<"dataPassingcsc="<<dataPassingcsc<<endl;
	for(int i=0;i<14;i++)cout<<"nPass["<<i<<"]="<<nPass[i]<<endl;
	
	TH1D * th2o=new TH1D("Nbtagjet","Nbtagjet",5,-0.5,4.5);
	th2o->SetBinContent(1,nPassB[0]);
	th2o->SetBinContent(2,nPassB[1]);
	th2o->SetBinContent(3,nPassB[2]);
	if(nameRoot!=2)th2o->SetBinContent(4,nPassB[3]);
	if(nameRoot!=2)th2o->SetBinContent(5,nPassB[4]);
	
	TH1D * th2s=new TH1D("NbtagjetS","NbtagjetS",5,-0.5,4.5);
	th2s->SetBinContent(1,nPassB[0]);
	th2s->SetBinContent(2,nPassB[1]);
	th2s->SetBinContent(3,nPassB[2]);
	if(nameRoot!=2)th2s->SetBinContent(4,nPassB[3]);
	if(nameRoot!=2)th2s->SetBinContent(5,nPassB[4]);
	
	TFile* outFile = new TFile(Form("root_files_btaggedEff/%s_2.root",st2.data()),"recreate");
	th2o->Write();
	th2s->Scale(2245.87*Xsec/total);
	th2s->Write();
	outFile->Close();
	
}

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
			for(int it=0; it< nsize; it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				cout<<it<<"="<<thisTrig<<endl;
				 //std::cout << thisTrig << " : " << results << std::endl;
				
				
			}
		}
		cout<<"entries="<<total<<endl;
		}
}



void HHbbbbBtagEff(int a){
	//setNCUStyle();
	
	//setNCUStyle();
	
	string st1[40]={
		/*0-11*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph.root",
	    "/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph.root",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/BulkGravTohhTohbbhbb/BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph.root",
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
		/*22-30*/
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060112/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060112/0001/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060145/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060215/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160306_060251/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130001/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130035/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130118/0000/NCUGlobalTuples_",
		"/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/f81355f/QCD_HTBinned/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160308_130153/0000/NCUGlobalTuples_",
		
	};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	string  fileName[40]={"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD100_1","QCD100_2","QCD200","QCD300","QCD500","QCD700","QCD1000","QCD1500","QCD2000"};
	double eff1,eff2,eff3,eff4,eff5;
	TH1D* th1[30];
	int aa0[40]={
		1,1,1,1,1,
		1,1,1,1,1,
		1,
		1,1,1,1,1,
		1,1,1,1,1,
		1,
		1,1000,1,1,1,
		1,1,1,1,
	};
	
	int aa[40]={
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		2,2,2,2,2,
		2,2,2,2,2,
		2,
		1000,1904,432,467,449,
		346,120,90,45,
		
	};
	double Xsec[40]={
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,0.1,0.1,0.1,0.1,
		0.1,
		27990000,27990000,1712000,347700,32100,
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
	//HHbbbbAnalyzerBaseC(1,412,data1,dataName1,1,2);	
	HHbbbbAnalyzerBaseC(1,501,data2,Form("%s1",dataName2.data()),1,2);	
	HHbbbbAnalyzerBaseC(501,1000,data2,Form("%s2",dataName2.data()),1,2);	
	//HHbbbbAnalyzerBaseC(1000,1094,data2_2,Form("%s3",dataName2.data()),1,2);	
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
	if(!(a==22||a==23||a==24||a==25||a==26||a==100))HHbbbbAnalyzerBaseC(aa0[a],aa[a],st1[a],fileName[a],Xsec[a],sigName);
	//dataPrinter(aa[a],st1[a],fileName[a],Xsec[a],sig);
	
}