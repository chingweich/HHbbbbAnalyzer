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
#include <TF1.h>
//math 
#include <cmath>
#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "untuplizer.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"


void miniTreeBase(int wMs,int wM, string st,string st2,string option="") {
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0;
	
	
	//option-----------------------------------------------------------
	int JESOption=0;
	if(option.find("JESUp")!= std::string::npos)JESOption=1;
	if(option.find("JESDown")!= std::string::npos)JESOption=2;
	if(option.find("BtagUp")!= std::string::npos)JESOption=3;
	if(option.find("BtagDown")!= std::string::npos)JESOption=4;
	if(option.find("tau21Up")!= std::string::npos)JESOption=5;
	if(option.find("tau21Down")!= std::string::npos)JESOption=6;
	cout<<"JESOption = "<<JESOption<<endl;
	
	int nameRoot=1;
	if(st2.find("QCD")!= std::string::npos)nameRoot=0;
	if(st2.find("bGen")!= std::string::npos)nameRoot=0;
	if(st2.find("bEnriched")!= std::string::npos)nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	
	//-------------------------btag
	string btagsystematicsType="central";
	if(JESOption==3)btagsystematicsType="up";
	else if(JESOption==4)btagsystematicsType="down";
	BTagCalibration calib("CSVv2L", "subjet_CSVv2_ichep.csv");
	BTagCalibrationReader LF(BTagEntry::OP_LOOSE,"central", {"up", "down"});  
	BTagCalibrationReader HFC(BTagEntry::OP_LOOSE, "central", {"up", "down"});      // other sys types
	BTagCalibrationReader HF(BTagEntry::OP_LOOSE,"central",{"up", "down"});      // other sys types
	LF.load(calib, BTagEntry::FLAV_UDSG,    // btag flavour
            "incl");               // measurement type
	HFC.load(calib, BTagEntry::FLAV_C,    // btag flavour
            "lt");               // measurement type
	HF.load(calib, BTagEntry::FLAV_B,    // btag flavour
            "lt");               // measurement type
	
	TFile *f1;
	if(nameRoot==2)f1=TFile::Open("btagEffSource/data.root");
	else if (nameRoot!=2 && (JESOption==0||JESOption==3||JESOption==4||JESOption==5||JESOption==6))f1=TFile::Open(Form("btagEffSource/%s.root",st2.data()));
	else if (nameRoot!=2 && JESOption==1)f1=TFile::Open(Form("btagEffSource/%s_JESUp.root",st2.data()));
	else if (nameRoot!=2 && JESOption==2)f1=TFile::Open(Form("btagEffSource/%s_JESDown.root",st2.data()));
	TH1D* th1[6];
	string btaggingEff[6]={"effD_b","effN_b","effD_c","effN_c","effD_l","effN_l"};
	for(int i=0;i<6;i++){
		th1[i]=(TH1D*)f1->FindObjectAny(Form("%s_1d",btaggingEff[i].data()));
		if(i==1||i==3||i==5)th1[i]->Divide(th1[i-1]);
	}
	//-------------------------btag
	
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
	
	//NCUtuple loop----------------------------------------------------------------------------------------
	for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		//Get ntuple----------------------------------------------------------------------------------------
		cout<<st.data()<<endl;
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		TDirectory * dir;
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		
		TFile* outFile = new TFile(Form("%s.root",st2.data()),"recreate");
		TTree* hh4b=new TTree("hh4b","hh4b");
	
		Int_t SelectedEvent_runno;                        hh4b->Branch("SelectedEvent_runno",&SelectedEvent_runno,"SelectedEvent_runno/I");
		Int_t SelectedEvent_lumisec;                      hh4b->Branch("SelectedEvent_lumisec",&SelectedEvent_lumisec,"SelectedEvent_lumisec/I");
		Int_t SelectedEvent_evtno;                        hh4b->Branch("SelectedEvent_evtno",&SelectedEvent_evtno,"SelectedEvent_evtno/I");
		Int_t SelectedEvent_bcno;                         hh4b->Branch("SelectedEvent_bcno",&SelectedEvent_bcno,"SelectedEvent_bcno/I");
		Int_t SelectedEvent_time;                         hh4b->Branch("SelectedEvent_time",&SelectedEvent_time,"SelectedEvent_time/I");
		
		Double_t SelectedEvent_btagsf;                    hh4b->Branch("SelectedEvent_btagsf",&SelectedEvent_btagsf,"SelectedEvent_btagsf/D");
		Double_t SelectedEvent_btagsf_bcUp;               hh4b->Branch("SelectedEvent_btagsf_bcUp",&SelectedEvent_btagsf_bcUp,"SelectedEvent_btagsf_bcUp/D");
		Double_t SelectedEvent_btagsf_bcDown;             hh4b->Branch("SelectedEvent_btagsf_bcDown",&SelectedEvent_btagsf_bcDown,"SelectedEvent_btagsf_bcDown/D");
		Double_t SelectedEvent_btagsf_lUp;                hh4b->Branch("SelectedEvent_btagsf_lUp",&SelectedEvent_btagsf_lUp,"SelectedEvent_btagsf_lUp/D");
		Double_t SelectedEvent_btagsf_lDown;              hh4b->Branch("SelectedEvent_btagsf_lDown",&SelectedEvent_btagsf_lDown,"SelectedEvent_btagsf_lDown/D");
	
		
		Int_t SelectedEvent_nsubjetsBTaggedCSVL;          hh4b->Branch("SelectedEvent_nsubjetsBTaggedCSVL",&SelectedEvent_nsubjetsBTaggedCSVL,"SelectedEvent_nsubjetsBTaggedCSVL/I");
		Int_t HJets_njets;                                hh4b->Branch("HJets_njets",&HJets_njets,"HJets_njets/I");
		Int_t SelectedEvent_npv;                          hh4b->Branch("SelectedEvent_npv",&SelectedEvent_npv,"SelectedEvent_npv/I");
		Int_t SelectedEvent_npuTrue;                      hh4b->Branch("SelectedEvent_npuTrue",&SelectedEvent_npuTrue,"SelectedEvent_npuTrue/I");
		Int_t SelectedEvent_ht;                           hh4b->Branch("SelectedEvent_ht",&SelectedEvent_ht,"SelectedEvent_ht/I");
		
		Double_t SelectedEvent_htHat;                     hh4b->Branch("SelectedEvent_htHat",&SelectedEvent_htHat,"SelectedEvent_htHat/D");
		Double_t SelectedEvent_deta_leading2hjets;        hh4b->Branch("SelectedEvent_deta_leading2hjets",&SelectedEvent_deta_leading2hjets,"SelectedEvent_deta_leading2hjets/D");
		Double_t SelectedEvent_pt_leading2hjets;          hh4b->Branch("SelectedEvent_pt_leading2hjets",&SelectedEvent_pt_leading2hjets,"SelectedEvent_pt_leading2hjets/D");
		Double_t SelectedEvent_eta_leading2hjets;         hh4b->Branch("SelectedEvent_eta_leading2hjets",&SelectedEvent_eta_leading2hjets,"SelectedEvent_eta_leading2hjets/D");
		Double_t SelectedEvent_phi_leading2hjets;         hh4b->Branch("SelectedEvent_phi_leading2hjets",&SelectedEvent_phi_leading2hjets,"SelectedEvent_phi_leading2hjets/D");
		Double_t SelectedEvent_y_leading2hjets;           hh4b->Branch("SelectedEvent_y_leading2hjets",&SelectedEvent_y_leading2hjets,"SelectedEvent_y_leading2hjets/D");
		
		Double_t SelectedEvent_evtwt;                     hh4b->Branch("SelectedEvent_evtwt",&SelectedEvent_evtwt,"SelectedEvent_evtwt/D");
		Double_t SelectedEvent_evtwtPV;                   hh4b->Branch("SelectedEvent_evtwtPV",&SelectedEvent_evtwtPV,"SelectedEvent_evtwtPV/D");
		Double_t SelectedEvent_evtwtPVLow;                hh4b->Branch("SelectedEvent_evtwtPVLow",&SelectedEvent_evtwtPVLow,"SelectedEvent_evtwtPVLow/D");
		Double_t SelectedEvent_evtwtPVHigh;               hh4b->Branch("SelectedEvent_evtwtPVHigh",&SelectedEvent_evtwtPVHigh,"SelectedEvent_evtwtPVHigh/D");
		Double_t SelectedEvent_minv_leading2hjets;        hh4b->Branch("SelectedEvent_minv_leading2hjets",&SelectedEvent_minv_leading2hjets,"SelectedEvent_minv_leading2hjets/D");
		Double_t SelectedEvent_minv_leading2hjets_subtr;  hh4b->Branch("SelectedEvent_minv_leading2hjets_subtr",&SelectedEvent_minv_leading2hjets_subtr,"SelectedEvent_minv_leading2hjets_subtr/D");
	
		
		Int_t HJets_nsubjets[HJets_njets];                hh4b->Branch("HJets_nsubjets",&HJets_nsubjets,"HJets_nsubjets[HJets_njets]/I");
	
		Double_t HJets_Pt[HJets_njets];                   hh4b->Branch("HJets_Pt",&HJets_Pt,"HJets_Pt[HJets_njets]/D");
		Double_t HJets_Eta[HJets_njets];                  hh4b->Branch("HJets_Eta",&HJets_Eta,"HJets_Eta[HJets_njets]/D");
		Double_t HJets_Phi[HJets_njets];                  hh4b->Branch("HJets_Phi",&HJets_Phi,"HJets_Phi[HJets_njets]/D");
		Double_t HJets_Energy[HJets_njets];               hh4b->Branch("HJets_Energy",&HJets_Energy,"HJets_Energy[HJets_njets]/D");
		Double_t HJets_Mass[HJets_njets];                 hh4b->Branch("HJets_Mass",&HJets_Mass,"HJets_Mass[HJets_njets]/D");
		Double_t HJets_tau1[HJets_njets];                 hh4b->Branch("HJets_tau1",&HJets_tau1,"HJets_tau1[HJets_njets]/D");
		Double_t HJets_tau2[HJets_njets];                 hh4b->Branch("HJets_tau2",&HJets_tau2,"HJets_tau2[HJets_njets]/D");
		Double_t HJets_tau3[HJets_njets];                 hh4b->Branch("HJets_tau3",&HJets_tau3,"HJets_tau3[HJets_njets]/D");
		Double_t HJets_hadFlavourSubjet0[HJets_njets];    hh4b->Branch("HJets_hadFlavourSubjet0",&HJets_hadFlavourSubjet0,"HJets_hadFlavourSubjet0[HJets_njets]/D");
		Double_t HJets_hadFlavourSubjet1[HJets_njets];    hh4b->Branch("HJets_hadFlavourSubjet1",&HJets_hadFlavourSubjet1,"HJets_hadFlavourSubjet1[HJets_njets]/D");
		Double_t HJets_ptSubjet0[HJets_njets];            hh4b->Branch("HJets_ptSubjet0",&HJets_ptSubjet0,"HJets_ptSubjet0[HJets_njets]/D");
		Double_t HJets_ptSubjet1[HJets_njets];            hh4b->Branch("HJets_ptSubjet1",&HJets_ptSubjet1,"HJets_ptSubjet1[HJets_njets]/D");
		Double_t HJets_etaSubjet0[HJets_njets];           hh4b->Branch("HJets_etaSubjet0",&HJets_etaSubjet0,"HJets_etaSubjet0[HJets_njets]/D");
		Double_t HJets_etaSubjet1[HJets_njets];           hh4b->Branch("HJets_etaSubjet1",&HJets_etaSubjet1,"HJets_etaSubjet1[HJets_njets]/D");
		Double_t HJets_phiSubjet0[HJets_njets];           hh4b->Branch("HJets_phiSubjet0",&HJets_phiSubjet0,"HJets_phiSubjet0[HJets_njets]/D");
		Double_t HJets_phiSubjet1[HJets_njets];           hh4b->Branch("HJets_phiSubjet1",&HJets_etaSubjet1,"HJets_phiSubjet1[HJets_njets]/D");
		Double_t HJets_energySubjet0[HJets_njets];        hh4b->Branch("HJets_energySubjet0",&HJets_energySubjet0,"HJets_energySubjet0[HJets_njets]/D");
		Double_t HJets_energySubjet1[HJets_njets];        hh4b->Branch("HJets_energySubjet1",&HJets_energySubjet1,"HJets_energySubjet1[HJets_njets]/D");
		Double_t HJets_CSVIVFv2[HJets_njets];             hh4b->Branch("HJets_CSVIVFv2",&HJets_CSVIVFv2,"HJets_CSVIVFv2[HJets_njets]/D");
		Double_t HJets_MassPruned[HJets_njets];           hh4b->Branch("HJets_MassPruned",&HJets_MassPruned,"HJets_MassPruned[HJets_njets]/D");
		Double_t HJets_MassSoftDrop[HJets_njets];         hh4b->Branch("HJets_MassSoftDrop",&HJets_MassSoftDrop,"HJets_MassSoftDrop[HJets_njets]/D");
		Double_t HJets_csvSubjet0[HJets_njets];           hh4b->Branch("HJets_csvSubjet0",&HJets_csvSubjet0,"HJets_csvSubjet0[HJets_njets]/D");
		Double_t HJets_csvSubjet1[HJets_njets];           hh4b->Branch("HJets_csvSubjet1",&HJets_csvSubjet1,"HJets_csvSubjet1[HJets_njets]/D");
		Double_t HJets_groomedMassCorr[HJets_njets];      hh4b->Branch("HJets_groomedMassCorr",&HJets_groomedMassCorr,"HJets_groomedMassCorr[HJets_njets]/D");
		Double_t HJets_nhf[HJets_njets];                  hh4b->Branch("HJets_nhf",&HJets_nhf,"HJets_nhf[HJets_njets]/D");
		Double_t HJets_chf[HJets_njets];                  hh4b->Branch("HJets_chf",&HJets_chf,"HJets_chf[HJets_njets]/D");
		Double_t HJets_emf[HJets_njets];                  hh4b->Branch("HJets_emf",&HJets_emf,"HJets_emf[HJets_njets]/D");
		Double_t HJets_phf[HJets_njets];                  hh4b->Branch("HJets_phf",&HJets_phf,"HJets_phf[HJets_njets]/D");
		Double_t HJets_muf[HJets_njets];                  hh4b->Branch("HJets_muf",&HJets_muf,"HJets_muf[HJets_njets]/D");
		Double_t HJets_hadFlavour[HJets_njets];           hh4b->Branch("HJets_hadFlavour",&HJets_hadFlavour,"HJets_hadFlavour[HJets_njets]/D");
		Double_t HJets_doublesv[HJets_njets];             hh4b->Branch("HJets_doublesv",&HJets_doublesv,"HJets_doublesv[HJets_njets]/D");
		Int_t HJets_nsubjetsBTaggedCSVL[HJets_njets];     hh4b->Branch("HJets_nsubjetsBTaggedCSVL",&HJets_nsubjetsBTaggedCSVL,"HJets_nsubjetsBTaggedCSVL[HJets_njets]/I");
		Int_t HJets_Index[HJets_njets];                   hh4b->Branch("HJets_Index",&HJets_Index,"HJets_Index[HJets_njets]/I");
		Int_t HJets_ParentIndex[HJets_njets];             hh4b->Branch("HJets_ParentIndex",&HJets_ParentIndex,"HJets_ParentIndex[HJets_njets]/I");
	
	
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			Int_t nVtx        = data.GetInt("nVtx");
			
			//0. has a good vertex
			if(nVtx<1)continue;nPass[0]++;
			//1.trigger
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize(); it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
				//if(trigResult[it]==1)cout<<it<<","<<thisTrig<<","<<trigResult[it]<<endl;
				if( ((thisTrig.find("HLT_PFHT800")!= std::string::npos||
						//thisTrig.find("HLT_PFHT650")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos||
						thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos
						) && results==1)){
					passTrigger=true;
					break;
				}
			}
			//if(!passTrigger)continue;
			nPass[1]++;
		
			int nFATJet         = data.GetInt("FATnJet");
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			//2.nJets
			if(nFATJet<2)continue;nPass[2]++;
			TLorentzVector* thisJet ,* thatJet;
			thisJet= (TLorentzVector*)fatjetP4->At(0);
			thatJet = (TLorentzVector*)fatjetP4->At(1);
			//3. Pt 
			if(thisJet->Pt()<200||thatJet->Pt()<200)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			if(FATjetPassIDTight[0]==0||FATjetPassIDTight[1]==0)continue;
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuEF = data.GetPtrFloat("FATjetMuEF");
			if(FATjetMuEF[0]>0.8||FATjetMuEF[1]>0.8)continue;
			if(FATjetCEmEF[0]>0.9||FATjetCEmEF[1]>0.9)continue;
			nPass[4]++;
			//5. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4||fabs(thatJet->Eta())>2.4)continue;
			nPass[5]++;
			//6. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			//if(dEta>1.3)continue;
			nPass[6]++;
			//7. Mjj-----------------------------------------
			float mjj = (*thisJet+*thatJet).M();
			float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			if(mjjRed<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			nPass[8]++;
			//9.-----------------------------------------
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			Float_t*  fatjetTau3 = data.GetPtrFloat("FATjetTau3");
			Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
			Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
			
			
			int nbtag=0;
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			if(subjetSDCSV[0][0]>0.46)nbtag++;
			if(subjetSDCSV[0][1]>0.46)nbtag++;
			if(subjetSDCSV[1][0]>0.46)nbtag++;
			if(subjetSDCSV[1][1]>0.46)nbtag++;
			nPass[10]++;
			
			
			Long64_t runId        = data.GetLong64("runId");
			Long64_t lumiSection        = data.GetLong64("lumiSection");
			Long64_t eventId        = data.GetLong64("eventId");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJet);
			//vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			Float_t*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Int_t*  FATjetHadronFlavor = data.GetPtrInt("FATjetHadronFlavor");
			Int_t*  FATnSubSDJet = data.GetPtrInt("FATnSubSDJet");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Float_t*  FATjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
			
			double sf[2][2],subjetPt[2][2],subjetEta[2][2],eff[2][2],btaggingscaleFactor=1;
			TLorentzVector* subjetP4[2][2];
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					sf[i][j]=1;
					subjetP4[i][j]=new TLorentzVector(0,0,0,0);
					subjetP4[i][j]->SetPxPyPzE(subjetSDPx[i][j],subjetSDPy[i][j],subjetSDPz[i][j],subjetSDE[i][j]);
					subjetPt[i][j]=subjetP4[i][j]->Pt();
					subjetEta[i][j]=subjetP4[i][j]->Eta();
				}
			}
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					
					//get btagging eff------------------------------------------------------------
					int getPtBin=1;
					if(subjetPt[i][j]<140)getPtBin=1;
					else if (140<=subjetPt[i][j] && subjetPt[i][j]<180)getPtBin=2;
					else if (180<=subjetPt[i][j] && subjetPt[i][j]<240)getPtBin=3;
					else getPtBin=4;
					if(FATsubjetSDHadronFlavor[i][j]==5)eff[i][j]=th1[1]->GetBinContent(getPtBin);
					else if(FATsubjetSDHadronFlavor[i][j]==4)eff[i][j]=th1[3]->GetBinContent(getPtBin);
					else {
						int temp=0;
						if(subjetPt[i][j]>=900)temp=10;
						else temp=ceil(subjetPt[i][j]/100);
						
						bool checkBinContentIfZero=0;
						while(checkBinContentIfZero==0){
							if(th1[4]->GetBinContent(temp)==0){
								temp--;
							}
							else checkBinContentIfZero=1;
						}
						eff[i][j]=th1[5]->GetBinContent(temp);
					}
					//Get SF from csv------------------------------------------------------------
					if(FATsubjetSDHadronFlavor[i][j]==5){
						sf[i][j]=HF.eval_auto_bounds("central",BTagEntry::FLAV_B, subjetEta[i][j],subjetPt[i][j]); 
					}
					else if(FATsubjetSDHadronFlavor[i][j]==4){
						sf[i][j]=HFC.eval_auto_bounds("central",BTagEntry::FLAV_C, subjetEta[i][j],subjetPt[i][j]); 
					}
					else {
						sf[i][j]=LF.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, subjetEta[i][j],subjetPt[i][j]); 
					}
					//get tot. btagging SF
					if(subjetSDCSV[i][j]>=0.46)btaggingscaleFactor*=sf[i][j];
					else btaggingscaleFactor*=((1-eff[i][j]*sf[i][j])/(1-eff[i][j]));
					
				}
			}
			
			
			Float_t*  ADDjet_DoubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
			Float_t  FATjet_DoubleSV[2]={0};
			
			Int_t ADDnJet        = data.GetInt("ADDnJet");
			bool matchThis=0,matchThat=0;
			TClonesArray* ADDjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
			for(int i=0;i<ADDnJet;i++){
				TLorentzVector* thisAddJet ;
				thisAddJet= (TLorentzVector*)ADDjetP4->At(i);
				if(!matchThis && thisAddJet->DeltaR(*thisJet)<0.8){
					matchThis=1;
					FATjet_DoubleSV[0]=ADDjet_DoubleSV[i];
					continue;
				}
				if(!matchThat && thisAddJet->DeltaR(*thatJet)<0.8){
					matchThat=1;
					FATjet_DoubleSV[1]=ADDjet_DoubleSV[i];
				}
				if(matchThis&& matchThat){
					//cout<<"match"<<FATjet_DoubleSV[0]<<","<<FATjet_DoubleSV[0]<<endl;
					break;
				}
				
			}
			
			
			
			Float_t  ptHat = data.GetFloat("ptHat");
			
			
			Float_t ntrue= data.GetFloat("pu_nTrueInt");
			double PU_weight[3]={1,1,1};
			
			if(nameRoot!=2){
				if(ntrue<51){
					PU_weight[0] = LumiWeights_central.weight(ntrue);
					PU_weight[1]= LumiWeights_up.weight(ntrue);
					PU_weight[2] = LumiWeights_down.weight(ntrue);
				}
				else {
					PU_weight[0] = LumiWeights_central.weight(50);
					PU_weight[1] = LumiWeights_up.weight(50);
					PU_weight[2]= LumiWeights_down.weight(50);
				}
			}
			
			Float_t mcWeight= data.GetFloat("mcWeight");
			
			//?
			SelectedEvent_bcno=1;
			SelectedEvent_time=1;
			
			SelectedEvent_btagsf_bcUp=1;
			SelectedEvent_btagsf_bcDown=1;
			SelectedEvent_btagsf_lUp=1;
			SelectedEvent_btagsf_lDown=1;
			//?
			
			SelectedEvent_btagsf=btaggingscaleFactor;
			SelectedEvent_npuTrue=(int)ntrue;
			Float_t HT= data.GetFloat("HT");
			SelectedEvent_minv_leading2hjets=mjj;
			SelectedEvent_minv_leading2hjets_subtr=mjjRed;
			SelectedEvent_ht=(int)HT;
			SelectedEvent_pt_leading2hjets=(*thisJet+*thatJet).Pt();
			SelectedEvent_eta_leading2hjets=(*thisJet+*thatJet).Eta();
			SelectedEvent_phi_leading2hjets=(*thisJet+*thatJet).Phi();
			SelectedEvent_y_leading2hjets=(*thisJet+*thatJet).Rapidity();
			SelectedEvent_htHat=ptHat;
			SelectedEvent_npv=nVtx;
			SelectedEvent_deta_leading2hjets=thisJet->Eta()-thatJet->Eta();
			SelectedEvent_runno=runId;
			SelectedEvent_lumisec=lumiSection;
			SelectedEvent_evtno=eventId;
			SelectedEvent_nsubjetsBTaggedCSVL=nbtag;
			SelectedEvent_evtwt=mcWeight;
			SelectedEvent_evtwtPV=PU_weight[0];
			SelectedEvent_evtwtPVHigh=PU_weight[1];
			SelectedEvent_evtwtPVLow=PU_weight[2];
			
			HJets_njets=2;
			
			Float_t*  FATjetCHadEF = data.GetPtrFloat("FATjetCHadEF");
			Float_t*  FATjetPhoEF = data.GetPtrFloat("FATjetPhoEF");
			Float_t*  FATjetNHadEF = data.GetPtrFloat("FATjetNHadEF");
			
			
			for(int i=0;i<2;i++){
				TLorentzVector* thisJetInLoop;
				thisJetInLoop= (TLorentzVector*)fatjetP4->At(i);
				HJets_Pt[i]=thisJetInLoop->Pt();
				HJets_Eta[i]=thisJetInLoop->Eta();
				HJets_Phi[i]=thisJetInLoop->Phi();
				HJets_Energy[i]=thisJetInLoop->E();
				HJets_Mass[i]=thisJetInLoop->M();
				HJets_tau1[i]=fatjetTau1[i];
				HJets_tau2[i]=fatjetTau2[i];
				HJets_tau3[i]=fatjetTau3[i];
				HJets_hadFlavourSubjet0[i]=FATsubjetSDHadronFlavor[i][0];
				HJets_hadFlavourSubjet1[i]=FATsubjetSDHadronFlavor[i][1];
				HJets_ptSubjet0[i]=subjetP4[i][0]->Pt();
				HJets_ptSubjet1[i]=subjetP4[i][1]->Pt();
				HJets_etaSubjet0[i]=subjetP4[i][0]->Eta();
				HJets_etaSubjet1[i]=subjetP4[i][1]->Eta();
				HJets_phiSubjet0[i]=subjetP4[i][0]->Phi();
				HJets_phiSubjet1[i]=subjetP4[i][1]->Phi();
				HJets_energySubjet0[i]=subjetP4[i][0]->E();
				HJets_energySubjet1[i]=subjetP4[i][1]->E();
				HJets_CSVIVFv2[i]=FATjetCISVV2[i];
				HJets_MassPruned[i]=FATjetPRmass[i];
				HJets_MassSoftDrop[i]=FATjetSDmass[i];
				HJets_csvSubjet0[i]=subjetSDCSV[i][0];
				HJets_csvSubjet1[i]=subjetSDCSV[i][1];
				HJets_groomedMassCorr[i]=fatjetPRmassL2L3Corr[i];
				HJets_muf[i]=FATjetMuEF[i];
				HJets_emf[i]=FATjetCEmEF[i];
				HJets_nhf[i]=FATjetNHadEF[i];
				HJets_phf[i]=FATjetPhoEF[i];
				HJets_chf[i]=FATjetCHadEF[i];
				HJets_hadFlavour[i]=FATjetHadronFlavor[i];
				HJets_nsubjets[i]=FATnSubSDJet[i];
				HJets_doublesv[i]=FATjet_DoubleSV[i];
				if(HJets_csvSubjet0[i]>0.46)HJets_nsubjetsBTaggedCSVL[i]++;
				if(HJets_csvSubjet1[i]>0.46)HJets_nsubjetsBTaggedCSVL[i]++;
				
				
				//?
				HJets_Index[i]=1;
				HJets_ParentIndex[i]=1;
			}
			
			
			hh4b ->Fill();
			
			
		
		}//end event loop----------------------------------------------------------------------------------------
		hh4b->Write();
		outFile->Close();
	}	//end ntuple loop----------------------------------------------------------------------------------------
	
	
	
}

void miniTree(){
	miniTreeBase(1,2,"/home/neutron/HHbbbbAnalyzer/80x/mass/Bulk1200.root","miniTreeTest");
	
}


  