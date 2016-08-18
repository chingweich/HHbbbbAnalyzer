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


void miniTreeBase(int wMs,int wM, string st,string st2,int option=0) {
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0;
	
	
	
	//not yet added
	/*
	Double_t SelectedEvent_minv_leading2hjets;        hh4b->Branch("SelectedEvent_minv_leading2hjets",&SelectedEvent_minv_leading2hjets,"SelectedEvent_minv_leading2hjets/D");
	Double_t SelectedEvent_minv_leading2hjets_subtr;  hh4b->Branch("SelectedEvent_minv_leading2hjets_subtr",&SelectedEvent_minv_leading2hjets_subtr,"SelectedEvent_minv_leading2hjets_subtr/D");
	Double_t SelectedEvent_btagsf;                    hh4b->Branch("SelectedEvent_btagsf",&SelectedEvent_btagsf,"SelectedEvent_btagsf/D");
	Double_t SelectedEvent_btagsf_bcUp;               hh4b->Branch("SelectedEvent_btagsf_bcUp",&SelectedEvent_btagsf_bcUp,"SelectedEvent_btagsf_bcUp/D");
	Double_t SelectedEvent_btagsf_bcDown;             hh4b->Branch("SelectedEvent_btagsf_bcDown",&SelectedEvent_btagsf_bcDown,"SelectedEvent_btagsf_bcDown/D");
	Double_t SelectedEvent_btagsf_lUp;                hh4b->Branch("SelectedEvent_btagsf_lUp",&SelectedEvent_btagsf_lUp,"SelectedEvent_btagsf_lUp/D");
	Double_t SelectedEvent_btagsf_lDown;              hh4b->Branch("SelectedEvent_btagsf_lDown",&SelectedEvent_btagsf_lDown,"SelectedEvent_btagsf_lDown/D");
	Double_t SelectedEvent_y_leading2hjets;           hh4b->Branch("SelectedEvent_y_leading2hjets",&SelectedEvent_y_leading2hjets,"SelectedEvent_y_leading2hjets/D");
	Int_t HJets_Index[HJets_njets];                   hh4b->Branch("HJets_Index",&HJets_Index,"HJets_Index[HJets_njets]/I");
	Int_t HJets_nconsts[HJets_njets];                 hh4b->Branch("HJets_nconsts",&HJets_nconsts,"HJets_nconsts[HJets_njets]/I");
	Int_t HJets_nsubjetsBTaggedCSVL[HJets_njets];     hh4b->Branch("HJets_nsubjetsBTaggedCSVL",&HJets_nsubjetsBTaggedCSVL,"HJets_nsubjetsBTaggedCSVL[HJets_njets]/I");
 
  SelectedEvent_lhewts_               : vector<pair<int,double> >
  SelectedEvent_lhewts.first          : pair<int,double>
  SelectedEvent_lhewts.second         : pair<int,double>
 
	*/
	
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
	
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
		
		TFile* outFile = new TFile(Form("%s.root",st2.data()),"recreate");
		TTree* hh4b=new TTree("hh4b","hh4b");
	
		Int_t SelectedEvent_runno;                        hh4b->Branch("SelectedEvent_runno",&SelectedEvent_runno,"SelectedEvent_runno/I");
		Int_t SelectedEvent_lumisec;                      hh4b->Branch("SelectedEvent_lumisec",&SelectedEvent_lumisec,"SelectedEvent_lumisec/I");
		Int_t SelectedEvent_evtno;                        hh4b->Branch("SelectedEvent_evtno",&SelectedEvent_evtno,"SelectedEvent_evtno/I");
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
		Double_t SelectedEvent_evtwt;                     hh4b->Branch("SelectedEvent_evtwt",&SelectedEvent_evtwt,"SelectedEvent_evtwt/D");
		Double_t SelectedEvent_evtwtPV;                   hh4b->Branch("SelectedEvent_evtwtPV",&SelectedEvent_evtwtPV,"SelectedEvent_evtwtPV/D");
		Double_t SelectedEvent_evtwtPVLow;                hh4b->Branch("SelectedEvent_evtwtPVLow",&SelectedEvent_evtwtPVLow,"SelectedEvent_evtwtPVLow/D");
		Double_t SelectedEvent_evtwtPVHigh;               hh4b->Branch("SelectedEvent_evtwtPVHigh",&SelectedEvent_evtwtPVHigh,"SelectedEvent_evtwtPVHigh/D");
	
		Int_t HJets_nsubjets[HJets_njets];                hh4b->Branch("HJets_nsubjets",&HJets_nsubjets,"HJets_nsubjets[HJets_njets]/I");
	
		Double_t HJets_Pt[HJets_njets];                   hh4b->Branch("HJets_Pt",&HJets_Pt,"HJets_Pt[HJets_njets]/D");
		Double_t HJets_Eta[HJets_njets];                  hh4b->Branch("HJets_Eta",&HJets_Eta,"HJets_Eta[HJets_njets]/D");
		Double_t HJets_Phi[HJets_njets];                  hh4b->Branch("HJets_Phi",&HJets_Phi,"HJets_Phi[HJets_njets]/D");
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
			TLorentzVector* subjetP4[2][2];
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					subjetP4[i][j]=new TLorentzVector(0,0,0,0);
					subjetP4[i][j]->SetPxPyPzE(subjetSDPx[i][j],subjetSDPy[i][j],subjetSDPz[i][j],subjetSDE[i][j]);
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
			
			SelectedEvent_npuTrue=(int)ntrue;
			Float_t HT= data.GetFloat("HT");
			SelectedEvent_ht=(int)HT;
			SelectedEvent_pt_leading2hjets=(*thisJet+*thatJet).Pt();
			SelectedEvent_eta_leading2hjets=(*thisJet+*thatJet).Eta();
			SelectedEvent_phi_leading2hjets=(*thisJet+*thatJet).Phi();
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


  