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
#include "../untuplizer.h"
//#include "ApplyJER.h"
//#include "jetEnergyScale.h"

#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"

float getPUPPIweight(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(1.00626, -1.06161, 0.07999,1.20454 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters( 1.05807,-5.91971e-05,2.296e-07, -1.98795e-10,6.67382e-14,  -7.80604e-18);
  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters( 1.26638,-0.000658496,  9.73779e-07,-5.93843e-10, 1.61619e-13, -1.6272e-17);
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  totalWeight = genCorr * recoCorr;
  return totalWeight;
}


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
	/*
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
	*/
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
		Double_t LeadingDiJets_deta_leading2hjets;        hh4b->Branch("LeadingDiJets_deta_leading2hjets",&LeadingDiJets_deta_leading2hjets,"LeadingDiJets_deta_leading2hjets/D");
		Double_t LeadingDiJets_pt_leading2hjets;          hh4b->Branch("LeadingDiJets_pt_leading2hjets",&LeadingDiJets_pt_leading2hjets,"LeadingDiJets_pt_leading2hjets/D");
		Double_t LeadingDiJets_eta_leading2hjets;         hh4b->Branch("LeadingDiJets_eta_leading2hjets",&LeadingDiJets_eta_leading2hjets,"LeadingDiJets_eta_leading2hjets/D");
		Double_t LeadingDiJets_phi_leading2hjets;         hh4b->Branch("LeadingDiJets_phi_leading2hjets",&LeadingDiJets_phi_leading2hjets,"LeadingDiJets_phi_leading2hjets/D");
		Double_t LeadingDiJets_y_leading2hjets;           hh4b->Branch("LeadingDiJets_y_leading2hjets",&LeadingDiJets_y_leading2hjets,"LeadingDiJets_y_leading2hjets/D");
		Double_t LeadingDiJets_energy_leading2hjets;    hh4b->Branch("LeadingDiJets_energy_leading2hjets",&LeadingDiJets_energy_leading2hjets,"LeadingDiJets_energy_leading2hjets/D");
		Double_t LeadingDiJets_minv_leading2hjets;        hh4b->Branch("LeadingDiJets_minv_leading2hjets",&LeadingDiJets_minv_leading2hjets,"LeadingDiJets_minv_leading2hjets/D");
		Double_t LeadingDiJets_minv_leading2hjets_subtr;  hh4b->Branch("LeadingDiJets_minv_leading2hjets_subtr",&LeadingDiJets_minv_leading2hjets_subtr,"LeadingDiJets_minv_leading2hjets_subtr/D");
		Int_t LeadingDiJets_doublebcategory;                   hh4b->Branch("LeadingDiJets_doublebcategory",&LeadingDiJets_doublebcategory,"LeadingDiJets_doublebcategory/I");
	
	
		Double_t SelectedEvent_evtwt;                     hh4b->Branch("SelectedEvent_evtwt",&SelectedEvent_evtwt,"SelectedEvent_evtwt/D");
		Double_t SelectedEvent_evtwtPV;                   hh4b->Branch("SelectedEvent_evtwtPV",&SelectedEvent_evtwtPV,"SelectedEvent_evtwtPV/D");
		Double_t SelectedEvent_evtwtPVLow;                hh4b->Branch("SelectedEvent_evtwtPVLow",&SelectedEvent_evtwtPVLow,"SelectedEvent_evtwtPVLow/D");
		Double_t SelectedEvent_evtwtPVHigh;               hh4b->Branch("SelectedEvent_evtwtPVHigh",&SelectedEvent_evtwtPVHigh,"SelectedEvent_evtwtPVHigh/D");
		
		
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
		
		Double_t AK8Jets_Pt[HJets_njets];                    hh4b->Branch("AK8Jets_Pt",&AK8Jets_Pt,"AK8Jets_Pt[HJets_njets]/D");
		Double_t AK8Jets_Eta[HJets_njets];                  hh4b->Branch("AK8Jets_Eta",&AK8Jets_Eta,"AK8Jets_Eta[HJets_njets]/D");
		Double_t AK8Jets_Phi[HJets_njets];                  hh4b->Branch("AK8Jets_Phi",&AK8Jets_Phi,"AK8Jets_Phi[HJets_njets]/D");
		Double_t AK8Jets_Energy[HJets_njets];           hh4b->Branch("AK8Jets_Energy",&AK8Jets_Energy,"AK8Jets_Energy[HJets_njets]/D");
		Double_t AK8Jets_Mass[HJets_njets];                 hh4b->Branch("AK8Jets_Mass",&AK8Jets_Mass,"AK8Jets_Mass[HJets_njets]/D");
		Double_t AK8Jets_MassPruned[HJets_njets];           hh4b->Branch("AK8Jets_MassPruned",&AK8Jets_MassPruned,"AK8Jets_MassPruned[HJets_njets]/D");
		Double_t AK8Jets_MassSoftDrop[HJets_njets];         hh4b->Branch("AK8Jets_MassSoftDrop",&AK8Jets_MassSoftDrop,"AK8Jets_MassSoftDrop[HJets_njets]/D");
		Double_t AK8Jets_tau1[HJets_njets];                 hh4b->Branch("AK8Jets_tau1",&AK8Jets_tau1,"AK8Jets_tau1[HJets_njets]/D");
		Double_t AK8Jets_tau2[HJets_njets];                 hh4b->Branch("AK8Jets_tau2",&AK8Jets_tau2,"AK8Jets_tau2[HJets_njets]/D");
		Double_t AK8Jets_tau3[HJets_njets];                 hh4b->Branch("AK8Jets_tau3",&AK8Jets_tau3,"AK8Jets_tau3[HJets_njets]/D");
		Double_t AK8Jets_ptSubjet0[HJets_njets];            hh4b->Branch("AK8Jets_ptSubjet0",&AK8Jets_ptSubjet0,"AK8Jets_ptSubjet0[HJets_njets]/D");
		Double_t AK8Jets_ptSubjet1[HJets_njets];            hh4b->Branch("AK8Jets_ptSubjet1",&AK8Jets_ptSubjet1,"AK8Jets_ptSubjet1[HJets_njets]/D");
		Double_t AK8Jets_etaSubjet0[HJets_njets];           hh4b->Branch("AK8Jets_etaSubjet0",&AK8Jets_etaSubjet0,"AK8Jets_etaSubjet0[HJets_njets]/D");
		Double_t AK8Jets_etaSubjet1[HJets_njets];           hh4b->Branch("AK8Jets_etaSubjet1",&AK8Jets_etaSubjet1,"AK8Jets_etaSubjet1[HJets_njets]/D");
		Double_t AK8Jets_phiSubjet0[HJets_njets];           hh4b->Branch("AK8Jets_phiSubjet0",&AK8Jets_phiSubjet0,"AK8Jets_phiSubjet0[HJets_njets]/D");
		Double_t AK8Jets_phiSubjet1[HJets_njets];           hh4b->Branch("AK8Jets_phiSubjet1",&AK8Jets_phiSubjet1,"AK8Jets_phiSubjet1[HJets_njets]/D");
		Double_t AK8Jets_energySubjet0[HJets_njets];        hh4b->Branch("AK8Jets_energySubjet0",&AK8Jets_energySubjet0,"AK8Jets_energySubjet0[HJets_njets]/D");
		Double_t AK8Jets_energySubjet1[HJets_njets];        hh4b->Branch("AK8Jets_energySubjet1",&AK8Jets_energySubjet1,"AK8Jets_energySubjet1[HJets_njets]/D");
		Double_t AK8Jets_hadFlavour[HJets_njets];           hh4b->Branch("AK8Jets_hadFlavour",&AK8Jets_hadFlavour,"AK8Jets_hadFlavour[HJets_njets]/D");
		Double_t AK8Jets_hadFlavourSubjet0[HJets_njets];    hh4b->Branch("AK8Jets_hadFlavourSubjet0",&AK8Jets_hadFlavourSubjet0,"AK8Jets_hadFlavourSubjet0[HJets_njets]/D");
		Double_t AK8Jets_hadFlavourSubjet1[HJets_njets];    hh4b->Branch("AK8Jets_hadFlavourSubjet1",&AK8Jets_hadFlavourSubjet1,"AK8Jets_hadFlavourSubjet1[HJets_njets]/D");
		Double_t AK8Jets_csvSubjet0[HJets_njets];           hh4b->Branch("AK8Jets_csvSubjet0",&AK8Jets_csvSubjet0,"AK8Jets_csvSubjet0[HJets_njets]/D");
		Double_t AK8Jets_csvSubjet1[HJets_njets];           hh4b->Branch("AK8Jets_csvSubjet1",&AK8Jets_csvSubjet1,"AK8Jets_csvSubjet1[HJets_njets]/D");
		
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			bool isData  = data.GetBool("isData"); // only apply trigger if it's data

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];

        if( (thisTrig.find("HLT_PFHT800_v")!= std::string::npos && results==1) ||
	    (thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1) || // for runH
	    (thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos && results==1) ||
	    (thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos && results==1) ||
	    (thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos && results==1) ||
	    (thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos && results==1) 
            )
          {
            passTrigger=true;
            break;
          }
      }

    if(!passTrigger && isData)continue;
    nPass[0]++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[1]++;

    //1. veto events containing 1 or 2 leptons

    Int_t nEle         = data.GetInt("nEle");
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
    Float_t* eleSCEta                = data.GetPtrFloat("eleScEta");
    Float_t* eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
    Float_t* eleHoverE               = data.GetPtrFloat("eleHoverE");
    Float_t* eleEcalPFClusterIso     = data.GetPtrFloat("eleEcalPFClusterIso");
    Float_t* eleHcalPFClusterIso     = data.GetPtrFloat("eleHcalPFClusterIso");
    Float_t* eleDr03TkSumPt          = data.GetPtrFloat("eleDr03TkSumPt");
    Float_t* eledEtaAtVtx            = data.GetPtrFloat("eledEtaAtVtx");
    Float_t* eledPhiAtVtx            = data.GetPtrFloat("eledPhiAtVtx");
    Int_t*   eleCharge               = data.GetPtrInt("eleCharge");
    vector<bool>    &isMediumEle= *((vector<bool>*) data.GetPtr("eleIsPassMVAMedium"));
    vector<bool>     &isTightEle= *((vector<bool>*) data.GetPtr("eleIsPassMVATight"));

    // check if there is any single tight electron
    bool hasSingleTightElectron=false;

    for(int ie=0; ie < nEle; ie++){

      TLorentzVector* thisEle= (TLorentzVector*)eleP4->At(ie);
      if( fabs(thisEle->Eta())>2.4 )continue;
      if( thisEle->Pt() < 25 )continue;
      bool preSelect= 

	( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
	  eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||

	( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
	  eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
    
      if(!preSelect)continue;

      if( !isTightEle[ie] )continue;

	
      hasSingleTightElectron=true;
      break;
      
    }
    if( hasSingleTightElectron )continue;
    nPass[2]++;


    // check if there is any double loose electron
    bool hasLooseElectronPair=false;
    for(int ie=0; ie < nEle; ie++){
      for(int je=0; je < ie; je++){

	if(eleCharge[ie]*eleCharge[je]>0)continue;

   	TLorentzVector* thisEle= (TLorentzVector*)eleP4->At(ie);

    	if( fabs(thisEle->Eta())>2.4 )continue;
    	if( thisEle->Pt() < 15 )continue;

   	TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(je);
    	if( fabs(thatEle->Eta())>2.4 )continue;
    	if( thatEle->Pt() < 15 )continue;

	bool preSelect_ie= 

	( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
	  eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||

	( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
	  eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
	
	if(!preSelect_ie)continue;


	bool preSelect_je= 

	( fabs(eleSCEta[je]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[je] < 0.012 && 
	  eleHoverE[je] < 0.09 && (eleEcalPFClusterIso[je]/thatEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[je]/thatEle->Pt()) < 0.25 && (eleDr03TkSumPt[je]/thatEle->Pt()) < 0.18 
	  && fabs(eledEtaAtVtx[je]) < 0.0095 && fabs(eledPhiAtVtx[je]) < 0.065 ) ||

	( fabs(eleSCEta[je]) > 1.5660 && fabs(eleSCEta[je]) < 2.5 && eleSigmaIEtaIEtaFull5x5[je] < 0.033 && 
	  eleHoverE[je] <0.09 && (eleEcalPFClusterIso[je]/thatEle->Pt()) < 0.45 && (eleHcalPFClusterIso[je]/thatEle->Pt()) < 0.28 && (eleDr03TkSumPt[je]/thatEle->Pt()) < 0.18 );
	
	if(!preSelect_je)continue;


	if( !isMediumEle[ie] )continue;
	if( !isMediumEle[je] )continue;
	
	hasLooseElectronPair=true;
	break;
      
      }
    }

    if( hasLooseElectronPair )continue;
    nPass[3]++;
    
    
    // veto of single tight muon or two oppositely-charged muon pair
    Int_t nMu          = data.GetInt("nMu");
    TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
    Int_t*   muCharge        = data.GetPtrInt("muCharge");
    vector<bool>    &isLooseMuon= *((vector<bool>*) data.GetPtr("isLooseMuon"));
    vector<bool>    &isTightMuon= *((vector<bool>*) data.GetPtr("isTightMuon"));
    Float_t*          muIso1  = data.GetPtrFloat("muChHadIso");
    Float_t*          muIso2  = data.GetPtrFloat("muNeHadIso");
    Float_t*          muIso3  = data.GetPtrFloat("muGamIso");
    Float_t*          muIso4  = data.GetPtrFloat("muPUPt");
    
    // check if there is any single tight muon
    bool hasSingleTightMuon=false;

    for(int im=0; im < nMu; im++){

      TLorentzVector* thisMu = (TLorentzVector*)muP4->At(im);

      if( fabs(thisMu->Eta())>2.4 )continue;
      if( thisMu->Pt() < 25 )continue;
      if( !isTightMuon[im] )continue;

      float iso =  (muIso1[im] + max(0., muIso2[im]+muIso3[im] - 0.5*muIso4[im]))/thisMu->Pt();
      if( iso > 0.15)continue;
	
      hasSingleTightMuon=true;
      break;
      
    }
    if( hasSingleTightMuon )continue;
    nPass[4]++;


    // check if there is any double loose muons
    bool hasLooseMuonPair=false;
    for(int im=0; im < nMu; im++){
      for(int jm=0; jm < im; jm++){

	if(muCharge[im]*muCharge[jm]>0)continue;
	if( !isLooseMuon[im] )continue;
	if( !isLooseMuon[jm] )continue;

   	TLorentzVector* thisMu = (TLorentzVector*)muP4->At(im);

    	if( fabs(thisMu->Eta())>2.4 )continue;
    	if( thisMu->Pt() < 15 )continue;

   	TLorentzVector* thatMu = (TLorentzVector*)muP4->At(jm);
    	if( fabs(thatMu->Eta())>2.4 )continue;
    	if( thatMu->Pt() < 15 )continue;

	float iso1=  (muIso1[im] + max(0., muIso2[im]+muIso3[im] - 0.5*muIso4[im]))/thisMu->Pt();
	if( iso1> 0.25)continue;

	float iso2=  (muIso1[jm] + max(0., muIso2[jm]+muIso3[jm] - 0.5*muIso4[jm]))/thatMu->Pt();
	if( iso2> 0.25)continue;
	
	hasLooseMuonPair=true;
	break;
      
      }
    }
    if( hasLooseMuonPair )continue;

    nPass[5]++;

    // apply jet requirement

    int nJets         = data.GetInt("FATnJet");
    TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");

    Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
    Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
    Float_t*  FATjetPuppiTau3 = data.GetPtrFloat("FATjetPuppiTau3");
    Float_t*  fatjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
   
    vector<bool>    &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));

    Float_t*  fatjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
    Float_t*  fatjetMuoEF = data.GetPtrFloat("FATjetMuoEF");

    int nADDJets         = data.GetInt("ADDnJet");
    TClonesArray* addjetP4   = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doublesv = data.GetPtrFloat("ADDjet_DoubleSV");

    if(nJets<2)continue;
    nPass[6]++;

    if(nADDJets<2)continue;
    nPass[7]++;

    Int_t nGoodJets=0;
    const float dRMax=0.8;
    int addJetIndex[2]={-1,-1};
    
    float thea_mass[2]={-1,-1};

    for(int ij=0; ij<2; ij++)
      {
	if( !passFatJetTightID[ij] )continue;
	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if( thisJet->Pt()<300 )continue;
	if( fabs(thisJet->Eta())>2.4 )continue;
	float tau21_puppi = FATjetPuppiTau2[ij]/FATjetPuppiTau1[ij];
	if( tau21_puppi > 0.6 )continue;

	TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
	float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
	float raw_mass = fatjetPuppiSDmass[ij];
	thea_mass[ij] = raw_mass*thea_corr;

	if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;


	// now look for add jets

	for(int k=0; k < nADDJets; k++){
	  TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
	  if(thisJet->DeltaR(*thatJet)<dRMax && addJetIndex[ij]<0)
	    {
	      addJetIndex[ij]=k;
	      break;
	    }
	}

	nGoodJets++;
      } // end of loop over two leading jets


    if(nGoodJets<2)continue;
    nPass[8]++;

    if(addJetIndex[0]==addJetIndex[1])continue;
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    nPass[9]++;

    // check if any jet fail loose double b-tag

    const float DBLoose = 0.3;
    const float DBTight = 0.8;

    if(addjet_doublesv[addJetIndex[0]]< DBLoose)continue;
    if(addjet_doublesv[addJetIndex[1]]< DBLoose)continue;


    nPass[10]++;
    TLorentzVector* higgsJet[2];
    for(int i=0;i<2;i++)higgsJet[i] = (TLorentzVector*)fatjetP4->At(i);
    float dEta = fabs(higgsJet[0]->Eta()-higgsJet[1]->Eta());

    if(dEta>1.3)continue;
    nPass[11]++;
    float mjj = (*higgsJet[0]+*higgsJet[1]).M();

    float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
    if(mjjred<750)continue;
    nPass[12]++;
			
			
			Long64_t runId        = data.GetLong64("runId");
			Long64_t lumiSection        = data.GetLong64("lumiSection");
			Long64_t eventId        = data.GetLong64("eventId");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nJets);
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nJets);
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nJets);
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nJets);
			//vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			Float_t*  FATjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
			Int_t*  FATjetHadronFlavor = data.GetPtrInt("FATjetHadronFlavor");
			Int_t*  FATnSubSDJet = data.GetPtrInt("FATnSubSDJet");
			Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuoEF = data.GetPtrFloat("FATjetMuoEF");
			
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
			int nbtag=0;
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			if(subjetSDCSV[0][0]>0.46)nbtag++;
			if(subjetSDCSV[0][1]>0.46)nbtag++;
			if(subjetSDCSV[1][0]>0.46)nbtag++;
			if(subjetSDCSV[1][1]>0.46)nbtag++;
			
			
			int nFail=0,nLoose=0,nMed=0,nTight=0;
			if(addjet_doublesv[addJetIndex[0]]<0.3)nFail++;
			else if (addjet_doublesv[addJetIndex[0]]<0.6)nLoose++;
			else if (addjet_doublesv[addJetIndex[0]]<0.8)nMed++;
			else nTight++;
			
			if(addjet_doublesv[addJetIndex[1]]<0.3)nFail++;
			else if (addjet_doublesv[addJetIndex[1]]<0.6)nLoose++;
			else if (addjet_doublesv[addJetIndex[1]]<0.8)nMed++;
			else nTight++;
			
			
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
			LeadingDiJets_doublebcategory=-1;
			if(nFail==2)LeadingDiJets_doublebcategory=0;
			else if (nLoose==2)LeadingDiJets_doublebcategory=2;
			else if (nMed==2)LeadingDiJets_doublebcategory=4;
			else if (nTight==2)LeadingDiJets_doublebcategory=8;
			else if (nLoose+nMed==2)LeadingDiJets_doublebcategory=6;
			else if (nFail+nLoose==2)LeadingDiJets_doublebcategory=3;
			else if (nFail+nMed==2)LeadingDiJets_doublebcategory=5;
			SelectedEvent_btagsf=btaggingscaleFactor;
			SelectedEvent_npuTrue=(int)ntrue;
			Float_t HT= data.GetFloat("HT");
			LeadingDiJets_minv_leading2hjets=mjj;
			LeadingDiJets_minv_leading2hjets_subtr=mjjred;
			SelectedEvent_ht=(int)HT;
			LeadingDiJets_pt_leading2hjets=(*higgsJet[0]+*higgsJet[1]).Pt();
			LeadingDiJets_eta_leading2hjets=(*higgsJet[0]+*higgsJet[1]).Eta();
			LeadingDiJets_phi_leading2hjets=(*higgsJet[0]+*higgsJet[1]).Phi();
			LeadingDiJets_energy_leading2hjets=(*higgsJet[0]+*higgsJet[1]).E();
			LeadingDiJets_y_leading2hjets=(*higgsJet[0]+*higgsJet[1]).Rapidity();
			SelectedEvent_htHat=ptHat;
			SelectedEvent_npv=nVtx;
			LeadingDiJets_deta_leading2hjets=dEta;
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
				HJets_tau1[i]=FATjetPuppiTau1[i];
				HJets_tau2[i]=FATjetPuppiTau2[i];
				HJets_tau3[i]=FATjetPuppiTau3[i];
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
				HJets_muf[i]=FATjetMuoEF[i];
				HJets_emf[i]=FATjetCEmEF[i];
				HJets_nhf[i]=FATjetNHadEF[i];
				HJets_phf[i]=FATjetPhoEF[i];
				HJets_chf[i]=FATjetCHadEF[i];
				HJets_hadFlavour[i]=FATjetHadronFlavor[i];
				HJets_nsubjets[i]=FATnSubSDJet[i];
				HJets_doublesv[i]=addjet_doublesv[addJetIndex[i]];
				if(HJets_csvSubjet0[i]>0.46)HJets_nsubjetsBTaggedCSVL[i]++;
				if(HJets_csvSubjet1[i]>0.46)HJets_nsubjetsBTaggedCSVL[i]++;
				
				
				//?
				HJets_Index[i]=1;
				
			}
			
			
			hh4b ->Fill();
			
			
		
		}//end event loop----------------------------------------------------------------------------------------
		hh4b->Write();
		outFile->Close();
	}	//end ntuple loop----------------------------------------------------------------------------------------
	
	
	
}

void miniTree_v2(){
	miniTreeBase(1,2,"BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root","miniTreeTest");
	
}


  