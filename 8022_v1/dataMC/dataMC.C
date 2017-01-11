//Histogram 
#include <TH1.h>
//#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
//#include <TGraph.h>
//#include <TGraphErrors.h>

//vector ,string ,stream
#include <vector>
#include <string>
#include <iostream>
//#include <fstream>
//#include <sstream>

//root feature
//#include <TLegend.h>
//#include <TRandom.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TFile.h>
//#include <TCanvas.h>
//#include "TSystem.h"
//#include "TStyle.h"
#include <TClonesArray.h>
#include <TF1.h>
//math 
#include <cmath>
//#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../../untuplizer.h"
//#include "jetEnergyScale.h"

//#include "getPUPPIweight.h"
#include "standalone_LumiReWeighting.cc"
#include "standalone_LumiReWeighting.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"
//#include "HHbbbbBtagMakeEff_76.C"
//#include "HH4bBtagEffBase_80_v2.C"

using namespace std;

//data MC comparison
float getPUPPIweight(float puppipt, float puppieta ){
  TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(1.00626, -1.06161, 0.07999,1.20454 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters( 1.05807,-5.91971e-05,2.296e-07, -1.98795e-10,6.67382e-14,  -7.80604e-18);
  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters( 1.26638,-0.000658496,  9.73779e-07,-5.93843e-10, 1.61619e-13, -1.6272e-17);
  float genCorr  = 1.,recoCorr = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  return genCorr * recoCorr;
}

void skimTree(int w , string st,string st2){
	double nPass[30]={0},total=0,fixScaleNum[2]={0};

	vector<string> str;
	vector<double> up;vector<double> down;vector<double> nbins;
	str.push_back("pt"); up.push_back(2000);  down.push_back(0);nbins.push_back(200);
	//str.push_back("pt_j1"); up.push_back(2000);  down.push_back(0);nbins.push_back(200);
	str.push_back("eta"); up.push_back(3);  down.push_back(-3);nbins.push_back(60);
	//str.push_back("eta_j1"); up.push_back(3);  down.push_back(-3);nbins.push_back(60);
	str.push_back("puppiSDMassThea"); up.push_back(200);  down.push_back(40);nbins.push_back(16);
	//str.push_back("puppiSDMassThea_j1"); up.push_back(200);  down.push_back(40);nbins.push_back(16);
	str.push_back("doubleSV"); up.push_back(1);  down.push_back(-1);nbins.push_back(40);
	//str.push_back("doubleSV_j1"); up.push_back(1);  down.push_back(-1);nbins.push_back(40);
	str.push_back("puppiTau21"); up.push_back(1);  down.push_back(0);nbins.push_back(25);
	//str.push_back("puppiTau21_j1"); up.push_back(1);  down.push_back(0);nbins.push_back(25);
	
	str.push_back("nMu"); up.push_back(6);  down.push_back(0);nbins.push_back(6);
	str.push_back("nEle"); up.push_back(6);  down.push_back(0);nbins.push_back(6);
	str.push_back("nCh"); up.push_back(100);  down.push_back(0);nbins.push_back(25);
	str.push_back("MuEnFr"); up.push_back(1);  down.push_back(0);nbins.push_back(50);
	str.push_back("EmEnFr"); up.push_back(1);  down.push_back(0);nbins.push_back(50);
	str.push_back("nSv"); up.push_back(6);  down.push_back(0);nbins.push_back(6);
	str.push_back("svMass0"); up.push_back(10);  down.push_back(0);nbins.push_back(25);
	str.push_back("svMass1"); up.push_back(10);  down.push_back(0);nbins.push_back(25);
	
	str.push_back("deltaEta"); up.push_back(2);  down.push_back(0);nbins.push_back(40);
	str.push_back("totalMass"); up.push_back(5000);  down.push_back(0);nbins.push_back(250);
	str.push_back("totalMassRed"); up.push_back(5000);  down.push_back(0);nbins.push_back(250);
	
	
	str.push_back("nVtx"); up.push_back(60);  down.push_back(0);nbins.push_back(60);
	str.push_back("nVtxWeighted"); up.push_back(60);  down.push_back(0);nbins.push_back(60);
	
	
	
	TH1D* th1[35],* th_flavor[4][35];
	for(int i=0;i<13;i++ ){
		th1[i*2]=new TH1D(Form("%s_j0",str[i].data()),Form("%s_j0",str[i].data()),nbins[i],down[i],up[i]);
		th1[i*2+1]=new TH1D(Form("%s_j1",str[i].data()),Form("%s_j1",str[i].data()),nbins[i],down[i],up[i]);
	}
	for(int i=26;i<31;i++)th1[i]=new TH1D(Form("%s",str[i-13].data()),Form("%s",str[i-13].data()),nbins[i-13],down[i-13],up[i-13]);
	
	for(int i=0;i<31;i++){
		th1[i]->Sumw2();
		th_flavor[0][i]=(TH1D* )th1[i]->Clone(Form("%s_bb",th1[i]->GetTitle()));
		th_flavor[1][i]=(TH1D* )th1[i]->Clone(Form("%s_b",th1[i]->GetTitle()));
		th_flavor[2][i]=(TH1D* )th1[i]->Clone(Form("%s_cc",th1[i]->GetTitle()));
		th_flavor[3][i]=(TH1D* )th1[i]->Clone(Form("%s_udcsg",th1[i]->GetTitle()));
		//for(int j=0;j<4;j++)th_flavor[j][i]->Sumw2();
	}
	
	standalone_LumiReWeighting LumiWeights_central(0);
	TFile* f= TFile::Open(Form("%s%d.root",st.data(),w));
	TTree* tree=(TTree*)f->Get("treeMaker");
	TreeReader data(tree);
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
		if(jEntry%1000==0)cout<<jEntry<<" out of "<<data.GetEntriesFast() <<" events are processed"<<endl;
		data.GetEntry(jEntry);
		bool isData  = data.GetBool("isData"); // only apply trigger if it's data
		
		
		Float_t ntrue= data.GetFloat("pu_nTrueInt");
		double PU_weight=1;
		if(!isData){
			if(ntrue<51)PU_weight= LumiWeights_central.weight(ntrue);
			else PU_weight=  LumiWeights_central.weight(52);
			//if(PU_weight>2)cout<<jEntry<<","<<PU_weight<<endl;
		}
		fixScaleNum[0]+=PU_weight;
		
		std::string* trigName = data.GetPtrString("hlt_trigName");
		vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
		const Int_t nsize = data.GetPtrStringSize();

		bool passTrigger=false;
		for(int it=0; it< nsize; it++){
			std::string thisTrig= trigName[it];
			bool results = trigResult[it];
			if( (thisTrig.find("HLT_PFHT800_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1) || // for runH
				(thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos && results==1) ||
				(thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos && results==1)){
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
		/*
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
			bool preSelect= ( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
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

				bool preSelect_ie= ( fabs(eleSCEta[ie]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[ie] < 0.012 && 
									eleHoverE[ie] < 0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.37 && ( eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.25 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 
									&& fabs(eledEtaAtVtx[ie]) < 0.0095 && fabs(eledPhiAtVtx[ie]) < 0.065 ) ||
									( fabs(eleSCEta[ie]) > 1.5660 && fabs(eleSCEta[ie]) < 2.5 && eleSigmaIEtaIEtaFull5x5[ie] < 0.033 && 
									eleHoverE[ie] <0.09 && (eleEcalPFClusterIso[ie]/thisEle->Pt()) < 0.45 && (eleHcalPFClusterIso[ie]/thisEle->Pt()) < 0.28 && (eleDr03TkSumPt[ie]/thisEle->Pt()) < 0.18 );
	
				if(!preSelect_ie)continue;
				bool preSelect_je= ( fabs(eleSCEta[je]) < 1.4442 && eleSigmaIEtaIEtaFull5x5[je] < 0.012 && 
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
*/
		// apply jet requirement

		int nJets         = data.GetInt("FATnJet");
		TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
		TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");

		Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
		Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
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

		if( !passFatJetTightID[0] ||  !passFatJetTightID[1] )continue;
		TLorentzVector* higgsJet[2];
		for(int i=0;i<2;i++)higgsJet[i] = (TLorentzVector*)fatjetP4->At(i);
		if( higgsJet[0]->Pt()<300 || higgsJet[1]->Pt()<300)continue;
		if( fabs(higgsJet[0]->Eta())>2.4 || fabs(higgsJet[1]->Eta())>2.4)continue;
		float dEta = fabs(higgsJet[0]->Eta()-higgsJet[1]->Eta());
		if(dEta>1.3)continue;
		
		Int_t nGoodJets=0;
		const float dRMax=0.8;
		int addJetIndex[2]={-1,-1};
		float thea_mass[2]={-1,-1};
		for(int ij=0; ij<2; ij++){
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
			float tau21_puppi = FATjetPuppiTau2[ij]/FATjetPuppiTau1[ij];
			if( tau21_puppi > 0.6 )continue;
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
			float raw_mass = fatjetPuppiSDmass[ij];
			thea_mass[ij] = raw_mass*thea_corr;

			//if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;
			if(thea_mass[ij] < 50)continue;
			// now look for add jets
			for(int k=0; k < nADDJets; k++){
				TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
				if(thisJet->DeltaR(*thatJet)<dRMax && addJetIndex[ij]<0){
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

		//if(addjet_doublesv[addJetIndex[0]]< DBLoose)continue;
		//if(addjet_doublesv[addJetIndex[1]]< DBLoose)continue;


		nPass[10]++;
		
		
		nPass[11]++;
		float mjj = (*higgsJet[0]+*higgsJet[1]).M();

		float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
		if(mjjred<750)continue;
		nPass[12]++;
		
		int event_flavor=-1;
		Int_t* FATjetHadronFlavor        = data.GetPtrInt("FATjetHadronFlavor");
		vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
		if(FATjetHadronFlavor[0]==5 && FATsubjetSDHadronFlavor[0][0]==5 && FATsubjetSDHadronFlavor[0][1]==5)event_flavor=0;
		else if(FATjetHadronFlavor[1]==5 && FATsubjetSDHadronFlavor[1][0]==5 && FATsubjetSDHadronFlavor[1][1]==5)event_flavor=0;
		else if(FATjetHadronFlavor[0]==5 && (FATsubjetSDHadronFlavor[0][0]==5 || FATsubjetSDHadronFlavor[0][1]==5))event_flavor=1;
		else if(FATjetHadronFlavor[1]==5 && (FATsubjetSDHadronFlavor[1][0]==5 || FATsubjetSDHadronFlavor[1][1]==5))event_flavor=1;
		else if(FATjetHadronFlavor[0]==4 && FATsubjetSDHadronFlavor[0][0]==4 && FATsubjetSDHadronFlavor[0][1]==4)event_flavor=2;
		else if(FATjetHadronFlavor[1]==4 && FATsubjetSDHadronFlavor[1][0]==4 && FATsubjetSDHadronFlavor[1][1]==4)event_flavor=2;
		else if(FATjetHadronFlavor[0]==4 && (FATsubjetSDHadronFlavor[0][0]==4 || FATsubjetSDHadronFlavor[0][1]==4))event_flavor=2;
		else if(FATjetHadronFlavor[1]==4 && (FATsubjetSDHadronFlavor[1][0]==4 || FATsubjetSDHadronFlavor[1][1]==4))event_flavor=2;
		else event_flavor=3;
		
		Int_t* CMulti=data.GetPtrInt("FATjetCMulti");
		Int_t* EleMulti=data.GetPtrInt("FATjetEleMulti");
		Int_t* MuoMulti=data.GetPtrInt("FATjetMuoMulti");
		Float_t*  FATjetMuoEF = data.GetPtrFloat("FATjetMuoEF");
		Float_t*  EleEF = data.GetPtrFloat("FATjetEleEF");
		Int_t* nSV=data.GetPtrInt("ADDjet_nSV");
		vector<float>   * SVMass  =  data.GetPtrVectorFloat("ADDjet_SVMass");
		
		
		for(int i=0;i<2;i++){
			th_flavor[event_flavor][i]->Fill(higgsJet[i]->Pt(),PU_weight);
			th_flavor[event_flavor][2+i]->Fill(higgsJet[i]->Eta(),PU_weight);
			th_flavor[event_flavor][4+i]->Fill(thea_mass[i],PU_weight);
			th_flavor[event_flavor][6+i]->Fill(addjet_doublesv[addJetIndex[i]],PU_weight);
			th_flavor[event_flavor][8+i]->Fill(FATjetPuppiTau2[i]/FATjetPuppiTau1[i],PU_weight);
			th_flavor[event_flavor][10+i]->Fill(MuoMulti[i],PU_weight);
			th_flavor[event_flavor][12+i]->Fill(EleMulti[i],PU_weight);
			th_flavor[event_flavor][14+i]->Fill(CMulti[i],PU_weight);		
			th_flavor[event_flavor][16+i]->Fill(FATjetMuoEF[i],PU_weight);
			th_flavor[event_flavor][18+i]->Fill(EleEF[i],PU_weight);
			th_flavor[event_flavor][20+i]->Fill(nSV[addJetIndex[i]],PU_weight);
			if(nSV[addJetIndex[i]]>0)th_flavor[event_flavor][22+i]->Fill(SVMass[addJetIndex[i]][0],PU_weight);
			if(nSV[addJetIndex[i]]>1)th_flavor[event_flavor][24+i]->Fill(SVMass[addJetIndex[i]][1],PU_weight);
			
		}
		
		th_flavor[event_flavor][26]->Fill(dEta,PU_weight);
		th_flavor[event_flavor][27]->Fill(mjj,PU_weight);
		th_flavor[event_flavor][28]->Fill(mjjred,PU_weight);
		
		th_flavor[event_flavor][29]->Fill(nVtx);
		th_flavor[event_flavor][30]->Fill(nVtx,PU_weight);
		fixScaleNum[1]+=PU_weight;
	}
	TH1D * fixScale=new TH1D("fixScale","fixScale",2,-0.5,1.5);
	fixScale->SetBinContent(1,fixScaleNum[0]);
	fixScale->SetBinContent(2,fixScaleNum[1]);
	TFile* outFile= new TFile(Form("sf/%s_%d.root",st2.data(),w),"recreate");
	for(int i=0;i<31;i++){
		for(int j=0;j<4;j++)th_flavor[j][i]->Write();
	}
	fixScale->Write();
	outFile->Close();
}



void dataMC(int a,int b){

	string st1[40]={
		"/data7/chchen/Dec292016/JetHT/dataB/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataC/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataD/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataE/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataF/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataG/skimmed_NCUGlobalTuples_",
		"/data7/chchen/Dec292016/JetHT/dataH/skimmed_NCUGlobalTuples_",
		
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/C/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/D/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/E/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/F/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/G/skimmed/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/H/skimmed/NCUGlobalTuples_",
	};
	string  fileName[40]={
	"dataB",
	"dataC",
	"dataD",
	"dataE",
	"dataF",
	"dataG",
	"dataH",
	"dataH",
	};
	int aa[40]={1045,348,584,497,361,854,24};
	cout<<st1[a]<<b<<".root"<<endl;
	skimTree(b,st1[a],fileName[a]);	
	
}