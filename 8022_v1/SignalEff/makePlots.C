#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFrame.h>
#include <TLatex.h>
//#include "../macro/mkPlotsLivia/CMS_lumi.C"
#include <iostream>
#include <vector>
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
#include <fstream>
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
#include "../../untuplizer.h"


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

void makePlotBase(string st,double & db,double &db2){
	double nPass[30]={0},total=0,fixScaleNum[8]={0};

	TreeReader data(st.data());
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
		if(jEntry%1000==0)cout<<jEntry<<" out of "<<data.GetEntriesFast() <<" events are processed"<<endl;
		data.GetEntry(jEntry);
		bool isData  = data.GetBool("isData"); // only apply trigger if it's data

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
			//float tau21_puppi = FATjetPuppiTau2[ij]/FATjetPuppiTau1[ij];
			//if( tau21_puppi > 0.6 )continue;
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
			float raw_mass = fatjetPuppiSDmass[ij];
			thea_mass[ij] = raw_mass*thea_corr;

			if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;
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

		if(addjet_doublesv[addJetIndex[0]]< DBLoose)continue;
		if(addjet_doublesv[addJetIndex[1]]< DBLoose)continue;

		nPass[11]++;
		float mjj = (*higgsJet[0]+*higgsJet[1]).M();

		float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
		if(mjjred<750)continue;
		nPass[12]++;
		
		if(FATjetPuppiTau2[0]/FATjetPuppiTau1[0]>0.6 || FATjetPuppiTau2[1]/FATjetPuppiTau1[1]>0.6)continue;
		nPass[13]++;
		}
		db= double(nPass[12])/data.GetEntriesFast();
		db2= double(nPass[13])/data.GetEntriesFast();
	}
	

void makePlots(){
	double eff[10]={0},eff2[10]={0};
	int masspoint[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double masspointdb[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	
	for(int i=0;i<10;i++){
		makePlotBase(Form("/data7/chchen/Jan142017/B/B%d.root",masspoint[i]),eff[i],eff2[i]);
		cout<<eff[i]<<endl;
	}
	ofstream myfile;
	myfile.open ("effTable.txt");
	myfile<<"without Tau21 ";
	for(int i=0;i<10;i++)myfile<<" & "<<eff[i];
	myfile<<"////"<<endl;
	myfile<<"with Tau21 ";
	for(int i=0;i<10;i++)myfile<<" & "<<eff2[i];
	myfile<<"////"<<endl;
	TGraph* tg1 =new TGraph(10,masspointdb,eff);
	TGraph* tg2 =new TGraph(10,masspointdb,eff2);
	
	TFile* out=new TFile("eff.root","recreate");
	tg1->SetName("wo");
	tg2->SetName("with");
	tg1->Write();
	tg2->Write();
	out->Close();
	/*
	TCanvas * cboth = new TCanvas("cboth","",550,550);
	cboth->cd();
	gStyle->SetOptStat(0);
	 //gStyle->SetPaintTextFormat("2.1f");
	gStyle->SetPalette(57);
	gStyle->SetFrameLineWidth(3);
	//gStyle->SetPadRightMargin(0.109);
	 //gStyle->SetPadLeftMargin(0.13);
	 //cboth->SetLeftMargin(0.1);
	 //cboth->SetRightMargin(0.1);
	gStyle->SetTitleOffset(2, "Z");
	gStyle->SetNdivisions(605, "XYZ");
	
	tg1->GetXaxis()->SetTitle("m_{X} [GeV]");
	tg1->GetYaxis()->SetTitle("Efficiency");
	//limits->GetZaxis()->SetTitle("#sigma_{95\% CL} / #sigma_{th}");
	tg1->SetTitle("");
	tg1->SetMaximum(0.15);
	tg1->SetMarkerSize(1);
	tg1->SetMarkerStyle(20);
	tg1->SetMinimum(0.05);
	tg1->GetYaxis()->SetTitleOffset(1);
	//tg1->GetZaxis()->SetTitleOffset(0.65);
	// size of axis labels
	tg1->GetXaxis()->SetTitleSize(0.04);
	tg1->GetYaxis()->SetTitleSize(0.04);
	//limits->GetZaxis()->SetTitleSize(0.035);
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.03); 
	//limits->GetZaxis()->SetLabelSize(0.025);
	
	tg1->Draw("APL");
	cboth->Print("eff.pdf");
	*/
}