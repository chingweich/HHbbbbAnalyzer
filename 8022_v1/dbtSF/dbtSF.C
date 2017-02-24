#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
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
#include <TObject.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "../../HHbbbbAnalyzer/untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
//#include "../../setNCUStyle.C"


TH1D* makeEff(string input,double workingPoint=0.3){
	double bins[5]={250,350,430,840,2000};
	TH1D* th1=new TH1D("all","",4,bins);
	TH1D* th2=new TH1D("pass","",4,bins);
	TFile* tf1;
	tf1=TFile::Open(Form("%s",input.data()));
	TTree* tree;
	tf1->GetObject("mynewTree",tree);
	TreeReader data(tree);
				//data.Print();
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			Float_t  etadiff = data.GetFloat("etadiff");
			Float_t  jet1pt = data.GetFloat("jet1pt");
			Float_t  jet2pt = data.GetFloat("jet2pt");			
			//Float_t  dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
			Float_t  jet1bbtag = data.GetFloat("jet1bbtag");
			Float_t  jet2bbtag = data.GetFloat("jet2bbtag");
			Float_t  jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
			Float_t  jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
			if(jet1_puppi_msoftdrop_TheaCorr<105 ||jet1_puppi_msoftdrop_TheaCorr>135)continue;
			if(jet2_puppi_msoftdrop_TheaCorr<105 ||jet2_puppi_msoftdrop_TheaCorr>135)continue;
			//if(dijetmass_softdrop_corr<mjj*0.8||dijetmass_softdrop_corr>mjj*1.2)continue;
			//if(jet1bbtag<0.3 ||jet2bbtag<0.3)continue;
			//if(jet1bbtag>0.8 &&jet2bbtag>0.8)continue;
			if(etadiff>1.3)continue;
			th1->Fill(jet1pt);
			if(jet1bbtag>workingPoint)th2->Fill(jet1pt);
			th1->Fill(jet2pt);
			if(jet2bbtag>workingPoint)th2->Fill(jet2pt);
	}
	th2->Divide(th1);
	//tf1->Close();
	return th2;
}

void skimB(string input,string output,TH1D* th1[]){
  TFile* f;
  f =  TFile::Open(Form("%s",input.data()));
 // TDirectory * dir= (TDirectory*)f->Get(Form("%sNCUGlobalTuples_%d.root:/tree",st.data(),w));
  TTree* tree;
  f->GetObject("mynewTree",tree);
  //TTree* newtree = tree->CopyTree("AK8PuppijetSDmass[0]>50","",1000000000,0);
  TFile* outputfile = new TFile(Form("%s",output.data()),"RECREATE");
 TTree* mynewTree=new TTree("mynewTree","mynewTree");
		
		//Float_t ; mynewTree->Branch("",&,"/F");
		
		Float_t jet1pt; mynewTree->Branch("jet1pt",&jet1pt,"jet1pt/F");
		Float_t jet2pt; mynewTree->Branch("jet2pt",&jet2pt,"jet2pt/F");
		Float_t jet1eta; mynewTree->Branch("jet1eta",&jet1eta,"jet1eta/F");
		Float_t jet2eta; mynewTree->Branch("jet2eta",&jet2eta,"jet2eta/F");
		Float_t jet1phi; mynewTree->Branch("jet1phi",&jet1phi,"jet1phi/F");
		Float_t jet2phi; mynewTree->Branch("jet2phi",&jet2phi,"jet2phi/F");
		Float_t jet1mass; mynewTree->Branch("jet1mass",&jet1mass,"jet1mass/F");
		Float_t jet2mass; mynewTree->Branch("jet2mass",&jet2mass,"jet2mass/F");
		Float_t etadiff; mynewTree->Branch("etadiff",&etadiff,"etadiff/F");
		Float_t dijetmass; mynewTree->Branch("dijetmass",&dijetmass,"dijetmass/F");
		Float_t dijetmass_pruned_corr; mynewTree->Branch("dijetmass_pruned_corr",&dijetmass_pruned_corr,"dijetmass_pruned_corr/F");
		Float_t dijetmass_softdrop_corr; mynewTree->Branch("dijetmass_softdrop_corr",&dijetmass_softdrop_corr,"dijetmass_softdrop_corr/F");
		Float_t dijetmass_corr_punc; mynewTree->Branch("dijetmass_corr_punc",&dijetmass_corr_punc,"dijetmass_corr_punc/F");
		
		Float_t jet1pmass; mynewTree->Branch("jet1pmass",&jet1pmass,"jet1pmass/F");
		Float_t jet2pmass; mynewTree->Branch("jet2pmass",&jet2pmass,"jet2pmass/F");
		Float_t jet1bbtag; mynewTree->Branch("jet1bbtag",&jet1bbtag,"jet1bbtag/F");
		Float_t jet2bbtag; mynewTree->Branch("jet2bbtag",&jet2bbtag,"jet2bbtag/F");
		Float_t jet1_puppi_pt; mynewTree->Branch("jet1_puppi_pt",&jet1_puppi_pt,"jet1_puppi_pt/F");
		Float_t jet2_puppi_pt; mynewTree->Branch("jet2_puppi_pt",&jet2_puppi_pt,"jet2_puppi_pt/F");
		Float_t jet1_puppi_eta; mynewTree->Branch("jet1_puppi_eta",&jet1_puppi_eta,"jet1_puppi_eta/F");
		Float_t jet2_puppi_eta; mynewTree->Branch("jet2_puppi_eta",&jet2_puppi_eta,"jet2_puppi_eta/F");
		Float_t jet1_puppi_phi; mynewTree->Branch("jet1_puppi_phi",&jet1_puppi_phi,"jet1_puppi_phi/F");
		Float_t jet2_puppi_phi; mynewTree->Branch("jet2_puppi_phi",&jet2_puppi_phi,"jet2_puppi_phi/F");
		Float_t jet1_puppi_mass; mynewTree->Branch("jet1_puppi_mass",&jet1_puppi_mass,"jet1_puppi_mass/F");
		Float_t jet2_puppi_mass; mynewTree->Branch("jet2_puppi_mass",&jet2_puppi_mass,"jet2_puppi_mass/F");
		
		Float_t jet1_puppi_tau21; mynewTree->Branch("jet1_puppi_tau21",&jet1_puppi_tau21,"jet1_puppi_tau21/F");
		Float_t jet2_puppi_tau21; mynewTree->Branch("jet2_puppi_tau21",&jet2_puppi_tau21,"jet2_puppi_tau21/F");
		Float_t jet1_puppi_msoftdrop; mynewTree->Branch("jet1_puppi_msoftdrop",&jet1_puppi_msoftdrop,"jet1_puppi_msoftdrop/F");
		Float_t jet2_puppi_msoftdrop; mynewTree->Branch("jet2_puppi_msoftdrop",&jet2_puppi_msoftdrop,"jet2_puppi_msoftdrop/F");
		Float_t jet1_puppi_msoftdrop_TheaCorr; mynewTree->Branch("jet1_puppi_msoftdrop_TheaCorr",&jet1_puppi_msoftdrop_TheaCorr,"jet1_puppi_msoftdrop_TheaCorr/F");
		Float_t jet2_puppi_msoftdrop_TheaCorr; mynewTree->Branch("jet2_puppi_msoftdrop_TheaCorr",&jet2_puppi_msoftdrop_TheaCorr,"jet2_puppi_msoftdrop_TheaCorr/F");
		Float_t jet1_puppi_TheaCorr; mynewTree->Branch("jet1_puppi_TheaCorr",&jet1_puppi_TheaCorr,"jet1_puppi_TheaCorr/F");
		Float_t jet2_puppi_TheaCorr; mynewTree->Branch("jet2_puppi_TheaCorr",&jet2_puppi_TheaCorr,"jet2_puppi_TheaCorr/F");
		Float_t nTrueInt; mynewTree->Branch("nTrueInt",&nTrueInt,"nTrueInt/F");
		Float_t puWeights=1; mynewTree->Branch("puWeights",&puWeights,"puWeights/F");
		Float_t puWeightsUp=1; mynewTree->Branch("puWeightsUp",&puWeightsUp,"puWeightsUp/F");
		Float_t puWeightsDown=1; mynewTree->Branch("puWeightsDown",&puWeightsDown,"puWeightsDown/F");
		
		//Float_t isData; mynewTree->Branch("isData",&isData,"isData/F");
		
		Float_t SFTight=1; mynewTree->Branch("SFTight",&SFTight,"SFTight/F");
		Float_t SFTightup=1; mynewTree->Branch("SFTightup",&SFTightup,"SFTightup/F");
		Float_t SFTightdown=1; mynewTree->Branch("SFTightdown",&SFTightdown,"SFTightdown/F");
		Float_t SFLoose=1; mynewTree->Branch("SFLoose",&SFLoose,"SFLoose/F");
		Float_t SFLooseup=1; mynewTree->Branch("SFLooseup",&SFLooseup,"SFLooseup/F");
		Float_t SFLoosedown=1; mynewTree->Branch("SFLoosedown",&SFLoosedown,"SFLoosedown/F");
		
		Float_t trigWeightUp=1; mynewTree->Branch("trigWeightUp",&trigWeightUp,"trigWeightUp/F");
		Float_t trigWeightDown=1; mynewTree->Branch("trigWeightDown",&trigWeightDown,"trigWeightDown/F");
	 TreeReader data(tree);
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
	  data.GetEntry(jEntry);
	jet1pt = data.GetFloat("jet1pt");
	jet2pt = data.GetFloat("jet2pt");
	jet1eta = data.GetFloat("jet1eta");
	jet2eta = data.GetFloat("jet2eta");
	jet1phi = data.GetFloat("jet1phi");
	jet2phi = data.GetFloat("jet2phi");
	jet1mass = data.GetFloat("jet1mass");
	jet2mass = data.GetFloat("jet2mass");
	etadiff = data.GetFloat("etadiff");
	dijetmass = data.GetFloat("dijetmass");
	dijetmass_pruned_corr = data.GetFloat("dijetmass_pruned_corr");
	dijetmass_softdrop_corr = data.GetFloat("dijetmass_softdrop_corr");
	dijetmass_corr_punc = data.GetFloat("dijetmass_corr_punc");
	jet1pmass = data.GetFloat("jet1pmass");
	jet2pmass = data.GetFloat("jet2pmass");
	jet1bbtag = data.GetFloat("jet1bbtag");
	jet2bbtag = data.GetFloat("jet2bbtag");
	jet1_puppi_pt = data.GetFloat("jet1_puppi_pt");
	jet2_puppi_pt = data.GetFloat("jet2_puppi_pt");
	jet1_puppi_eta = data.GetFloat("jet1_puppi_eta");
	jet2_puppi_eta = data.GetFloat("jet2_puppi_eta");
	jet1_puppi_phi = data.GetFloat("jet1_puppi_phi");
	jet2_puppi_phi = data.GetFloat("jet2_puppi_phi");
	jet1_puppi_mass = data.GetFloat("jet1_puppi_mass");
	jet2_puppi_mass = data.GetFloat("jet2_puppi_mass");
	jet1_puppi_tau21 = data.GetFloat("jet1_puppi_tau21");
	jet2_puppi_tau21 = data.GetFloat("jet2_puppi_tau21");
	jet1_puppi_msoftdrop = data.GetFloat("jet1_puppi_msoftdrop");
	jet2_puppi_msoftdrop = data.GetFloat("jet2_puppi_msoftdrop");
	jet1_puppi_msoftdrop_TheaCorr = data.GetFloat("jet1_puppi_msoftdrop_TheaCorr");
	jet2_puppi_msoftdrop_TheaCorr = data.GetFloat("jet2_puppi_msoftdrop_TheaCorr");
	jet1_puppi_TheaCorr = data.GetFloat("jet1_puppi_TheaCorr");
	jet2_puppi_TheaCorr = data.GetFloat("jet2_puppi_TheaCorr");
	nTrueInt = data.GetFloat("nTrueInt");
	puWeights = data.GetFloat("puWeights");
	puWeightsUp = data.GetFloat("puWeightsUp");
	puWeightsDown = data.GetFloat("puWeightsDown");
	SFTight =1;
	SFTightup = 1;
	SFTightdown = 1;
	SFLoose =1;
	SFLooseup = 1;
	SFLoosedown = 1;
	
	for(int i=0;i<2;i++){
		double pt=0,bbtag=0;
		if(i==0){
			pt=jet1pt;
			bbtag=jet1bbtag;
		}
		else {
			pt=jet2pt;
			bbtag=jet2bbtag;
		}
		double SF=1,SFUp=1,SFDown=1,eff=1;
		if (pt<350){
			SF=0.96;
			SFUp=0.99;
			SFDown=0.94;
			eff=th1[0]->GetBinContent(1);
		}
		else if (pt <430){
			SF=1;
			SFUp=1.04;
			SFDown=0.97;
			eff=th1[0]->GetBinContent(2);
		}
		else if (pt <840){
			SF=1.01;
			SFUp=1.03;
			SFDown=0.97;
			eff=th1[0]->GetBinContent(3);
		}
		else {
			SF=1.01;
			SFUp=1.05;
			SFDown=0.93;
			eff=th1[0]->GetBinContent(4);
		}
		if(bbtag>0.3){
				SFLoose *=SF;
				SFLooseup *= SFUp;
				SFLoosedown *= SFDown;
		}
		else {
			SFLoose *=(1-SF*eff)/(1-eff);
			SFLooseup *=(1-SFUp*eff)/(1-eff);
			SFLoosedown *=(1-SFDown*eff)/(1-eff);
		}
	}
	
	for(int i=0;i<2;i++){
		double pt=0,bbtag=0;
		if(i==0){
			pt=jet1pt;
			bbtag=jet1bbtag;
		}
		else {
			pt=jet2pt;
			bbtag=jet2bbtag;
		}
		double SF=1,SFUp=1,SFDown=1,eff=1;
		if (pt<350){
			SF=0.92;
			SFUp=0.95;
			SFDown=0.89;
			eff=th1[1]->GetBinContent(1);
		}
		else if (pt <430){
			SF=1.01;
			SFUp=1.04;
			SFDown=0.97;
			eff=th1[1]->GetBinContent(2);
		}
		else if (pt <840){
			SF=0.92;
			SFUp=0.95;
			SFDown=0.87;
			eff=th1[1]->GetBinContent(3);
		}
		else {
			SF=0.92;
			SFUp=0.98;
			SFDown=0.82;
			eff=th1[1]->GetBinContent(4);
		}
		if(bbtag>0.8){
				SFTight *=SF;
				SFTightup *= SFUp;
				SFTightdown *= SFDown;
		}
		else {
			SFTight *=(1-SF*eff)/(1-eff);
			SFTightup *=(1-SFUp*eff)/(1-eff);
			SFTightdown *=(1-SFDown*eff)/(1-eff);
		}
	}
	
	
	  
	  mynewTree  ->Fill();
  }
 // outputfile->cd();
  mynewTree ->Write();
  //outputfile->Write();
    outputfile->Close();
}

void dbtSF(){
	TH1D* th1[2];
	int massP[13]={750,800,900,1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	for(int i=0;i<13;i++){
		th1[0]=makeEff(Form("../Summer2016GarvitonJECv3/BulkGrav_M-%d_0.root",massP[i]),0.3);
		th1[1]=makeEff(Form("../Summer2016GarvitonJECv3/BulkGrav_M-%d_0.root",massP[i]),0.8);
		skimB(Form("../Summer2016GarvitonJECv3/BulkGrav_M-%d_0.root",massP[i]),Form("BulkGrav_M-%d_0.root",massP[i]),th1);
	}
	int massP2[11]={750,800,900,1000,1200,1400,1600,2500,3000,3500,4500};
	for(int i=0;i<11;i++){
		th1[0]=makeEff(Form("../Summer2016RadionJECv3/Radion%d.root",massP2[i]),0.3);
		th1[1]=makeEff(Form("../Summer2016RadionJECv3/Radion%d.root",massP2[i]),0.8);
		skimB(Form("../Summer2016RadionJECv3/Radion%d.root",massP2[i]),Form("Radion%d.root",massP[i]),th1);
	}
	
}