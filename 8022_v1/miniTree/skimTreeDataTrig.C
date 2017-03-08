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
  puppisd_corrRECO_cen->SetParameters(1.09302,-0.000150068,3.44866e-07,-2.68100e-10,8.67440e-14,-1.00114e-17);
  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters( 1.27212,-0.000571640,8.37289e-07,-5.20433e-10,1.45375e-13,-1.50389e-17);
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  totalWeight = genCorr * recoCorr;
  return totalWeight;
}


void skimTreeDataTrig(int w , string st){
	int nPass[20]={0},total=0;
	
	//option-----------------------------------------------------------
	
	
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
	
		TreeReader data(st.data());
		total+=data.GetEntriesFast();
		//cout<<"here"<<endl;
		TFile* outFile = new TFile(Form("minitree/%d.root",w),"recreate");
		//TTree* hh4b=new TTree("hh4b","hh4b");
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
		
		bool HLT_PFHT800=1; mynewTree->Branch("HLT_PFHT800",&HLT_PFHT800,"HLT_PFHT800/O");
		bool HLT_PFHT900=1; mynewTree->Branch("HLT_PFHT900",&HLT_PFHT900,"HLT_PFHT900/O");
		bool HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=1; mynewTree->Branch("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",&HLT_PFHT650_WideJetMJJ900DEtaJJ1p5,"HLT_PFHT650_WideJetMJJ900DEtaJJ1p5/O");
		bool HLT_AK8PFJet360_TrimMass30=1; mynewTree->Branch("HLT_AK8PFJet360_TrimMass30",&HLT_AK8PFJet360_TrimMass30,"HLT_AK8PFJet360_TrimMass30/O");
		bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=1; mynewTree->Branch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",&HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20,"HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20/O");
		bool HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT650_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT650_TrimR0p1PT0p03Mass50/O");
		bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=1; mynewTree->Branch("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50",&HLT_AK8PFHT700_TrimR0p1PT0p03Mass50,"HLT_AK8PFHT700_TrimR0p1PT0p03Mass50/O");
		
		
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
			data.GetEntry(jEntry);
			
			bool isData  = data.GetBool("isData"); // only apply trigger if it's data

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    HLT_PFHT800=0;
    HLT_PFHT900=0;
    HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=0;
    HLT_AK8PFJet360_TrimMass30=0;
    HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=0;
    HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=0;
    HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=0;
    
    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];
		if(thisTrig.find("HLT_PFHT800_v")!= std::string::npos )HLT_PFHT800=results;
		if(thisTrig.find("HLT_PFHT900_v")!= std::string::npos )HLT_PFHT900=results;
		if(thisTrig.find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v")!= std::string::npos )HLT_PFHT650_WideJetMJJ900DEtaJJ1p5=results;
		if(thisTrig.find("HLT_AK8PFJet360_TrimMass30_v")!= std::string::npos )HLT_AK8PFJet360_TrimMass30=results;
		if(thisTrig.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v")!= std::string::npos )HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=results;
		if(thisTrig.find("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v")!= std::string::npos )HLT_AK8PFHT650_TrimR0p1PT0p03Mass50=results;
		if(thisTrig.find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v")!= std::string::npos )HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=results;
      }

    //if(!passTrigger && isData)continue;
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
	if( tau21_puppi > 0.55 )continue;

	TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
	float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
	float raw_mass = fatjetPuppiSDmass[ij];
	thea_mass[ij] = raw_mass*thea_corr;

	//if(thea_mass[ij] < 105 || thea_mass[ij] > 135)continue;


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

    //if(addjet_doublesv[addJetIndex[0]]< DBLoose)continue;
    //if(addjet_doublesv[addJetIndex[1]]< DBLoose)continue;


    nPass[10]++;
    TLorentzVector* higgsJet[2];
    for(int i=0;i<2;i++)higgsJet[i] = (TLorentzVector*)fatjetP4->At(i);
    float dEta = fabs(higgsJet[0]->Eta()-higgsJet[1]->Eta());

    //if(dEta>1.3)continue;
    nPass[11]++;
    float mjj = (*higgsJet[0]+*higgsJet[1]).M();

    float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
    //if(mjjred<750)continue;
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
			//Float_t*  FATjetPuppiSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr");
			Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
			Float_t*  FATjetSDmass = data.GetPtrFloat("FATjetSDmass");
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuoEF = data.GetPtrFloat("FATjetMuoEF");
			
			Float_t  ptHat = data.GetFloat("ptHat");

			Float_t ntrue= data.GetFloat("pu_nTrueInt");
			double PU_weight[3]={1,1,1};
			
			
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
			puWeights=PU_weight[0] ;
			puWeightsUp=PU_weight[1];
			puWeightsDown=PU_weight[2];
			
			nTrueInt=ntrue;
			//isData=isData;
			jet1pt=higgsJet[0]->Pt();
			jet1eta=higgsJet[0]->Eta();
			jet1phi=higgsJet[0]->Phi();
			jet1mass=higgsJet[0]->M();
			jet2pt=higgsJet[1]->Pt();
			jet2eta=higgsJet[1]->Eta();
			jet2phi=higgsJet[1]->Phi();
			jet2mass=higgsJet[1]->M();
			
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(0);
			TLorentzVector* thatJetPuppi = (TLorentzVector*)puppijetP4->At(1);
			jet1_puppi_pt=thisJetPuppi->Pt();
			jet1_puppi_eta=thisJetPuppi->Eta();
			jet1_puppi_phi=thisJetPuppi->Phi();
			jet1_puppi_mass=thisJetPuppi->M();
			jet2_puppi_pt=thatJetPuppi->Pt();
			jet2_puppi_eta=thatJetPuppi->Eta();
			jet2_puppi_phi=thatJetPuppi->Phi();
			jet2_puppi_mass=thatJetPuppi->M();
			
			etadiff=dEta;
			dijetmass=mjj;
			dijetmass_softdrop_corr=mjjred;
			
			dijetmass_pruned_corr=mjj;
			dijetmass_corr_punc=mjj;
			
			jet1pmass=fatjetPRmassL2L3Corr[0];
			jet2pmass=fatjetPRmassL2L3Corr[1];
			jet1bbtag=addjet_doublesv[addJetIndex[0]];
			jet2bbtag=addjet_doublesv[addJetIndex[1]];
		
			jet1_puppi_tau21=FATjetPuppiTau2[0]/FATjetPuppiTau1[0];
			jet2_puppi_tau21=FATjetPuppiTau2[1]/FATjetPuppiTau1[1];
			
			jet1_puppi_msoftdrop=fatjetPuppiSDmass[0];
			jet2_puppi_msoftdrop=fatjetPuppiSDmass[1];
		
			jet1_puppi_msoftdrop_TheaCorr=thea_mass[0];
			jet2_puppi_msoftdrop_TheaCorr=thea_mass[1];
			
			jet1_puppi_TheaCorr=thea_mass[0]/fatjetPuppiSDmass[0];
			jet2_puppi_TheaCorr=thea_mass[1]/fatjetPuppiSDmass[1];
			
			mynewTree ->Fill();
			
		}//end event loop----------------------------------------------------------------------------------------
		mynewTree->Write();
		
		TH1F* th1f= new  TH1F("CountWeighted","",1,0,1);
		th1f->SetBinContent(1,total);
		th1f->Write();
		
		outFile->Close();
		for(int i=0;i<13;i++)cout<<"npass["<<i<<"]="<<nPass[i]<<endl;
	
	
	
	
}
