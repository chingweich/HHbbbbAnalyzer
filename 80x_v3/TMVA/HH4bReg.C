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
#include <cstdlib>
#include <map>
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

//math 
#include <cmath>
//#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "../untuplizer.h"
//#include "jetEnergyScale.h"

#include "../standalone_LumiReWeighting.cc"
#include "../standalone_LumiReWeighting.h"

#include "../BTagCalibrationStandalone.h"
#include "../BTagCalibrationStandalone.cpp"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "HH4bCategoryBase_80.C"
//#include "TMVARegressionApplication.C"

#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace TMVA;

void TMVARegressionApplication( int wMs,int wM, string st,string st2,string option="",TString myMethodList = "" ) 
{
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0; 
   Use["KNN"]             = 0;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		        = 0;
   // 
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"]             = 0; 
   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   //Float_t var1, var2;
   //reader->AddVariable( "var1", &var1 );
   //reader->AddVariable( "var2", &var2 );
   Float_t pt_AK8MatchedToHbb,eta_AK8MatchedToHbb,nsv_AK8MatchedToHbb,sv0mass_AK8MatchedToHbb,sv1mass_AK8MatchedToHbb,
   nch_AK8MatchedToHbb,nmu_AK8MatchedToHbb,nel_AK8MatchedToHbb,muenfr_AK8MatchedToHbb,emenfr_AK8MatchedToHbb;
   reader->AddVariable( "pt_AK8MatchedToHbb", &pt_AK8MatchedToHbb );
   reader->AddVariable( "eta_AK8MatchedToHbb", &eta_AK8MatchedToHbb );
   reader->AddVariable( "nsv_AK8MatchedToHbb", &nsv_AK8MatchedToHbb );
   reader->AddVariable( "sv0mass_AK8MatchedToHbb", &sv0mass_AK8MatchedToHbb );
   reader->AddVariable( "sv1mass_AK8MatchedToHbb", &sv1mass_AK8MatchedToHbb );
   reader->AddVariable( "nch_AK8MatchedToHbb", &nch_AK8MatchedToHbb );
   reader->AddVariable( "nmu_AK8MatchedToHbb", &nmu_AK8MatchedToHbb );
   reader->AddVariable( "nel_AK8MatchedToHbb", &nel_AK8MatchedToHbb );
   reader->AddVariable( "muenfr_AK8MatchedToHbb", &muenfr_AK8MatchedToHbb );
   reader->AddVariable( "emenfr_AK8MatchedToHbb", &emenfr_AK8MatchedToHbb );

   
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t spec1,spec2;
   reader->AddSpectator( "spec1:=n_pv",  &spec1 );
   reader->AddSpectator( "spec2:=mpruned_AK8MatchedToHbb",  &spec2 );

   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVARegression";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
     TH1* hists[100];
   Int_t nhists = -1;
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
         TH1* h = new TH1F( it->first.c_str(), TString(it->first) + " method", 100, -100, 600 );
         if (it->second) hists[++nhists] = h;
      }
      nhists++;
   
  //1=signal ,0=QCD ,2=data
	int nameRoot=1;
	if((st2.find("QCD")!= std::string::npos)||
	(st2.find("bGen")!= std::string::npos)||
	(st2.find("bEnriched")!= std::string::npos))nameRoot=0;
	if(st2.find("data")!= std::string::npos)nameRoot=2;
	cout<<"nameRoot = "<<nameRoot<<endl;
	
	//option-----------------------------------------------------------
	
	int JESOption=0;
	if(option.find("JESUp")!= std::string::npos)JESOption=1;
	if(option.find("JESDown")!= std::string::npos)JESOption=2;
	cout<<"JESOption = "<<JESOption<<endl;
	
	
	//using for Btag Eff -----------------------------------------------------------------------------
	string btagsystematicsType="central";
	if(JESOption==3)btagsystematicsType="up";
	else if(JESOption==4)btagsystematicsType="down";
	BTagCalibration calib("CSVv2L", "../subjet_CSVv2_ichep.csv");
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
	if(nameRoot==2)f1=TFile::Open("../btagEffSource/data.root");
	else if (nameRoot!=2 && (JESOption==0||JESOption==3||JESOption==4||JESOption==5||JESOption==6))f1=TFile::Open(Form("../btagEffSource/%s.root",st2.data()));
	else if (nameRoot!=2 && JESOption==1)f1=TFile::Open(Form("btagEffSource/%s_JESUp.root",st2.data()));
	else if (nameRoot!=2 && JESOption==2)f1=TFile::Open(Form("btagEffSource/%s_JESDown.root",st2.data()));
	TH1D* th1[6];
	string btaggingEff[6]={"effD_b","effN_b","effD_c","effN_c","effD_l","effN_l"};
	for(int i=0;i<6;i++){
		th1[i]=(TH1D*)f1->FindObjectAny(Form("%s_1d",btaggingEff[i].data()));
		if(i==1||i==3||i==5)th1[i]->Divide(th1[i-1]);
	}
	
	standalone_LumiReWeighting LumiWeights_central(0),LumiWeights_up(1),LumiWeights_down(-1);
  
  
   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *f;
	TTree *tree;
	int nPass[20]={0};
	int total=0;
	double fixScaleNum[2]={0};
	double xBinsForHeavyFlavor[5]={30,140,180,240,3000};
	double xBinsForLightFlavor[11]={20,100,200,300,400,500,600,700,800,900,3000};
	
	TH1D* th5[30];
	th5[0]=new TH1D("cat0","cat0",4000,1000,5000);
	th5[1]=new TH1D("cat1","cat1",4000,1000,5000);
	th5[2]=new TH1D("cat2","cat2",4000,1000,5000);
	
	th5[3]=new TH1D("TT","TT",4000,1000,5000);
	th5[4]=new TH1D("TM","TM",4000,1000,5000);
	th5[5]=new TH1D("TL","TL",4000,1000,5000);
	th5[6]=new TH1D("MM","MM",4000,1000,5000);
	th5[7]=new TH1D("ML","ML",4000,1000,5000);
	th5[8]=new TH1D("LL","LL",4000,1000,5000);
	
	th5[9]=new TH1D("DBTLL","DBTLL",4000,1000,5000);
	th5[10]=new TH1D("DBTML","DBTML",4000,1000,5000);
	th5[11]=new TH1D("DBTTL","DBTTL",4000,1000,5000);
	
	th5[12]=new TH1D("DBTLL1","DBTLL1",4000,1000,5000);
	th5[13]=new TH1D("DBTML1","DBTML1",4000,1000,5000);
	th5[14]=new TH1D("DBTTL1","DBTTL1",4000,1000,5000);
	
	th5[15]=new TH1D("TMe","TMe",4000,1000,5000);
	th5[16]=new TH1D("TLe","TLe",4000,1000,5000);
	th5[17]=new TH1D("MMe","MMe",4000,1000,5000);
	th5[18]=new TH1D("MLe","MLe",4000,1000,5000);
	th5[19]=new TH1D("LLe","LLe",4000,1000,5000);

	for(int i=0;i<20;i++){
		th5[i]->Sumw2();
	}
	
		
		
		TH1D* th6[6];
		th6[0]= new TH1D("regMass_j0","",60,90,150);
		th6[1]= new TH1D("regMass_j1","",60,90,150);
		th6[2]= new TH1D("prMass_j0","",60,90,150);
		th6[3]= new TH1D("prMass_j1","",60,90,150);
		th6[4]= new TH1D("regMjj","",60,90,150);
		th6[5]= new TH1D("prMjj","",60,90,150);
   
for(int i=0;i<6;i++){
		th6[i]->Sumw2();
	}
   
   for (int w=wMs;w<wM;w++){
		if(w%20==0)cout<<w<<endl;
		
		if (nameRoot!=1)f = TFile::Open(Form("%s%d.root",st.data(),w));
		else f = TFile::Open(st.data());
		if (!f || !f->IsOpen())continue;
		
		TDirectory * dir;
		if (nameRoot!=1)dir = (TDirectory*)f->Get(Form("%s%d.root:/tree",st.data(),w));
		else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		
		dir->GetObject("treeMaker",tree);
		TreeReader data(tree);
		total+=data.GetEntriesFast();
		for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
			data.GetEntry(jEntry);
			
			
			double PU_weight[3]={1,1,1};
			Float_t ntrue= data.GetFloat("pu_nTrueInt");
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
			fixScaleNum[0]+=PU_weight[0];
			
			Int_t nVtx        = data.GetInt("nVtx");
			//0. has a good vertex
			if(nVtx<1)continue;
			nPass[0]++;
			
			//1.trigger
			std::string* trigName = data.GetPtrString("hlt_trigName");
		 	vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
			bool passTrigger=false;
			for(int it=0; it< data.GetPtrStringSize();it++){
				std::string thisTrig= trigName[it];
				bool results = trigResult[it];
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
			if(!passTrigger && nameRoot==2)continue;
			nPass[1]++;

			const int nFATJet=data.GetInt("FATnJet");
			//2.nJets
			if(nFATJet<2)continue;nPass[2]++;
			TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
			float*  FATjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp"); 
			float*  FATjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown"); 
			TLorentzVector* thisJet ,* thatJet;
			if(JESOption==0){
				thisJet=(TLorentzVector*)fatjetP4->At(0);
				thatJet=(TLorentzVector*)fatjetP4->At(1);
			}
			else if (JESOption==1){
				TLorentzVector thisJetScaleUp= (*((TLorentzVector*)fatjetP4->At(0)))*(1+FATjetCorrUncUp[0]);
				TLorentzVector thatJetScaleUp= (*((TLorentzVector*)fatjetP4->At(1)))*(1+FATjetCorrUncUp[1]);
				thisJet= &thisJetScaleUp;
				thatJet= &thatJetScaleUp;
			}
			else if (JESOption==2){
				TLorentzVector thisJetScaleDown= (*((TLorentzVector*)fatjetP4->At(0)))*(1-FATjetCorrUncDown[0]);
				TLorentzVector thatJetScaleDown= (*((TLorentzVector*)fatjetP4->At(1)))*(1-FATjetCorrUncDown[1]);
				thisJet= &thisJetScaleDown;
				thatJet= &thatJetScaleDown;
			}
			
			//3. Pt 
			if(thisJet->Pt()<200)continue;
			if(thatJet->Pt()<200)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
			if(FATjetPassIDTight[0]==0)continue;
			if(FATjetPassIDTight[1]==0)continue;
			Float_t*  FATjetCEmEF = data.GetPtrFloat("FATjetCEmEF");
			Float_t*  FATjetMuEF = data.GetPtrFloat("FATjetMuoEF");
			if(FATjetMuEF[0]>0.8)continue;
			if(FATjetCEmEF[0]>0.9)continue;
			if(FATjetMuEF[1]>0.8)continue;
			if(FATjetCEmEF[1]>0.9)continue;
			nPass[4]++;
			//5. Eta-----------------------------------------
			if(fabs(thisJet->Eta())>2.4)continue;
			if(fabs(thatJet->Eta())>2.4)continue;
			nPass[5]++;
			//6. DEta-----------------------------------------
			float dEta = fabs(thisJet->Eta()-thatJet->Eta());
			if(dEta>1.3)continue;
			nPass[6]++;
			//7. Mjj-----------------------------------------
			float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			if(mjjRed<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
			if(fatjetPRmassL2L3Corr[0]<90||fatjetPRmassL2L3Corr[0]>150)continue;
			if(fatjetPRmassL2L3Corr[1]<90||fatjetPRmassL2L3Corr[1]>150)continue;
			nPass[8]++;
			//9.-----------------------------------------
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			
			
		
	   Int_t* FATjet_nSV=data.GetPtrInt("FATjet_nSV");
	   vector<float>   *FATjet_SVMass  =  data.GetPtrVectorFloat("FATjet_SVMass");
	   int nEle= data.GetInt("nEle");
	   int nMu=data.GetInt("nMu");
	   Float_t*  FATjetEleEF = data.GetPtrFloat("FATjetEleEF");
	Int_t* FATjetCMulti=data.GetPtrInt("FATjetCMulti");
	Int_t* FATjetEleMulti=data.GetPtrInt("FATjetEleMulti");
	Int_t* FATjetMuoMulti=data.GetPtrInt("FATjetMuoMulti");
			 pt_AK8MatchedToHbb=thisJet->Pt();
		eta_AK8MatchedToHbb=thisJet->Eta();
		nsv_AK8MatchedToHbb=FATjet_nSV[0];
		sv0mass_AK8MatchedToHbb=FATjet_SVMass[0][0];
		sv1mass_AK8MatchedToHbb=FATjet_SVMass[0][1];
		nmu_AK8MatchedToHbb=FATjetMuoMulti[0];
		nel_AK8MatchedToHbb=FATjetEleMulti[0];
		muenfr_AK8MatchedToHbb=FATjetMuEF[0];
		nch_AK8MatchedToHbb=FATjetCMulti[0];
		emenfr_AK8MatchedToHbb=FATjetEleEF[0];
		
		Float_t*  FATjetPRmass = data.GetPtrFloat("FATjetPRmass");
		
		double varTemp[2];
		
		 for (Int_t ih=0; ih<nhists; ih++) {
         TString title = hists[ih]->GetTitle();
         Float_t val = (reader->EvaluateRegression( title ))[0];
         hists[ih]->Fill( val );   
	   th6[0]->Fill(FATjetPRmass[0]*val);
	   th6[2]->Fill(fatjetPRmassL2L3Corr[0]);
	   varTemp[0]=val;
	 //  cout<<val<<endl;
      }
	
	 pt_AK8MatchedToHbb=thatJet->Pt();
		eta_AK8MatchedToHbb=thatJet->Eta();
		nsv_AK8MatchedToHbb=FATjet_nSV[1];
		sv0mass_AK8MatchedToHbb=FATjet_SVMass[1][0];
		sv1mass_AK8MatchedToHbb=FATjet_SVMass[1][1];
		nmu_AK8MatchedToHbb=FATjetMuoMulti[1];
		nel_AK8MatchedToHbb=FATjetEleMulti[1];
		muenfr_AK8MatchedToHbb=FATjetMuEF[1];
		nch_AK8MatchedToHbb=FATjetCMulti[1];
		emenfr_AK8MatchedToHbb=FATjetEleEF[1];
		
		 for (Int_t ih=0; ih<nhists; ih++) {
         TString title = hists[ih]->GetTitle();
         Float_t val = (reader->EvaluateRegression( title ))[0];
         hists[ih]->Fill( val );   
	   th6[1]->Fill(FATjetPRmass[1]*val);
	   th6[3]->Fill(fatjetPRmassL2L3Corr[1]);
	   varTemp[1]=val;
	   //cout<<val<<endl;
      }
	
			if(fatjetPRmassL2L3Corr[0]<105||fatjetPRmassL2L3Corr[0]>135)continue;
			if(fatjetPRmassL2L3Corr[1]<105||fatjetPRmassL2L3Corr[1]>135)continue;
			th6[4]->Fill(mjjRed);
			th6[4]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M());
			
			
			//10.btag
			vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV");
			int nbtag=0;
			if(subjetSDCSV[0][0]>0.46)nbtag++;
			if(subjetSDCSV[0][1]>0.46)nbtag++;
			if(subjetSDCSV[1][0]>0.46)nbtag++;
			if(subjetSDCSV[1][1]>0.46)nbtag++;
			
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx");
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy");
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz");
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE");
			vector<Int_t>   *FATsubjetSDHadronFlavor =  data.GetPtrVectorInt("FATsubjetSDHadronFlavor");
			
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
			
			//if(tau21_1>0.75 || tau21_2>0.75) continue;
			
			//if(nbtag==3 && ((tau21_1>0.6 && tau21_2<0.6)||(tau21_1<0.6 && tau21_2>0.6)))th1[2]->Fill(mjjRed);
			
			if(tau21_1>0.6 || tau21_2>0.6) continue;
			if(nbtag==4)th5[0]->Fill(mjjRed,btaggingscaleFactor*PU_weight[0]);
			if(nbtag==3)th5[1]->Fill(mjjRed,btaggingscaleFactor*PU_weight[0]);
			if(nbtag==2)th5[2]->Fill(mjjRed,btaggingscaleFactor*PU_weight[0]);
			nPass[9]++;
			
			 fixScaleNum[1]+=PU_weight[0];
			
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
				if(matchThis&& matchThat)break;
			}
			
			int nDSVTight=0,nDSVMed=0,nDSVLoose=0;
			if(FATjet_DoubleSV[0]>0.8)nDSVTight++;
			else if (FATjet_DoubleSV[0]>0.6)nDSVMed++;
			else if (FATjet_DoubleSV[0]>0.3)nDSVLoose++;
			
			if(FATjet_DoubleSV[1]>0.8)nDSVTight++;
			else if (FATjet_DoubleSV[1]>0.6)nDSVMed++;
			else if (FATjet_DoubleSV[1]>0.3)nDSVLoose++;
			
			if(nDSVTight==2)th5[3]->Fill(mjjRed,PU_weight[0]);
			if(nDSVTight==1 && nDSVMed==1)th5[4]->Fill(mjjRed,PU_weight[0]);
			if(nDSVTight==1 && nDSVLoose==1)th5[5]->Fill(mjjRed,PU_weight[0]);
			if(nDSVMed==2)th5[6]->Fill(mjjRed,PU_weight[0]);
			if(nDSVMed==1 && nDSVLoose==1)th5[7]->Fill(mjjRed,PU_weight[0]);
			if(nDSVLoose==2)th5[8]->Fill(mjjRed,PU_weight[0]);
			
			if(nDSVTight==2)continue;
			
			if(nDSVTight==1 && nDSVMed==1)th5[15]->Fill(mjjRed,PU_weight[0]);
			if(nDSVTight==1 && nDSVLoose==1)th5[16]->Fill(mjjRed,PU_weight[0]);
			if(nDSVMed==2)th5[17]->Fill(mjjRed,PU_weight[0]);
			if(nDSVMed==1 && nDSVLoose==1)th5[18]->Fill(mjjRed,PU_weight[0]);
			if(nDSVLoose==2)th5[19]->Fill(mjjRed,PU_weight[0]);
			
			if((FATjet_DoubleSV[0]>0.3&& subjetSDCSV[1][0]>0.46 && subjetSDCSV[1][1]>0.46)||
			(FATjet_DoubleSV[1]>0.3&& subjetSDCSV[0][0]>0.46 && subjetSDCSV[0][1]>0.46))th5[9]->Fill(mjjRed,PU_weight[0]);
			if((FATjet_DoubleSV[0]>0.6&& subjetSDCSV[1][0]>0.46 && subjetSDCSV[1][1]>0.46)||
			(FATjet_DoubleSV[1]>0.6&& subjetSDCSV[0][0]>0.46 && subjetSDCSV[0][1]>0.46))th5[10]->Fill(mjjRed,PU_weight[0]);
			if((FATjet_DoubleSV[0]>0.8&& subjetSDCSV[1][0]>0.46 && subjetSDCSV[1][1]>0.46)||
			(FATjet_DoubleSV[1]>0.8&& subjetSDCSV[0][0]>0.46 && subjetSDCSV[0][1]>0.46))th5[11]->Fill(mjjRed,PU_weight[0]);
			
			if((FATjet_DoubleSV[0]>0.3&& (subjetSDCSV[1][0]>0.46 || subjetSDCSV[1][1]>0.46))||
			(FATjet_DoubleSV[1]>0.3&& (subjetSDCSV[0][0]>0.46 || subjetSDCSV[0][1]>0.46)))th5[12]->Fill(mjjRed,PU_weight[0]);
			if((FATjet_DoubleSV[0]>0.6&& (subjetSDCSV[1][0]>0.46 || subjetSDCSV[1][1]>0.46))||
			(FATjet_DoubleSV[1]>0.6&& (subjetSDCSV[0][0]>0.46 || subjetSDCSV[0][1]>0.46)))th5[13]->Fill(mjjRed,PU_weight[0]);
			if((FATjet_DoubleSV[0]>0.8&& (subjetSDCSV[1][0]>0.46 || subjetSDCSV[1][1]>0.46))||
			(FATjet_DoubleSV[1]>0.8&& (subjetSDCSV[0][0]>0.46 || subjetSDCSV[0][1]>0.46)))th5[14]->Fill(mjjRed,PU_weight[0]);
		}
	}	
   
   
  

 

   TH1D * fixScale=new TH1D("fixScale","fixScale",2,-0.5,1.5);
	fixScale->SetBinContent(1,fixScaleNum[0]);
	fixScale->SetBinContent(2,fixScaleNum[1]);
	
	TFile* outFile ;
	if(JESOption==0)outFile= new TFile(Form("category/%s.root",st2.data()),"recreate");
	else if(JESOption==1)outFile= new TFile(Form("category/%s_JESUp.root",st2.data()),"recreate");
	else if(JESOption==2)outFile= new TFile(Form("category/%s_JESDown.root",st2.data()),"recreate");
	fixScale->Write();
	for(int i=0;i<20;i++){
		th5[i]->Write();
	}
	for(int i=0;i<6;i++){
		th6[i]->Write();
	}
	outFile->Close();
  
   delete reader;
    
   
}

void HH4bReg(int a){

	string st1[40]={
		/*0-11*/
		"/data7/chchen/8011_regression/B/B1000.root",
		"/data7/chchen/8011_regression/B/B1200.root",
		"/data7/chchen/8011_regression/B/B1400.root",
		"/data7/chchen/8011_regression/B/B1600.root",
		"/data7/chchen/8011_regression/B/B1800.root",
		"/data7/chchen/8011_regression/B/B2000.root",
		"/data7/chchen/8011_regression/B/B2500.root",
		"/data7/chchen/8011_regression/B/B3000.root",
		"/data7/chchen/8011_regression/B/B3500.root",
		"/data7/chchen/8011_regression/B/B4000.root",
		"/data7/chchen/8011_regression/B/B4500.root",
		/*12-21*/
		"",//"/data7/chchen/8011_regression/R/R1000.root",
		"",//"/data7/chchen/8011_regression/R/R1200.root",
		"",//"/data7/chchen/8011_regression/R/R1400.root",
		"",//"/data7/chchen/8011_regression/R/R1600.root",
		"",//"/data7/chchen/8011_regression/R/R1800.root",
		"",//"/data7/chchen/8011_regression/R/R2000.root",
		"",//"/data7/chchen/8011_regression/R/R2500.root",
		"",//"/data7/chchen/8011_regression/R/R3000.root",
		"",//"/data7/chchen/8011_regression/R/R3500.root",
		"",//"/data7/chchen/8011_regression/R/R4000.root",
		"",//"/data7/chchen/8011_regression/R/R4500.root",
		/*22-32*/
		"/data7/chchen/8011_regression/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/8011_regression/QCD/2000_1/NCUGlobalTuples_",	
		
		"",//"/data7/chchen/8011_regression/QCDHTbGen/700/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbGen/1000/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbGen/1500/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbGen/2000/NCUGlobalTuples_",
		
		"",//"/data7/chchen/8011_regression/QCDHTbEnriched/700/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbEnriched/1000/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbEnriched/1500/NCUGlobalTuples_",
		"",//"/data7/chchen/8011_regression/QCDHTbEnriched/2000/NCUGlobalTuples_",
	};
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2",
	"bGen700","bGen1000","bGen1500","bGen2000",
	"bEnriched700","bEnriched1000","bEnriched1500","bEnriched2000",
	};
	int aa[40]={119,225,40,80,31,62,17,33,27,12,5,3,9,3,2,3};
	string dataPathB="/data7/chchen/8011_regression/JetHT_runB/NCUGlobalTuples_";
	string dataPathC="/data7/chchen/8011_regression/JetHT_runC/NCUGlobalTuples_";
	string dataPathD="/data7/chchen/8011_regression/JetHT_runD/NCUGlobalTuples_";
	
	if(a==38)TMVARegressionApplication(1,400,dataPathB,"data1");	
	else if(a==39)TMVARegressionApplication(400,800,dataPathB,"data2");	
	else if(a==40)TMVARegressionApplication(800,1048,dataPathB,"data3");	
	else if(a==41)TMVARegressionApplication(1,348,dataPathC,"data4");	
	else if(a==42)TMVARegressionApplication(1,300,dataPathD,"data5");	
	else if(a==43)TMVARegressionApplication(300,584,dataPathD,"data6");
	else if (a>21)TMVARegressionApplication(1,aa[a-22],st1[a],fileName[a],"");
	else TMVARegressionApplication(1,2,st1[a],fileName[a],"");
}