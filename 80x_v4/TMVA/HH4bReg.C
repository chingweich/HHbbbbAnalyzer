//Histogram 
#include <TH1.h>
//#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TF1.h"
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

#define  nWidth 5
#define  nBmin 11

using namespace std;
using namespace TMVA;

float getPUPPIweight(float puppipt, float puppieta,TF1 *tf1[]){
	float genCorr  = 1.;
	float recoCorr = 1.;
	float totalWeight = 1.;
	
  
  
	genCorr =  tf1[0]->Eval( puppipt );
	if( fabs(puppieta)  <= 1.3 ) recoCorr = tf1[1]->Eval( puppipt );
	else if( fabs(puppieta) > 1.3 ) recoCorr = tf1[2]->Eval( puppipt );
  
	totalWeight = genCorr * recoCorr;

	return totalWeight;
	
}

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
	
	TFile *f3;
	f3=TFile::Open("puppiCorr.root");
	TF1* tf1[3];
	tf1[0]=(TF1 *) f3->FindObjectAny("puppiJECcorr_gen");
	tf1[1]=(TF1 *) f3->FindObjectAny("puppiJECcorr_reco_0eta1v3");
	tf1[2]=(TF1 *) f3->FindObjectAny("puppiJECcorr_reco_1v3eta2v5");
	
	
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
	
	TH1D* th2d[6];
	th2d[0]=new TH1D("unpr","unpr",4000,1000,5000);	
	th2d[1]=new TH1D("pr","pr",4000,1000,5000);	
	th2d[2]=new TH1D("sd","sd",4000,1000,5000);	
	
	th2d[3]=new TH1D("unprTM","unprTM",4000,1000,5000);	
	th2d[4]=new TH1D("prTM","prTM",4000,1000,5000);	
	th2d[5]=new TH1D("sdTM","sdTM",4000,1000,5000);	
	
		
	
	//int nWidth=5,nBmin=11;
	 int width [nWidth]={20,25,30,35,40};
	 int bmin[nBmin]={95,100,105,110,115,120,125,130,135,140,145};
	 
	 TH1D* th3d[6][nWidth][nBmin];
	 
	 for(int i=0;i<nWidth;i++){
		 for(int j=0;j<nBmin;j++){
			 th3d[0][i][j]=(TH1D*) th2d[0]->Clone(Form("%s_%d_%d",th2d[0]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[1][i][j]=(TH1D*) th2d[1]->Clone(Form("%s_%d_%d",th2d[1]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[2][i][j]=(TH1D*) th2d[2]->Clone(Form("%s_%d_%d",th2d[2]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 
			 th3d[3][i][j]=(TH1D*) th2d[3]->Clone(Form("%s_%d_%d",th2d[3]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[4][i][j]=(TH1D*) th2d[4]->Clone(Form("%s_%d_%d",th2d[4]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[5][i][j]=(TH1D*) th2d[5]->Clone(Form("%s_%d_%d",th2d[5]->GetTitle(),bmin[j],width[i]+bmin[j]));
		 }
	 }
  
   for(int i=0;i<3;i++)th2d[i]->Sumw2();
   
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
			//if(fatjetPRmassL2L3Corr[0]<90||fatjetPRmassL2L3Corr[0]>150)continue;
			//if(fatjetPRmassL2L3Corr[1]<90||fatjetPRmassL2L3Corr[1]>150)continue;
			nPass[8]++;
			//9.-----------------------------------------
			Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
			Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
			double tau21_1=(fatjetTau2[0]/fatjetTau1[0]),tau21_2=(fatjetTau2[1]/fatjetTau1[1]);
			
			
    
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
	   
	   varTemp[1]=val;
	   //cout<<val<<endl;
      }
	
	
	
	
	
	TClonesArray* FATjetPuppiP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");
	thisJet=(TLorentzVector*)FATjetPuppiP4->At(0);
	thatJet=(TLorentzVector*)FATjetPuppiP4->At(1);
	Float_t*  FATjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
	Float_t  FATjetPuppiSDmassThea[2] ={0};
			FATjetPuppiSDmassThea[0]=FATjetPuppiSDmass[0]*getPUPPIweight(thisJet->Pt(),thisJet->Eta(),tf1);
			FATjetPuppiSDmassThea[1]=FATjetPuppiSDmass[1]*getPUPPIweight(thatJet->Pt(),thatJet->Eta(),tf1);
			
			//cout<<thisJet->Pt()<<","<<thisJet->Eta()<<endl;
			
			//th6[4]->Fill(mjjRed);
			//th6[4]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M());
			
			
			
			 fixScaleNum[1]+=PU_weight[0];
			 
			
			//th2d[1]->Fill(fatjetPRmassL2L3Corr[0]*varTemp[0],fatjetPRmassL2L3Corr[1]*varTemp[1]);
	            //th2d[0]->Fill(fatjetPRmassL2L3Corr[0],fatjetPRmassL2L3Corr[1]);
			//th2d[2]->Fill(FATjetPuppiSDmassThea[0],FATjetPuppiSDmassThea[1]);
			
			//if( !((nDSVTight==2)||(nDSVTight==1 && nDSVMed==1)))continue;
			
			 for(int i=0;i<nWidth;i++){
				for(int j=0;j<nBmin;j++){
					if(width[i]+bmin[j]>166)continue;
					
					if(fatjetPRmassL2L3Corr[0]>bmin[j]&& fatjetPRmassL2L3Corr[0]<width[i]+bmin[j]
					&&fatjetPRmassL2L3Corr[1]>bmin[j] && fatjetPRmassL2L3Corr[1]<width[i]+bmin[j] && nDSVTight==2) 
					th3d[0][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
					if(fatjetPRmassL2L3Corr[0]*varTemp[0]>bmin[j] && fatjetPRmassL2L3Corr[0]*varTemp[0]<width[i]+bmin[j]
					&&fatjetPRmassL2L3Corr[1]*varTemp[1]>bmin[j] && fatjetPRmassL2L3Corr[1]*varTemp[1]<width[i]+bmin[j]  && nDSVTight==2)
					th3d[1][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
					if(FATjetPuppiSDmassThea[0]>bmin[j]&& FATjetPuppiSDmassThea[0]<width[i]+bmin[j]
					&&FATjetPuppiSDmassThea[1]>bmin[j]&& FATjetPuppiSDmassThea[1]<width[i]+bmin[j] && nDSVTight==2)  
					th3d[2][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
					if(fatjetPRmassL2L3Corr[0]>bmin[j] && fatjetPRmassL2L3Corr[0]<width[i]+bmin[j]
					&&fatjetPRmassL2L3Corr[1]>bmin[j] && fatjetPRmassL2L3Corr[1]<width[i]+bmin[j] && (nDSVTight==1 && nDSVMed==1)) 
					th3d[3][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
					if(fatjetPRmassL2L3Corr[0]*varTemp[0]>bmin[j] && fatjetPRmassL2L3Corr[0]*varTemp[0]<width[i]+bmin[j]
					&&fatjetPRmassL2L3Corr[1]*varTemp[1]>bmin[j]&& fatjetPRmassL2L3Corr[1]*varTemp[1]<width[i]+bmin[j]  && (nDSVTight==1 && nDSVMed==1))
					th3d[4][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
					if(FATjetPuppiSDmassThea[0]>bmin[j]&& FATjetPuppiSDmassThea[0]<width[i]+bmin[j]
					&&FATjetPuppiSDmassThea[1]>bmin[j]&& FATjetPuppiSDmassThea[1]<width[i]+bmin[j] && (nDSVTight==1 && nDSVMed==1))  
					th3d[5][i][j]->Fill((*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M(),PU_weight[0]);
					
				}
			}
			
		}
	}	
   
   //cout<<th2d->GetCorrelationFactor()<<endl;
  

 

   TH1D * fixScale=new TH1D("fixScale","fixScale",2,-0.5,1.5);
	fixScale->SetBinContent(1,fixScaleNum[0]);
	fixScale->SetBinContent(2,fixScaleNum[1]);
	
	TFile* outFile ;
	if(JESOption==0)outFile= new TFile(Form("category/%s.root",st2.data()),"recreate");
	else if(JESOption==1)outFile= new TFile(Form("category/%s_JESUp.root",st2.data()),"recreate");
	else if(JESOption==2)outFile= new TFile(Form("category/%s_JESDown.root",st2.data()),"recreate");
	fixScale->Write();
	
	for(int i=0;i<nWidth;i++){
		for(int j=0;j<nBmin;j++){
			if(width[i]+bmin[j]>166)continue;
			th3d[0][i][j]->Write();
			th3d[1][i][j]->Write();
			th3d[2][i][j]->Write();
			th3d[3][i][j]->Write();
			th3d[4][i][j]->Write();
			th3d[5][i][j]->Write();
		}
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