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
   reader->AddSpectator( "spec2:=msoftdrop_AK8MatchedToHbb",  &spec2 );

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
	
	
	TFile *f;
	TTree *tree;
	int nPass[20]={0};
	int total=0;
	
	TH1D* th1;
	th1=new TH1D("mass","mass",150,0,150);
	
	TH1D* th3;
	th3=new TH1D("mass","mass",1100,200,2400);
	
	double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	
	TH1D* th2[6][14];
	
	for(int i=0;i<14;i++){
		th2[0][i]=(TH1D*)th1->Clone(Form("genBarelMass%.0f",ptBins[i]));
		th2[1][i]=(TH1D*)th1->Clone(Form("genEndcapMass%.0f",ptBins[i]));
		th2[2][i]=(TH1D*)th1->Clone(Form("recoBarelMass%.0f",ptBins[i]));
		th2[3][i]=(TH1D*)th1->Clone(Form("recoEndcapMass%.0f",ptBins[i]));
		th2[4][i]=(TH1D*)th3->Clone(Form("ptBarel%.0f",ptBins[i]));
		th2[5][i]=(TH1D*)th3->Clone(Form("ptEndcap%.0f",ptBins[i]));
	}

	for(int i=0;i<14;i++){
		th2[0][i]->Sumw2();
		th2[1][i]->Sumw2();
		th2[2][i]->Sumw2();
		th2[3][i]->Sumw2();
		th2[4][i]->Sumw2();
		th2[5][i]->Sumw2();
	}
	
	//---------------------------------
	
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
			if(jEntry%2)continue;
			data.GetEntry(jEntry);
			
			
			
			Int_t nGenPar        = data.GetInt("nGenPar");
			Int_t* genParId      = data.GetPtrInt("genParId");
			Int_t* genParSt      = data.GetPtrInt("genParSt");
			Int_t* genMomParId   = data.GetPtrInt("genMomParId");
			Int_t* genDa1      = data.GetPtrInt("genDa1");
			Int_t* genDa2      = data.GetPtrInt("genDa2");

			int genHIndex[2]={-1,-1};
			int genbIndex[2][2]={{-1,-1},
						{-1,-1}};		       

			for(int ig=0; ig < nGenPar; ig++){

				if(genParId[ig]!=25)continue;

				if(genHIndex[0]<0)
				{
				genHIndex[0]=ig;
				genbIndex[0][0]=genDa1[ig];
				genbIndex[0][1]=genDa2[ig];
				}

				else if(genHIndex[1]<0)
				{
				genHIndex[1]=ig;
				genbIndex[1][0]=genDa1[ig];
				genbIndex[1][1]=genDa2[ig];
				}

			}    

			if(genHIndex[0]<0 || genHIndex[1]<0)continue;
			if(genbIndex[0][0]<0 || genbIndex[0][1]<0)continue;
			if(genbIndex[1][0]<0 || genbIndex[1][1]<0)continue;

			nPass[0]++;

			if(genHIndex[0]==genHIndex[1])continue;
			nPass[1]++;

			TLorentzVector genH_l4[2];
			TLorentzVector genb_l4[2][2];
			TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
			
			

			for(int ih=0; ih<2; ih++)
			{
				genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));
				for(int ib=0; ib<2; ib++)
				{
				genb_l4[ih][ib] = *((TLorentzVector*)genParP4->At(genbIndex[ih][ib]));
				}
				}


			
			
			TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");

			// check matching first    
			bool findAK8Match=false;
			const float dRMax=0.4;
			const float dRbMax=0.8;
			int matchedHAK8JetIndex[2]={-1,-1};
		      int AK8nJet=data.GetInt("AK8PuppinJet");
			if(AK8nJet<2)continue;
			bool matchb=1;
			for(int ij=0; ij<AK8nJet; ij++)
			{
				TLorentzVector* thisJet = (TLorentzVector*)AK8PuppijetP4->At(ij);

				for(int jj=0; jj<AK8nJet; jj++)
				{

				if(ij==jj)continue;
				TLorentzVector* thatJet = (TLorentzVector*)AK8PuppijetP4->At(jj);
	    
				if(thisJet->DeltaR(genH_l4[0])<dRMax && 
					(!matchb || (matchb && 
						thisJet->DeltaR(genb_l4[0][0])<dRbMax && 
						thisJet->DeltaR(genb_l4[0][1])<dRbMax)) &&
					thatJet->DeltaR(genH_l4[1])<dRMax &&
					(!matchb || (matchb &&
						thatJet->DeltaR(genb_l4[1][0])<dRbMax &&
						thatJet->DeltaR(genb_l4[1][1])<dRbMax)))
					{
					
					if(ij<jj){
					matchedHAK8JetIndex[0]=ij;
					matchedHAK8JetIndex[1]=jj;
					}
					else
					{
					matchedHAK8JetIndex[0]=jj;
					matchedHAK8JetIndex[1]=ij;
					}
					findAK8Match=true;
					break;
					}

				if(findAK8Match)break;
	
				}	

				if(findAK8Match)break;

				}
	
			if(!findAK8Match)continue;
			
			
				
			
			
			Float_t*  AK8PuppijetGenSDmass = data.GetPtrFloat("AK8PuppijetGenSDmass");
			Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");
			int* AK8PuppinSubSDJet=data.GetPtrInt("AK8PuppinSubSDJet");
			
		
			Int_t* AK8Puppijet_nSV=data.GetPtrInt("AK8Puppijet_nSV");
			vector<float>   *AK8Puppijet_SVMass  =  data.GetPtrVectorFloat("AK8Puppijet_SVMass");
			int nEle= data.GetInt("nEle");
			int nMu=data.GetInt("nMu");
			Float_t*  AK8PuppijetEleEF = data.GetPtrFloat("AK8PuppijetEleEF");
			Float_t*  AK8PuppijetMuoEF = data.GetPtrFloat("AK8PuppijetMuoEF");
			Int_t* AK8PuppijetCMulti=data.GetPtrInt("AK8PuppijetCMulti");
			Int_t* AK8PuppijetEleMulti=data.GetPtrInt("AK8PuppijetEleMulti");
			Int_t* AK8PuppijetMuoMulti=data.GetPtrInt("AK8PuppijetMuoMulti");
			Int_t nVtx        = data.GetInt("nVtx");
			for(int i=0; i<2;i++){
		
				
				int AK8jet=matchedHAK8JetIndex[i];
				
				
				
				TLorentzVector* thisAK8Jet = (TLorentzVector*)AK8PuppijetP4->At(AK8jet);
				
				if(thisAK8Jet->Pt()<200)continue;
				if(fabs(thisAK8Jet->Eta())>2.4)continue;
				if(AK8PuppinSubSDJet[AK8jet]!=2)continue;
				
				
				pt_AK8MatchedToHbb=thisAK8Jet->Pt();
				eta_AK8MatchedToHbb=thisAK8Jet->Eta();
				nsv_AK8MatchedToHbb=AK8Puppijet_nSV[AK8jet];
				sv0mass_AK8MatchedToHbb=AK8Puppijet_SVMass[AK8jet][0];
				sv1mass_AK8MatchedToHbb=AK8Puppijet_SVMass[AK8jet][1];
				nmu_AK8MatchedToHbb=AK8PuppijetMuoMulti[AK8jet];
				nel_AK8MatchedToHbb=AK8PuppijetEleMulti[AK8jet];
				muenfr_AK8MatchedToHbb=AK8PuppijetMuoEF[AK8jet];
				nch_AK8MatchedToHbb=AK8PuppijetCMulti[AK8jet];
				emenfr_AK8MatchedToHbb=AK8PuppijetEleEF[AK8jet];
				spec1=nVtx;
				spec2=AK8PuppijetSDmass[AK8jet];
				Float_t val ;
				for (Int_t ih=0; ih<nhists; ih++) {
				TString title = hists[ih]->GetTitle();
				val= (reader->EvaluateRegression( title ))[0];
				}
				
				
				TLorentzVector test0= (*((TLorentzVector*)AK8PuppijetP4->At(AK8jet)))*(val );
				thisAK8Jet=&test0;
				if(thisAK8Jet->Pt()<200)continue;
				if(fabs(thisAK8Jet->Eta())>2.4)continue;
				th1->Fill(AK8PuppijetGenSDmass[AK8jet]*val);
				
				for(int j=0;j<13;j++){
					if(thisAK8Jet->Pt()>ptBins[j] && thisAK8Jet->Pt()<ptBins[j+1]){
						
						if(fabs(thisAK8Jet->Eta())<1.3){
							th2[0][j]->Fill(AK8PuppijetGenSDmass[AK8jet]*val);
							th2[2][j]->Fill(AK8PuppijetSDmass[AK8jet]*val);
							th2[4][j]->Fill(thisAK8Jet->Pt());
						}
						else{
							th2[1][j]->Fill(AK8PuppijetGenSDmass[AK8jet]*val);
							th2[3][j]->Fill(AK8PuppijetSDmass[AK8jet]*val);
							th2[5][j]->Fill(thisAK8Jet->Pt());
						}
							
					}
				}
				
				if(thisAK8Jet->Pt()>ptBins[13]){
					if(fabs(thisAK8Jet->Eta())<1.3){
							th2[0][13]->Fill(AK8PuppijetGenSDmass[AK8jet]*val);
							th2[2][13]->Fill(AK8PuppijetSDmass[AK8jet]*val);
							th2[4][13]->Fill(thisAK8Jet->Pt());
						}
						else{
							th2[1][13]->Fill(AK8PuppijetGenSDmass[AK8jet]*val);
							th2[3][13]->Fill(AK8PuppijetSDmass[AK8jet]*val);
							th2[5][13]->Fill(thisAK8Jet->Pt());
						}
				}
				
			}
			
		}
	}	
	cout<<"entries="<<total<<endl;	
	
	
	TFile* outFile ;
	outFile= new TFile(Form("corr2/%s.root",st2.data()),"recreate");
	th1->Write();
	for(int i=0;i<14;i++){
		th2[0][i]->Write();
		th2[1][i]->Write();
		th2[2][i]->Write();
		th2[3][i]->Write();
		th2[4][i]->Write();
		th2[5][i]->Write();
	}
	outFile->Close();
  
		
	
	
	
  
   delete reader;
    
   
}

void HH4bReg(){

		string st1[40]={
		/*0-11*/
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1200.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1400.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1600.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B1800.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B2000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B2500.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B3000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B3500.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B4000.root",
		"/data7/chchen/AK8PuppijetGenSDmass/B/B4500.root",
		/*12-21*/
		"/data7/chchen/80x_dsv/R/R1000.root",
		"/data7/chchen/80x_dsv/R/R1200.root",
		"/data7/chchen/80x_dsv/R/R1400.root",
		"/data7/chchen/80x_dsv/R/R1600.root",
		"/data7/chchen/80x_dsv/R/R1800.root",
		"/data7/chchen/80x_dsv/R/R2000.root",
		"/data7/chchen/80x_dsv/R/R2500.root",
		"/data7/chchen/80x_dsv/R/R3000.root",
		"/data7/chchen/80x_dsv/R/R3500.root",
		"/data7/chchen/80x_dsv/R/R4000.root",
		"/data7/chchen/80x_dsv/R/R4500.root",
		/*22-32*/
		"/data7/chchen/80x_dsv/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCD/2000_1/NCUGlobalTuples_",	
		
		"/data7/chchen/80x_dsv/QCDHTbGen/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbGen/2000/NCUGlobalTuples_",
		
		"/data7/chchen/80x_dsv/QCDHTbEnriched/700/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/1000/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/1500/NCUGlobalTuples_",
		"/data7/chchen/80x_dsv/QCDHTbEnriched/2000/NCUGlobalTuples_",
	};
	string  fileName[40]={
	"B1000","B1200","B1400","B1600","B1800","B2000","B2500","B3000","B3500","B4000","B4500",
	"R1000","R1200","R1400","R1600","R1800","R2000","R2500","R3000","R3500","R4000","R4500",
	"QCD700_1","QCD700_2","QCD1000_1","QCD1000_2","QCD1500_1","QCD1500_2","QCD2000_1","QCD2000_2",
	"bGen700","bGen1000","bGen1500","bGen2000",
	"bEnriched700","bEnriched1000","bEnriched1500","bEnriched2000",
	};
	int aa[40]={119,225,40,80,31,62,17,33,27,12,5,3,9,3,2,3};
	string dataPathB="/data7/chchen/80x_dsv/JetHT_runB/NCUGlobalTuples_";
	string dataPathC="/data7/chchen/80x_dsv/JetHT_runC/NCUGlobalTuples_";
	string dataPathD="/data7/chchen/80x_dsv/JetHT_runD/NCUGlobalTuples_";
	
	for(int i=0;i<11;i++)TMVARegressionApplication(1,2,st1[i],fileName[i]);
}