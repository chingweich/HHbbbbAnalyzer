//Histogram 
#include <TH1.h>
//#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TF1.h"
#include <TGraph.h>
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

float getPUPPIweight(float puppipt, float puppieta ){

double barrleMean[13]={1.16593,1.14661,1.14849,1.14509,1.14044,1.14157,1.14284,1.14445,1.14922,1.15536,1.16216,1.17066,1.16181};
double endcapMean[13]={1.32643,1.24548,1.2552,1.24958,1.24539,1.23287,1.21979,1.22366,1.22922,1.22335,1.22981,1.23398,1.23844};
double barrlex[13]={325.254,364.602,407.173,452.286,540.273,627.553,714.669,799.949,883.213,1097.31,1310.44,1728.48,1937.44};
double endcapx[13]={300.721,337.376,374.014,408.777,472.758,531.384,589.311,637.263,683.63,804.23,910.163,1081.75,1176.51};
TGraph* tg1=new TGraph(13,barrlex,barrleMean);
TGraph* tg2=new TGraph(13,endcapx,endcapMean);

  float totalWeight=1;
   if( fabs(puppieta)  <= 1.3 ){
	return tg1->Eval(puppipt);
   }
else{
	 return tg2->Eval(puppipt);
}
  // file->Close();
  //return totalWeight;
}

float getPUPPIweightOnRegressed(float puppipt, float puppieta ){

double barrleMean[13]={0.983602,1.01681,1.0416,1.06712,1.09811,1.09607,1.05727,1.07846,1.10783,1.09147,1.12839};
double endcapMean[13]={1.09957,1.1156,1.12732,1.1481,1.16446,1.15879,1.14886,1.15286,1.16249,1.17772,1.17944};
double barrlex[13]={382.737,416.886,450.274,487.441,568.772,665.847,761.498,842.68,919.298,1152.37,1345.38};
double endcapx[13]={355.246,390.192,421.469,449.83,508.275,567.493,624.998,673.756,715.198,830.788,933.179};
TGraph* tg1=new TGraph(13,barrlex,barrleMean);
TGraph* tg2=new TGraph(13,endcapx,endcapMean);

  float totalWeight=1;
   if( fabs(puppieta)  <= 1.3 ){
	return tg1->Eval(puppipt);
   }
else{
	 return tg2->Eval(puppipt);
}
  // file->Close();
  //return totalWeight;
}

float getPUPPIweight_o(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters(
   				      1.26638,
   				      -0.000658496,
   				      9.73779e-07,
   				      -5.93843e-10,
   				      1.61619e-13,
   				      -1.6272e-17);

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else
    if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );

  totalWeight = genCorr * recoCorr;
  // file->Close();
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
	
	TH1D* th2d[8];
	th2d[0]=new TH1D("1a","1a",4000,1000,5000);	
	th2d[1]=new TH1D("2a","2a",4000,1000,5000);	
	th2d[2]=new TH1D("2b","2b",4000,1000,5000);	
	
	
	th2d[3]=new TH1D("1aL","1aL",4000,1000,5000);	
	th2d[4]=new TH1D("2aL","2aL",4000,1000,5000);	
	th2d[5]=new TH1D("2bL","2bL",4000,1000,5000);	
		
	
		
	
	//int nWidth=5,nBmin=11;
	 int width [nWidth]={20,25,30,35,40};
	 int bmin[nBmin]={95,100,105,110,115,120,125,130,135,140,145};
	 
	 TH1D* th3d[8][nWidth][nBmin];
	 
	 for(int i=0;i<nWidth;i++){
		 for(int j=0;j<nBmin;j++){
			 th3d[0][i][j]=(TH1D*) th2d[0]->Clone(Form("%s_%d_%d",th2d[0]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[1][i][j]=(TH1D*) th2d[1]->Clone(Form("%s_%d_%d",th2d[1]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[2][i][j]=(TH1D*) th2d[2]->Clone(Form("%s_%d_%d",th2d[2]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 
			 th3d[3][i][j]=(TH1D*) th2d[3]->Clone(Form("%s_%d_%d",th2d[3]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[4][i][j]=(TH1D*) th2d[4]->Clone(Form("%s_%d_%d",th2d[4]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 th3d[5][i][j]=(TH1D*) th2d[5]->Clone(Form("%s_%d_%d",th2d[5]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 //th3d[6][i][j]=(TH1D*) th2d[6]->Clone(Form("%s_%d_%d",th2d[6]->GetTitle(),bmin[j],width[i]+bmin[j]));
			 //th3d[7][i][j]=(TH1D*) th2d[7]->Clone(Form("%s_%d_%d",th2d[7]->GetTitle(),bmin[j],width[i]+bmin[j]));
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

			const int nAK8Jet=data.GetInt("AK8PuppinJet");
			//2.nJets
			if(nAK8Jet<2)continue;nPass[2]++;
			int* AK8PuppinSubSDJet=data.GetPtrInt("AK8PuppinSubSDJet");
			if(AK8PuppinSubSDJet[0]!=2||AK8PuppinSubSDJet[1]!=2)continue;
			TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");
			float*  AK8PuppijetCorrUncUp = data.GetPtrFloat("AK8PuppijetCorrUncUp"); 
			float*  AK8PuppijetCorrUncDown = data.GetPtrFloat("AK8PuppijetCorrUncDown"); 
			TLorentzVector* thisJet ,* thatJet;
			/*
			vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPx");
			vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPy");
			vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("AK8PuppisubjetSDPz");
			vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("AK8PuppisubjetSDE");
			TLorentzVector* subjetP4[2][2];
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					subjetP4[i][j]=new TLorentzVector(0,0,0,0);
					subjetP4[i][j]->SetPxPyPzE(subjetSDPx[i][j],subjetSDPy[i][j],subjetSDPz[i][j],subjetSDE[i][j]);
				}
			}
			*/
			
			if(JESOption==0){
				thisJet=(TLorentzVector*)AK8PuppijetP4->At(0);
				thatJet=(TLorentzVector*)AK8PuppijetP4->At(1);
			}
			else if (JESOption==1){
				TLorentzVector thisJetScaleUp= (*((TLorentzVector*)AK8PuppijetP4->At(0)))*(1+AK8PuppijetCorrUncUp[0]);
				TLorentzVector thatJetScaleUp= (*((TLorentzVector*)AK8PuppijetP4->At(1)))*(1+AK8PuppijetCorrUncUp[1]);
				thisJet= &thisJetScaleUp;
				thatJet= &thatJetScaleUp;
			}
			else if (JESOption==2){
				TLorentzVector thisJetScaleDown= (*((TLorentzVector*)AK8PuppijetP4->At(0)))*(1-AK8PuppijetCorrUncDown[0]);
				TLorentzVector thatJetScaleDown= (*((TLorentzVector*)AK8PuppijetP4->At(1)))*(1-AK8PuppijetCorrUncDown[1]);
				thisJet= &thisJetScaleDown;
				thatJet= &thatJetScaleDown;
			}
			
			//3. Pt 
			if(thisJet->Pt()>99998 ||thatJet->Pt()>99998 )continue;
			if(thisJet->Pt()<300)continue;
			if(thatJet->Pt()<300)continue;
			nPass[3]++;
			//4tightId-----------------------------------------
			vector<bool>    &AK8PuppijetPassIDTight = *((vector<bool>*) data.GetPtr("AK8PuppijetPassIDTight"));
			if(AK8PuppijetPassIDTight[0]==0)continue;
			if(AK8PuppijetPassIDTight[1]==0)continue;
			Float_t*  AK8PuppijetCEmEF = data.GetPtrFloat("AK8PuppijetCEmEF");
			Float_t*  AK8PuppijetMuoEF = data.GetPtrFloat("AK8PuppijetMuoEF");
			if(AK8PuppijetMuoEF[0]>0.8)continue;
			if(AK8PuppijetCEmEF[0]>0.9)continue;
			if(AK8PuppijetMuoEF[1]>0.8)continue;
			if(AK8PuppijetCEmEF[1]>0.9)continue;
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
			//float mjjRed = (*thisJet+*thatJet).M()+250-thisJet->M()-thatJet->M();
			//if(mjjRed<1000)continue;
			nPass[7]++;
			//8. fatjetPRmassL2L3Corr-----------------------------------------
			nPass[8]++;
			//9.-----------------------------------------
		
			
    
			Float_t*  AK8Puppijet_DoubleSV = data.GetPtrFloat("AK8Puppijet_DoubleSV");
			
			
			int nDSVTight=0,nDSVMed=0,nDSVLoose=0;
			if(AK8Puppijet_DoubleSV[0]>0.8)nDSVTight++;
			if(AK8Puppijet_DoubleSV[0]>0.6)nDSVMed++;
			if(AK8Puppijet_DoubleSV[0]>0.3)nDSVLoose++;
			
			if(AK8Puppijet_DoubleSV[1]>0.8)nDSVTight++;
			if(AK8Puppijet_DoubleSV[1]>0.6)nDSVMed++;
			if(AK8Puppijet_DoubleSV[1]>0.3)nDSVLoose++;
			
			double varTemp[2];
			
			Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");
			
			Int_t* AK8Puppijet_nSV=data.GetPtrInt("AK8Puppijet_nSV");
			vector<float>   *AK8Puppijet_SVMass  =  data.GetPtrVectorFloat("AK8Puppijet_SVMass");
			int nEle= data.GetInt("nEle");
			int nMu=data.GetInt("nMu");
			Float_t*  AK8PuppijetEleEF = data.GetPtrFloat("AK8PuppijetEleEF");
			//Float_t*  AK8PuppijetMuoEF = data.GetPtrFloat("AK8PuppijetMuoEF");
			Int_t* AK8PuppijetCMulti=data.GetPtrInt("AK8PuppijetCMulti");
			Int_t* AK8PuppijetEleMulti=data.GetPtrInt("AK8PuppijetEleMulti");
			Int_t* AK8PuppijetMuoMulti=data.GetPtrInt("AK8PuppijetMuoMulti");
			
			for(int i=0; i<2;i++){
		
				TLorentzVector* thisAK8Jet ;
				
				if(i==1)thisAK8Jet=thatJet;
				else thisAK8Jet=thisJet;
				
				
				pt_AK8MatchedToHbb=thisAK8Jet->Pt();
				eta_AK8MatchedToHbb=thisAK8Jet->Eta();
				nsv_AK8MatchedToHbb=AK8Puppijet_nSV[i];
				sv0mass_AK8MatchedToHbb=AK8Puppijet_SVMass[i][0];
				sv1mass_AK8MatchedToHbb=AK8Puppijet_SVMass[i][1];
				nmu_AK8MatchedToHbb=AK8PuppijetMuoMulti[i];
				nel_AK8MatchedToHbb=AK8PuppijetEleMulti[i];
				muenfr_AK8MatchedToHbb=AK8PuppijetMuoEF[i];
				nch_AK8MatchedToHbb=AK8PuppijetCMulti[i];
				emenfr_AK8MatchedToHbb=AK8PuppijetEleEF[i];
				spec1=nVtx;
				spec2=AK8PuppijetSDmass[i];
				Float_t val ;
				for (Int_t ih=0; ih<nhists; ih++) {
				TString title = hists[ih]->GetTitle();
				val= (reader->EvaluateRegression( title ))[0];
				}
				varTemp[i]=val;
			}
			
			
			double PUPPIweight[2]={0};
			PUPPIweight[0]=getPUPPIweight(thisJet->Pt(),thisJet->Eta());
			PUPPIweight[1]=getPUPPIweight(thatJet->Pt(),thatJet->Eta());
			
			
	
			double Mjja= ((*thisJet)+(*thatJet)).M()+250
									-((*thisJet)).M()-((*thatJet)).M();
									
			TLorentzVector  thisJetReg, thatJetReg;
			thisJetReg=(*thisJet)*varTemp[0];
			thatJetReg=(*thatJet)*varTemp[1];
			
			double Mjjb= (thisJetReg+thatJetReg).M()+250
									-(thisJetReg).M()-(thatJetReg).M();
			
			double PUPPIweightOnRegressed[2]={0};			
			PUPPIweightOnRegressed[0]=getPUPPIweightOnRegressed(thisJetReg.Pt(),thisJetReg.Eta());
			PUPPIweightOnRegressed[1]=getPUPPIweightOnRegressed(thatJetReg.Pt(),thatJetReg.Eta());
			
			TLorentzVector  thisJetRegHCorr, thatJetRegHCorr;
			thisJetRegHCorr=(thisJetReg)*PUPPIweightOnRegressed[0];
			thatJetRegHCorr=(thatJetReg)*PUPPIweightOnRegressed[1];
			double Mjjcase3= (thisJetRegHCorr+thatJetRegHCorr).M()+250
									-(thisJetRegHCorr).M()-(thatJetRegHCorr).M();
	
			 for(int i=0;i<nWidth;i++){
				for(int j=0;j<nBmin;j++){
					if(width[i]+bmin[j]>166)continue;
					
					if(AK8PuppijetSDmass[0]*PUPPIweight[0]>bmin[j]&& AK8PuppijetSDmass[0]*PUPPIweight[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*PUPPIweight[1]>bmin[j] && AK8PuppijetSDmass[1]*PUPPIweight[1]<width[i]+bmin[j] && nDSVTight==2) 
					th3d[0][i][j]->Fill(Mjja,PU_weight[0]);
					
					if(AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]>bmin[j]&& AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]>bmin[j] && AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]<width[i]+bmin[j]&& nDSVTight==2)  
					th3d[1][i][j]->Fill(Mjja,PU_weight[0]);
					
					if(AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]>bmin[j]&& AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]>bmin[j] && AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]<width[i]+bmin[j] && nDSVTight==2)  
					th3d[2][i][j]->Fill(Mjjb,PU_weight[0]);
					
					
					if(AK8PuppijetSDmass[0]*PUPPIweight[0]>bmin[j]&& AK8PuppijetSDmass[0]*PUPPIweight[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*PUPPIweight[1]>bmin[j] && AK8PuppijetSDmass[1]*PUPPIweight[1]<width[i]+bmin[j]&& (nDSVLoose==2 && !(nDSVTight==2))) 
					th3d[3][i][j]->Fill(Mjja,PU_weight[0]);
					
					if(AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]>bmin[j]&& AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]>bmin[j] && AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]<width[i]+bmin[j]  && (nDSVLoose==2 && !(nDSVTight==2)))  
					th3d[4][i][j]->Fill(Mjja,PU_weight[0]);
					
					if(AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]>bmin[j]&& AK8PuppijetSDmass[0]*varTemp[0]*PUPPIweightOnRegressed[0]<width[i]+bmin[j]
					&&AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]>bmin[j] && AK8PuppijetSDmass[1]*varTemp[1]*PUPPIweightOnRegressed[1]<width[i]+bmin[j] && (nDSVLoose==2 && !(nDSVTight==2)))  
					th3d[5][i][j]->Fill(Mjjb,PU_weight[0]);
					
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
			//th3d[6][i][j]->Write();
			//th3d[7][i][j]->Write();
		}
	}
	outFile->Close();
  
   delete reader;
    
   
}

void HH4bRegCategory(int a){

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
		"/data7/chchen/8_0_20/QCD/700/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/700_1/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/1000/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/1000_1/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/1500/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/1500_1/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/2000/NCUGlobalTuples_",
		"/data7/chchen/8_0_20/QCD/2000_1/NCUGlobalTuples_",	
		
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