//Histogram 
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
//#include <TH2D.h>
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
#include <TSystem.h>
#include <TF1.h>
//math 
//#include <cmath>
//#include <algorithm>

//other including
//#include "setNCUStyle.C"
#include "untuplizer.h"
//#include "ApplyJER.h"
//#include "jetEnergyScale.h"

//#include "standalone_LumiReWeighting.cc"
//#include "standalone_LumiReWeighting.h"
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


int skimTree(int wMs,string st) {
	TFile *f;
	TTree *tree;
	int nPass[20]={0},total=0;
	
	vector<string> str;
	vector<double> up;
	vector<double> down;
	str.push_back("jet1pt"); up.push_back(800);  down.push_back(100);  
	str.push_back("jet2pt"); up.push_back(800);  down.push_back(100);  
	str.push_back("jet1eta"); up.push_back(3);  down.push_back(-3);  
	str.push_back("jet2eta"); up.push_back(3);  down.push_back(-3);  
	
	str.push_back("unjet1pt"); up.push_back(800);  down.push_back(100);  
	str.push_back("unjet2pt"); up.push_back(800);  down.push_back(100);  
	str.push_back("unjet1eta"); up.push_back(3);  down.push_back(-3);  
	str.push_back("unjet2eta"); up.push_back(3);  down.push_back(-3);  
	
	str.push_back("etadiff"); up.push_back(1.5);  down.push_back(0);  
	str.push_back("dijetmass"); up.push_back(1500);  down.push_back(700);  
	str.push_back("dijetmass_softdrop_corr"); up.push_back(1500);  down.push_back(700);  
	str.push_back("jet1bbtag"); up.push_back(1);  down.push_back(-1);  
	str.push_back("jet2bbtag"); up.push_back(1);  down.push_back(-1);  
	str.push_back("jet1_puppi_tau21"); up.push_back(1);  down.push_back(0);  
	str.push_back("jet2_puppi_tau21"); up.push_back(1);  down.push_back(0);  
	str.push_back("jet1_puppi_msoftdrop"); up.push_back(150);  down.push_back(40);  
	str.push_back("jet2_puppi_msoftdrop"); up.push_back(150);  down.push_back(40);  
	str.push_back("jet1_puppi_msoftdrop_TheaCorr"); up.push_back(150);  down.push_back(40);  
	str.push_back("jet2_puppi_msoftdrop_TheaCorr"); up.push_back(150);  down.push_back(40);  
	
	TH1D* th1[19];
	for(int i=0;i<19;i++ )th1[i]=new TH1D(Form("%s",str[i].data()),Form("%s",str[i].data()),100,down[i],up[i]);
	//cout<<"here="<<Form("%s%d.root",st.data(),wMs)<<endl;	
	f = TFile::Open(Form("%s",st.data()));
	if (!f || !f->IsOpen())return 0;
	TDirectory * dir;
	dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
		//else dir = (TDirectory*)f->Get(Form("%s:/tree",st.data()));
	dir->GetObject("treeMaker",tree);
		
		//tree=(TTree* )f->Get("treeMaker");
	TreeReader data(tree);
		//total+=data.GetEntriesFast();
	//cout<<"here"<<endl;	
	for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){//event loop----------------------------------------------------------------------------------------
		data.GetEntry(jEntry);
		int nJets         = data.GetInt("FATnJet");
		if(nJets<2)continue;
			
		TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
		TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");
		bool check99999=0;
		for(int ij=0; ij<2; ij++){
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			if(thisJetPuppi->Pt()>99998)check99999=1;
		}
		if(check99999)continue;
		TClonesArray* FATunCorrJetP4   = (TClonesArray*) data.GetPtrTObject("FATunCorrJetP4");
		TLorentzVector* thisFatJet = (TLorentzVector*)fatjetP4->At(0);
		TLorentzVector* thatFatJet = (TLorentzVector*)fatjetP4->At(1);
			
		TLorentzVector* thisUnJet = (TLorentzVector*)FATunCorrJetP4->At(0);
		TLorentzVector* thatUnJet = (TLorentzVector*)FATunCorrJetP4->At(1);
		
		
		if(thatFatJet->Pt()>99998)continue;
		
		th1[0]->Fill(thisFatJet->Pt());
		th1[1]->Fill(thatFatJet->Pt());
		th1[2]->Fill(thisFatJet->Eta());
		th1[3]->Fill(thatFatJet->Eta());
			
		th1[4]->Fill(thisUnJet->Pt());
		th1[5]->Fill(thatUnJet->Pt());
		th1[6]->Fill(thisUnJet->Eta());
		th1[7]->Fill(thatUnJet->Eta());
		
		th1[8]->Fill(fabs(thisFatJet->Eta()-thatFatJet->Eta()));
		
		const float dRMax=0.8;
		int addJetIndex[2]={-1,-1};
		float thea_mass[2]={-1,-1};
		int nADDJets         = data.GetInt("ADDnJet");
		TClonesArray* addjetP4   = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
		Float_t*  addjet_doublesv = data.GetPtrFloat("ADDjet_DoubleSV");
		 Float_t*  FATjetPuppiTau1 = data.GetPtrFloat("FATjetPuppiTau1");
		Float_t*  FATjetPuppiTau2 = data.GetPtrFloat("FATjetPuppiTau2");
		//Float_t*  FATjetPuppiTau3 = data.GetPtrFloat("FATjetPuppiTau3");
		Float_t*  fatjetPuppiSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
		for(int ij=0; ij<2; ij++){
			TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
			TLorentzVector* thisJetPuppi = (TLorentzVector*)puppijetP4->At(ij);
			//if(thisJetPuppi->Pt()>99998)continue;
			float thea_corr = getPUPPIweight(thisJetPuppi->Pt(),thisJetPuppi->Eta());
			float raw_mass = fatjetPuppiSDmass[ij];
			thea_mass[ij] = raw_mass*thea_corr;
			
			for(int k=0; k < nADDJets; k++){
				TLorentzVector* thatJet = (TLorentzVector*)addjetP4->At(k);
				if(thisJet->DeltaR(*thatJet)<dRMax && addJetIndex[ij]<0){
					addJetIndex[ij]=k;
					break;
				}
			}
		}
		
		float mjj = (*thisFatJet+*thatFatJet).M();
		float mjjred = mjj + 250 - thea_mass[0]-thea_mass[1];
		
		th1[9]->Fill(mjj);
		th1[10]->Fill(mjjred);
		th1[11]->Fill(addjet_doublesv[addJetIndex[0]]);
		th1[12]->Fill(addjet_doublesv[addJetIndex[0]]);
		th1[13]->Fill(FATjetPuppiTau2[0]/FATjetPuppiTau1[0]);
		th1[14]->Fill(FATjetPuppiTau2[1]/FATjetPuppiTau1[1]);
		th1[15]->Fill(fatjetPuppiSDmass[0]);
		th1[16]->Fill(fatjetPuppiSDmass[1]);
		th1[17]->Fill(thea_mass[0]);
		th1[18]->Fill(thea_mass[1]);
	}
	
	
	TFile* outFile=new TFile(Form("CUSV/%d.root",wMs),"recreate");
	for(int i=0;i<19;i++ )th1[i]->Write();
	outFile->Close();
	return 0;
}
/*
void checkUnSkimmedVariables(int a,int b){
	//miniTreeBase(1,2,"BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph.root","miniTreeTest");
	//miniTreeBase(1,2,input,output);
	
	
	string st1[40]={
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/B/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/C/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/D/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/E/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/F/NCUGlobalTuples_",
		"/data7/chchen/AK8subjetSDRawFactorNov2016/JetHT/G/NCUGlobalTuples_",
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
	};
	int aa[40]={1045,348,584,497,361,854,24};
	cout<<st1[a]<<b<<".root"<<endl;
	miniTreeBase(b,st1[a],fileName[a]);	
}
*/

  
