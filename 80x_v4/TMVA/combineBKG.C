#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>

TCanvas* c1;

void combineBKG(){
	
	TFile *f[40];
	
	f[0]=TFile::Open("QCD700.root");
	f[1]=TFile::Open("QCD1000.root");
	f[2]=TFile::Open("QCD1500.root");
	f[3]=TFile::Open("QCD2000.root");
	
	
	
	vector<string> h_name;
	h_name.push_back("regMass_j0");
	h_name.push_back("regMass_j1");
	h_name.push_back("prMass_j0");
	h_name.push_back("prMass_j1");
	h_name.push_back("unprMass_j0");
	h_name.push_back("unprMass_j1");
	h_name.push_back("reg2Mass_j0");
	h_name.push_back("reg2Mass_j1");
	h_name.push_back("SD_j0");
	h_name.push_back("SD_j1");
	
	TH1D* th1[4];
	double Xsec[4]={6831,1207,119.9,25.24};
	double fixNumber=1.02879/1.13007;
	
	
	TFile* output=new TFile("opt.root","recreate");
	
	for(unsigned int k=0;k<h_name.size();k++){
		for(int i=0;i<4;i++){
			th1[i]=(TH1D*)f[j]->FindObjectAny(Form("%s",h_name[k].data()));
			
			TH1D *th2=(TH1D* )f[i]->FindObjectAny("fixScale");
			th1[i]->Scale(fixNumber*12883.846147301*Xsec[i]/th2->GetBinContent(1));
		}
		th1[0]->Add(th1[1]);
		th1[0]->Add(th1[2]);
		th1[0]->Add(th1[3]);
		
		th1[0]->Write();
	}
	
	output->Close();
	
}
