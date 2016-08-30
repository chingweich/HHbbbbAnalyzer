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

void makeCategory(){
	
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("B%s.root",masspoint[i].data()));
		f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	f[22]=TFile::Open("QCD700.root");
	f[23]=TFile::Open("QCD1000.root");
	f[24]=TFile::Open("QCD1500.root");
	f[25]=TFile::Open("QCD2000.root");
	//f[11]=TFile::Open("data1.root");
	
	TH1D* th1[3][26];
	
	//TFile* output=new TFile("MassPlotFineBins_subtr_Moriond_Silver.root","recreate");
	for (int i=0;i<3;i++){
		for(int j=0;j<26;j++){
			//if(j>11)continue;
			//cout<<i<<","<<j<<endl;
			th1[i][j]=(TH1D*)f[j]->FindObjectAny(Form("cat%d",i));
			//cout<<i<<","<<j<<endl;
			if(j<11)th1[i][j]->SetName(Form("Graviton%s_cat%d",masspoint[j].data(),i));
			else if(j<22)th1[i][j]->SetName(Form("Radion%s_cat%d",masspoint[j-11].data(),i));
			else th1[i][j]->SetName(Form("Data_cat%d",i));
			//th1[i][j]->Write();
			//cout<<th1[i][j]->Integral()<<endl;
		}
	}
	
	double fixNumber=22295/30584.3;//7306/11857.3;
	double Xsec[4]={6831,1207,119.9,25.24};
	
	for (int i=0;i<3;i++){
		for(int j=0;j<4;j++){
			TH1D *th2=(TH1D* )f[j+22]->FindObjectAny("fixScale");
			th1[i][j+22]->Scale(fixNumber*12883.846147301*Xsec[j]/th2->GetBinContent(1));
			if(j==3){
				th1[i][j+22]->Add(th1[i][22]);
				th1[i][j+22]->Add(th1[i][23]);
				th1[i][j+22]->Add(th1[i][24]);
			}
		}
	}
	
	//for(int i=0;i<12;i++)f[i]->Close();
	TFile* output=new TFile("MassPlotFineBins_subtr_Moriond_Silver.root","recreate");
	for (int i=0;i<3;i++){
		for(int j=0;j<26;j++){
			cout<<i<<","<<j<<","<<th1[i][j]->Integral()<<endl;
			if(j>21)continue;
			th1[i][j]->Write();
			//if(i==2 && j==11)break;
			//cout<<i<<","<<j<<endl;
		}
		th1[i][25]->Write();
	}
	/*for(int j=0;j<12;j++){
		cout<<j<<endl;
		th1[0][j]->Write();
		th1[1][j]->Write();
		th1[2][j]->Write();
	}
	*/
	//th1[i][j]->Write();
	output->Close();
	
}
