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
#include "../setNCUStyle.C"

TCanvas* c1;

void draw(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("sf/B%s.root",masspoint[i].data()));
		f[i+11]=TFile::Open(Form("sf/R%s.root",masspoint[i].data()));
	}
	f[22]=TFile::Open("sf/QCD700.root");
	f[23]=TFile::Open("sf/QCD1000.root");
	f[24]=TFile::Open("sf/QCD1500.root");
	f[25]=TFile::Open("sf/QCD2000.root");
	f[26]=TFile::Open("sf/data.root");
	
	

	
	ofstream myfile;
	myfile.open ("latex_cutflow.txt");
	
	string ss[15]={
		"entries",
		"nVtx$>=$1",
		"triggers",
		"Pt$>$200GeV and tightId",
		"$|\\eta |<$2.4",
		"$|\\Delta \\eta |<$1.3",
		"Mjj$>$1000",
		"105$<$fatjetPRmassL2L3Corr$<$135",
		"$\\tau 21_1$ two LP at least one HP",
		"0-btagged",
		"1-btagged",
		"2-btagged",
		"3-btagged",
		"4-btagged",
		"HPHP3-btagged"
		
	};
	
	TH1D* ct[11];
	for(int i=0;i<11;i++){
		ct[i]=(TH1D*)f[i]->FindObjectAny("Nbtagjet");
	}
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"BulkgravBtagSF 74x selections &1000 & 1200 & 1400 & 1600 & 1800 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int j=0;j<6;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<11;i++){
		
			if (i<10 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	
	for(int i=0;i<11;i++){
		ct[i]=(TH1D*)f[i+11]->FindObjectAny("Nbtagjet");
	}
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"RadionBtagSF 76x selections &1000 & 1200 & 1400 & 1600 & 1800 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int j=0;j<6;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<11;i++){
		
			if (i<10 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	//cout<<"here";
	for(int i=0;i<5;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("Nbtagjet");
	}
	myfile<<""<<endl;
	double Xsec[4]={6831,1207,119.9,25.24};
	for(int i=0;i<4;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("cutflow");
	//	ct[i]->Scale(4513.309419004*Xsec[i]/ct[i]->GetBinContent(1));
	}
	ct[4]=(TH1D*)f[26]->FindObjectAny("cutflow");
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	//myfile<<"Bulkgrav 76x selections &1000 & 1200 & 1400 & 1600 & 1800 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500"<<"\\\\"<<endl;
	myfile<<"QCDHT 76x &700-1000 & 1000-1500 & 1500-2000 & 2000-inf&& QCDAll &&Data\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int j=0;j<15;j++){
		myfile<<ss[j]<<" & ";
		double temp=0;
		for(int i=0;i<5;i++){
			if(i<4)temp+=ct[i]->GetBinContent(j+1);
			if (i<4 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<<temp<<" & "<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	TH1D* cutflow1=new TH1D("cutflowEff","",14,0,14);
	TH1D* cutflow2=new TH1D("cutflowEff","",14,0,14);
	
	/*
	for(int i=0;i<11;i++){
		ct[i]=(TH1D*)f[i+11]->FindObjectAny("cutflow");
	}
	myfile<<endl;
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"Radion selections &1000 & 1200 & 1400 & 1600 & 1800 & 2000 & 2500 & 3000 & 3500 & 4000 & 4500"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<11;i++){
		
			if (i<10 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	myfile<<endl;
	myfile<<""<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"QCD and data selections &QCD 700-1000 & QCD 1000-1500 & QCD 1500-2000 & QCD 2000-inf. &data"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int i=0;i<5;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("cutflow");
	}
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		for(int i=0;i<5;i++){
		
			if (i<4 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< ct[i]->GetBinContent(j+1)<<"\\\\";
		}
		myfile<<endl;
	}
	myfile<<endl;
	myfile<<"\\hline"<<endl;
	myfile<<"QCD normalized selections &QCD 700-1000 & QCD 1000-1500 & QCD 1500-2000 & QCD 2000-inf. & total"<<"\\\\"<<endl;
	myfile<<"\\hline"<<endl;
	for(int i=0;i<4;i++){
		ct[i]=(TH1D*)f[i+22]->FindObjectAny("cutflowS");
	}
	for(int j=0;j<14;j++){
		myfile<<ss[j]<<" & ";
		double sum=0;
		for(int i=0;i<5;i++){
			if (i<4 )sum+=ct[i]->GetBinContent(j+1);
			if (i<4 )myfile<< ct[i]->GetBinContent(j+1)<<" & ";
			else myfile<< sum<<"\\\\";
		}
		myfile<<endl;
	}
	*/
}