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

#define  nWidth 4
#define  nBmin 4

TCanvas* c1;

void makeCategory(){
	
	int width [nWidth]={25,30,35,40};
	int bmin[nBmin]={100,105,110,115};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("MjjVC/B%s.root",masspoint[i].data()));
		f[i+11]=TFile::Open(Form("MjjVC/B%s.root",masspoint[i].data()));
	}
	f[22]=TFile::Open("MjjVC/data.root");
	
	for(int k=0;k<nWidth;k++){
		 for(int m=0;m<nBmin;m++){
			 
			 if(width[k]+bmin[m]>166)continue;
	if(width[k]>116)continue;
	
	
	//f[11]=TFile::Open("data1.root");
	
	TH1D* th1[3][26];
	
	
	
	//TFile* output=new TFile("MassPlotFineBins_subtr_Moriond_Silver.root","recreate");
	for (int i=0;i<3;i++){
		for(int j=0;j<23;j++){
			//if(j>11)continue;
			//cout<<i<<","<<j<<endl;
			//th1[i][j]=(TH1D*)f[j]->FindObjectAny(Form("cat%d",i));
			if(i==0 )th1[i][j]=(TH1D*)f[j]->FindObjectAny(Form("0a_%d_%d_woTau21",bmin[m],bmin[m]+width[k]));
			if(i==1 )th1[i][j]=(TH1D*)f[j]->FindObjectAny(Form("0aL_%d_%d_woTau21",bmin[m],bmin[m]+width[k]));
			if(i==2 )th1[i][j]=(TH1D*)th1[0][j]->Clone("copy");
			
			//cout<<i<<","<<j<<endl;
			if(j<11)th1[i][j]->SetName(Form("Graviton%s_cat%d",masspoint[j].data(),i));
			else if(j<22)th1[i][j]->SetName(Form("Radion%s_cat%d",masspoint[j-11].data(),i));
			else th1[i][j]->SetName(Form("QCD_cat%d",i));
			//th1[i][j]->Write();
			//cout<<th1[i][j]->Integral()<<endl;
		}
	}
	
	TFile* output=new TFile(Form("mass/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",bmin[m],bmin[m]+width[k]),"recreate");
	for (int i=0;i<3;i++){
		for(int j=0;j<23;j++){
			cout<<i<<","<<j<<","<<th1[i][j]->Integral()<<endl;
			//if(j>21)continue;
			th1[i][j]->Write();
			//if(i==2 && j==11)break;
			//cout<<i<<","<<j<<endl;
		}
		//th1[i][25]->Write();
	}
	
	output->Close();
	}}
}
