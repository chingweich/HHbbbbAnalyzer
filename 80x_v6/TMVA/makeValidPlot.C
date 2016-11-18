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
#include <TF1.h>
#include <string>
#include <sstream>
#include "../../setNCUStyle.C"

#define  nWidth 5
#define  nBmin 11

#define nMass 3
#define nCat 4

TCanvas* c1;

void makeValidPlot(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	TFile* tf1=TFile::Open("MjjVC/data.root");
	
	string validName[7]={
		"0a","0c",
		"1a","1c",
		"2a","2b","2d"
	};
	
	string tau21Name[2]={"withTau21","woTau21"};
	
	for(int i=0;i<7;i++){
		for(int j=0;j<2;j++){
			TH1D* th1,* th2;
			th1= (TH1D*) tf1->Get(Form("fill_%s_105_135_%s",validName[i].data(),tau21Name[j].data()));
			th2= (TH1D*) tf1->Get(Form("valid_%s_105_135_%s",validName[i].data(),tau21Name[j].data()));
			
			th1->Rebin(10);
			th2->Rebin(10);
			th1->GetXaxis()->SetRangeUser(1000,2500);	
			TLegend *leg = new TLegend(0.64, 0.68, 0.92, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  
			leg->AddEntry(th1,"1-antiTag");
			leg->AddEntry(th2,"2-antiTag* PFRatio");
			th1->SetTitle("DBT-TT(Pass_{j0},Fail_{j1})");
			th1->Draw("");
			th2->SetLineColor(2);
			th2->Draw("same");
			leg->Draw("same");
			c1->Print(Form("validPlots/%s/TT_%s.pdf",tau21Name[j].data(),validName[i].data()));
		}
		
	}
	
	for(int i=0;i<7;i++){
		for(int j=0;j<2;j++){
			TH1D* th1,* th2;
			th1= (TH1D*) tf1->Get(Form("fill_%sL_105_135_%s",validName[i].data(),tau21Name[j].data()));
			th2= (TH1D*) tf1->Get(Form("valid_%sL_105_135_%s",validName[i].data(),tau21Name[j].data()));
			
			th1->Rebin(10);
			th2->Rebin(10);
			th1->GetXaxis()->SetRangeUser(1000,2500);	
			TLegend *leg = new TLegend(0.64, 0.68, 0.92, 0.87);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  
			leg->AddEntry(th1,"1-antiTag");
			leg->AddEntry(th2,"2-antiTag* PFRatio");
			th1->SetTitle("DBT-LL(Pass_{j0},Fail_{j1})");
			th1->Draw("");
			th2->SetLineColor(2);
			th2->Draw("same");
			leg->Draw("same");
			c1->Print(Form("validPlots/%s/LL_%s.pdf",tau21Name[j].data(),validName[i].data()));
		}
		
	}

	
}
