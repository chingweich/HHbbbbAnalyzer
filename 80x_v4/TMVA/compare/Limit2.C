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
#include"../../../setNCUStyle.C"

TCanvas* c1;

#define  nWidth 5
#define  nBmin 11

void Limit2(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",800,600);
	
	int width [nWidth]={20,25,30,35,40};
	int bmin[nBmin]={95,100,105,110,115,120,125,130,135,140,145};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	TLegend *leg = new TLegend(0.30, 0.68, 0.96, 0.92);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	
			 
			 TFile* tf1[3];
			 tf1[0]=TFile::Open(Form("../pr/massOpt/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",100,135));
			 tf1[1]=TFile::Open(Form("../sd/massOpt/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",105,145));
			 tf1[2]=TFile::Open(Form("../sd2/massOpt/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",100,135));
			 //if (!tf1 || !tf1->IsOpen())continue;
			 TGraphAsymmErrors* tg1[3];
			 tg1[0]=(TGraphAsymmErrors*)tf1[0]->Get("LimitExpectedCLs");
			 tg1[1]=(TGraphAsymmErrors*)tf1[1]->Get("LimitExpectedCLs");
			 tg1[2]=(TGraphAsymmErrors*)tf1[2]->Get("LimitExpectedCLs");
			 
			 tg1[0]->GetYaxis()->SetTitle("95% CLs on #sigma(X#rightarrowHH)#timesBR(HH#rightarrowb#bar{b}b#bar{b})[fb]");
			 tg1[0]->SetLineStyle(1);
			 tg1[1]->SetLineStyle(1);
			 tg1[2]->SetLineStyle(1);
			 tg1[0]->SetFillColor(0);
			 tg1[1]->SetFillColor(0);
			 tg1[2]->SetFillColor(0);
			 tg1[0]->SetLineColor(1);
			 tg1[1]->SetLineColor(2);
			 tg1[2]->SetLineColor(3);
			 leg->AddEntry(tg1[0],Form("L2L3pruned%dto%d",100,135));
			 leg->AddEntry(tg1[1],Form("PUPPISDmass+PUPPIjet%dto%d",105,145));
			 leg->AddEntry(tg1[2],Form("PUPPISDmass+PUPPISDjet%dto%d",100,135));
			 tg1[0]->SetMaximum(15);
			// cout<<k<<","<<m<<endl;
			tg1[0]->Draw("APL");
			tg1[1]->Draw("PL same");
			tg1[2]->Draw("PL same");
	
	leg->Draw("same");
	c1->Print("Limit2.pdf");
}
