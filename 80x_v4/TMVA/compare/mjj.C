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

void mjj(){
	TStyle* ts =setNCUStyle(true);
	ts->SetTitleOffset(2.5, "Y");
	ts->cd();
	
	c1 = new TCanvas("c1","",800,600);
	
	int width [nWidth]={20,25,30,35,40};
	int bmin[nBmin]={95,100,105,110,115,120,125,130,135,140,145};
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	TLegend *leg = new TLegend(0.12, 0.64, 0.8, 0.9);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	
			 
			 TFile* tf1[3];
			 tf1[0]=TFile::Open(Form("../pr/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",105,135));
			 tf1[1]=TFile::Open(Form("../sd/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",105,135));
			 tf1[2]=TFile::Open(Form("../sd2/MassPlotFineBins_subtr_Moriond_Silver%dto%d.root",105,135));
			 //if (!tf1 || !tf1->IsOpen())continue;
			 TH1D* tg1[4];
			 tg1[0]=(TH1D*)tf1[0]->Get("QCD_cat0");
			 tg1[1]=(TH1D*)tf1[0]->Get("Graviton1600_cat0");
			 tg1[2]=(TH1D*)tf1[1]->Get("Graviton1600_cat0");
			 tg1[3]=(TH1D*)tf1[2]->Get("Graviton1600_cat0");
			 
			 for(int i=0;i<4;i++){
tg1[i]->Scale(1/tg1[i]->Integral());
tg1[i]->SetLineColor(i+1);

			 }
			 
			 tg1[0]->SetYTitle("");
			 tg1[0]->SetLineStyle(1);
			 tg1[1]->SetLineStyle(1);
			 tg1[2]->SetLineStyle(1);
			 tg1[0]->SetFillColor(0);
			 tg1[1]->SetFillColor(0);
			 tg1[2]->SetFillColor(0);
			 
			 leg->AddEntry(tg1[0],Form("QCD,L2L3pruned%dto%d",105,135));
			 leg->AddEntry(tg1[1],Form("Graviton1.6TeV,L2L3pruned%dto%d",105,135));
			 leg->AddEntry((TObject*)0,Form("RMS=%f",tg1[1]->GetRMS()),"");
			 leg->AddEntry(tg1[2],Form("Graviton1.6TeV,PUPPISDThea%dto%d,M_{jj}(PUPPIjet)",105,135));
			  leg->AddEntry((TObject*)0,Form("RMS=%f",tg1[2]->GetRMS()),"");
leg->AddEntry(tg1[3],Form("Graviton1.6TeV,PUPPISDThea%dto%d,M_{jj}(PUPPISDjet)",105,135));
 leg->AddEntry((TObject*)0,Form("RMS=%f",tg1[3]->GetRMS()),"");
			 tg1[0]->SetMaximum(0.3);
			// cout<<k<<","<<m<<endl;
			
			for(int i=0;i<4;i++)tg1[i]->Rebin(20);
			tg1[0]->GetXaxis()->SetRangeUser(1000,1800);
			tg1[0]->SetTitle("");
			tg1[0]->SetXTitle("M^{reduced}_{jj}[GeV]");
			
			tg1[0]->Draw("");
			tg1[1]->Draw(" same");
			tg1[2]->Draw(" same");
			tg1[3]->Draw(" same");
			//c1->SetLogy(1);
	
	leg->Draw("same");
	c1->Print("mjj.pdf");
}
