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
#include "../../../setNCUStyle.C"
#include "TF1.h"

# define IntegratedLumi 27000

TCanvas* c1;

void draw(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
	string massName[3]={"raw","regression","HcorrReg"};
	for(int i=0;i<10;i++){
		TFile *f;
	
		f=TFile::Open(Form("B%s.root",masspoint[i].data()));
		TH1D* th1[3];
		
		
		TLegend *leg = new TLegend(0.15, 0.48, 0.46, 0.88);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
		
		for(int j=0;j<3;j++){
				th1[j]=(TH1D*)f->Get(Form("%s_leading",massName[j].data()));
				th1[j]->SetLineColor(j+1);
				th1[j]->GetXaxis()->SetRangeUser(50,150);
				th1[j]->SetTitle(Form("%s",masspoint[i].data()));
				if(j==0)th1[0]->Draw();
				else th1[j]->Draw("same");
				TF1 *tf1[4];
				tf1[0]=new TF1("fa1","gaus(125)",50,150);
				if(j<1)th1[j]->Fit(tf1[0],"","",70,140);
				else th1[j]->Fit(tf1[0],"","",85,150);
				
				
				leg->AddEntry(th1[j],Form("%s",massName[j].data()));
				leg->AddEntry((TObject*)0, Form("Mean=%f",tf1[0]->GetParameter(1)),"");
				leg->AddEntry((TObject*)0, Form("#sigma=%f",tf1[0]->GetParameter(2)),"");
		}
		leg->Draw("same");
		
		if(i==0)c1->Print("leading.pdf(");
		else if(i==9)c1->Print("leading.pdf)");
		else c1->Print("leading.pdf");
	}
	for(int i=0;i<10;i++){
		TFile *f;
	
		f=TFile::Open(Form("B%s.root",masspoint[i].data()));
		TH1D* th1[3];
		
		
	TLegend *leg = new TLegend(0.15, 0.48, 0.46, 0.88);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
		
		for(int j=0;j<3;j++){
				th1[j]=(TH1D*)f->Get(Form("%s_subleading",massName[j].data()));
				th1[j]->SetLineColor(j+1);
				th1[j]->GetXaxis()->SetRangeUser(50,150);
				th1[j]->SetTitle(Form("%s",masspoint[i].data()));
				if(j==0)th1[0]->Draw();
				else th1[j]->Draw("same");
				
				TF1 *tf1[4];
				tf1[0]=new TF1("fa1","gaus(125)",50,150);
				if(j<1)th1[j]->Fit(tf1[0],"","",70,140);
				else th1[j]->Fit(tf1[0],"","",85,150);
				
				leg->AddEntry(th1[j],Form("%s",massName[j].data()));
				leg->AddEntry((TObject*)0, Form("Mean=%f",tf1[0]->GetParameter(1)),"");
				leg->AddEntry((TObject*)0, Form("#sigma=%f",tf1[0]->GetParameter(2)),"");
		}
		leg->Draw("same");
		
		if(i==0)c1->Print("subleading.pdf(");
		else if(i==9)c1->Print("subleading.pdf)");
		else c1->Print("subleading.pdf");
	}
	
	}