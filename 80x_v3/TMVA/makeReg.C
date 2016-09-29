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
#include"../setNCUStyle.C"
#include "TF1.h"
TCanvas* c1;

void makeReg(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	int  masspointy[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	TFile *f[40];
	for(int i=0;i<11;i++){
		f[i]=TFile::Open(Form("category/B%s.root",masspoint[i].data()));
		//f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	//f[22]=TFile::Open("QCD700.root");
	//f[23]=TFile::Open("QCD1000.root");
	//f[24]=TFile::Open("QCD1500.root");
	//f[25]=TFile::Open("QCD2000.root");
	//f[11]=TFile::Open("data1.root");
	
	TH1D* th1[4];
	
	TLegend *leg = new TLegend(0.68, 0.65, 0.94, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  
	for(int i=0;i<11;i++){
		leg->Clear();
		th1[0]=(TH1D*)f[i]->FindObjectAny("regMass_j0");
		th1[3]=(TH1D*)f[i]->FindObjectAny("reg2Mass_j0");
		th1[1]=(TH1D*)f[i]->FindObjectAny("prMass_j0");
		th1[2]=(TH1D*)f[i]->FindObjectAny("unprMass_j0");
		
		TF1 *tf1[4];
		//tf1[0]=new TF1("fa0",Form("gaus(%f,125,10)",th1[0]->Integral()),90,150);
		//tf1[1]=new TF1("fa1",Form("gaus(%f,125,10)",th1[1]->Integral()),90,150);
		//tf1[2]=new TF1("fa2",Form("gaus(%f,125,10)",th1[2]->Integral()),90,150);
		tf1[0]=new TF1("fa2","gaus(1000)",105,135);
		tf1[1]=new TF1("fa2","gaus(1000)",105,135);
		tf1[2]=new TF1("fa2","gaus(1000)",100,130);
		tf1[3]=new TF1("fa2","gaus(1000)",105,135);
		
		
		
		th1[1]->SetLineColor(2);
		th1[2]->SetLineColor(3);
		
		for(int i=0;i<4;i++){
			th1[i]->SetLineColor(i+1);
			tf1[i]->SetLineColor(i+1);
			th1[i]->SetLineWidth(3);
			tf1[i]->SetLineWidth(3);
		}
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		
		th1[0]->Fit(tf1[0],"","",102,132);
		th1[1]->Fit(tf1[1],"","",105,135);
		th1[3]->Fit(tf1[3],"","",105,135);
		th1[2]->Fit(tf1[2],"","",100,130);
		
		leg->AddEntry(th1[2],"uncorrected pruned");
		//leg->AddEntry((TObject*)0,Form("rms/mean=%f",th1[2]->GetRMS()/th1[2]->GetMean()));
		leg->AddEntry((TObject*)0,Form("#sigma/mean=%f",tf1[2]->GetParameter(2)/tf1[2]->GetParameter(1)));
		
		leg->AddEntry(th1[0],"regression");
		//leg->AddEntry((TObject*)0,Form("rms/mean=%f",th1[0]->GetRMS()/th1[0]->GetMean()));
		leg->AddEntry((TObject*)0,Form("#sigma/mean=%f",tf1[0]->GetParameter(2)/tf1[0]->GetParameter(1)));
		leg->AddEntry(th1[1],"L2L3 pruned");
		
		//leg->AddEntry((TObject*)0,Form("rms/mean=%f",th1[1]->GetRMS()/th1[1]->GetMean()));
		leg->AddEntry((TObject*)0,Form("#sigma/mean=%f",tf1[1]->GetParameter(2)/tf1[1]->GetParameter(1)));
		
		leg->AddEntry(th1[3],"regression with L2L3pruned");
		//leg->AddEntry((TObject*)0,Form("rms/mean=%f",th1[2]->GetRMS()/th1[2]->GetMean()));
		leg->AddEntry((TObject*)0,Form("#sigma/mean=%f",tf1[3]->GetParameter(2)/tf1[3]->GetParameter(1)));
		
		
		th1[1]->SetLineWidth(3);
		th1[2]->SetLineWidth(3);
		
		th1[0]->Draw("");
		th1[1]->Draw("same");
		th1[2]->Draw("same");
		th1[3]->Draw("same");
		leg->Draw("same");
		c1->Print(Form("plots/B%s.pdf",masspoint[i].data()));
	}
	
}
