#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
//#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "../setNCUStyle.C"
#include<TH2.h>
#include "TLine.h"
#include "TF1.h"

#define nMassp 10
#define nSt 11

TCanvas* c1;
void uncert_btag(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[11],* tf2[11];
	TFile* tf3[11],* tf4[11];
	TFile* tf5[11],* tf6[11];
	string st[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	for (int i=0;i<11;i++)tf1[i]=TFile::Open(Form("../sf/B%s.root",st[i].data()));
	for (int i=0;i<11;i++)tf3[i]=TFile::Open(Form("../sf/B%s_BtagUp.root",st[i].data()));
	for (int i=0;i<11;i++)tf5[i]=TFile::Open(Form("../sf/B%s_BtagDown.root",st[i].data()));
	for (int i=0;i<11;i++)tf2[i]=TFile::Open(Form("../sf/R%s.root",st[i].data()));
	double x[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	double xx[11]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
	double y[6][11];

	for (int i=0;i<11;i++){
		TH1D * th1= (TH1D *)tf1[i]->FindObjectAny("totalMass");
		y[0][i]=th1->GetMean();
		th1=(TH1D *)tf3[i]->FindObjectAny("totalMass");
		y[1][i]=th1->GetMean();
		th1=(TH1D *)tf5[i]->FindObjectAny("totalMass");
		y[2][i]=th1->GetMean();
		y[1][i]/=y[0][i];
		y[2][i]/=y[0][i];
		y[0][i]/=y[0][i];
	}
	
	double yy[6][10];
	for (int i=0;i<11;i++){
		if(i>=8){
			yy[0][i]=y[0][i+1];
			yy[1][i]=y[1][i+1];
			yy[2][i]=y[2][i+1];
		}
		else {
			yy[0][i]=y[0][i];
			yy[1][i]=y[1][i];
			yy[2][i]=y[2][i];
		}
	}
	
	//TGraph * tg1=new TGraph(11,x,y[0]);
	//TGraph * tg2=new TGraph(11,x,y[1]);
	//TGraph * tg3=new TGraph(11,x,y[2]);
	TGraph * tg1=new TGraph(10,xx,yy[0]);
	TGraph * tg2=new TGraph(10,xx,yy[1]);
	TGraph * tg3=new TGraph(10,xx,yy[2]);
	
	tg1->SetTitle("Btag uncertainty");
	tg2->SetTitle("Btag_up fit");
	tg3->SetTitle("Btag_down fit");
	tg1->SetLineColor(1);
	tg2->SetLineColor(2);
	tg3->SetLineColor(3);
	
	tg1->SetFillColor(kWhite);
	tg2->SetFillColor(kWhite);
	tg3->SetFillColor(kWhite);
	
	tg1->GetXaxis()->SetTitle("M_{X}[GeV]");
	tg1->GetYaxis()->SetTitle("Mean/Mean_{central}");	
	tg1->SetMaximum(1.001);
	tg1->SetMinimum(0.999);
	tg1->Draw("APL");
	tg2->Draw("PL,same");
	tg3->Draw("PL,same");
	
	TLegend* leg ;
	leg=new TLegend(0.391452,0.712447,0.590645,0.863966);
	leg->AddEntry(tg1,"Btag_central");
	leg->AddEntry(tg2,"Btag_up");
	leg->AddEntry(tg3,"Btag_down");
	
	leg->Draw("same");
	
	c1->Print("plot/all.pdf");
	
	TF1 *fa1 = new TF1("fa1","[0]+[1]*x",1000,4500); 
	TF1 *fa2 = new TF1("fa1","[0]+[1]*x",1000,4500); 
	tg2->Fit(fa1);
	fa1->SetLineColor(4);
	
	tg2->Draw("APL");
	fa1->Draw("same");
	c1->Print("plot/fit_up.pdf");
	
	tg3->Fit(fa2);
	fa2->SetLineColor(4);
	tg3->Draw("APL");
	fa2->Draw("same");
	c1->Print("plot/fit_down.pdf");

}

