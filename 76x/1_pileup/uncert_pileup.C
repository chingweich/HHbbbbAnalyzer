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
void uncert_pileup(){
	setNCUStyle();
	c1 = new TCanvas("c1","",1360,768);
	
	TFile* tf1[11],* tf2[11];
	string st[11]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
	
	tf1[i]=TFile::Open(Form("../sf/B%s.root",st[i].data()));
	tf2[i]=TFile::Open(Form("../sf/R%s.root",st[i].data()));
	double x[11]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	double y[6][11];

	for (int i=0;i<11;i++){
		TH1D * th1=tf1[i]->FindObjectAny("totalMass");
		y[0][i]=th1->GetMean();
		th1=tf1[i]->FindObjectAny("totalMass_pileup_up");
		y[1][i]=th1->GetMean();
		th1=tf1[i]->FindObjectAny("totalMass_pileup_down");
		y[2][i]=th1->GetMean();
	}

	
	TGraph * tg1=new TGraph(11,x,y[0]);
	TGraph * tg2=new TGraph(11,x,y[1]);
	TGraph * tg3=new TGraph(11,x,y[2]);
	//TGraph * tg1=new TGraph(11,x,y1);
	//TGraph * tg2=new TGraph(11,x,y2);
	//TGraph * tg3=new TGraph(11,x,y3);
	
	
	tg1->SetLineColor(1);
	tg2->SetLineColor(2);
	tg3->SetLineColor(3);
	
	tg1->SetFillColor(kWhite);
	tg2->SetFillColor(kWhite);
	tg3->SetFillColor(kWhite);
	
	tg1->GetXaxis()->SetTitle("M_{X}[GeV]");
	tg1->GetYaxis()->SetTitle("Mean/M_{X}[GeV](normalized)");	
	tg1->SetMaximum(1.05);
	tg1->Draw("APL");
	tg2->Draw("PL,same");
	tg3->Draw("PL,same");
	
	TLegend* leg ;
	leg=new TLegend(0.191452,0.562447,0.290645,0.783966);
	leg->AddEntry(tg1,"JES");
	leg->AddEntry(tg2,"JES_UP");
	leg->AddEntry(tg3,"JES_DO");
	
	leg->Draw("same");
	
	c1->Print("plot/all.pdf");
	
	TF1 *fa1 = new TF1("fa1","[0]+[1]*x",1200,4500); 
	TF1 *fa2 = new TF1("fa1","[0]+[1]*x",1200,4500); 
	tg2->Fit(fa1);
	tg2->Draw("APL");
	fa1->Draw("same");
	c1->Print("plot/fit_up.pdf");
	
	tg3->Fit(fa2);
	tg3->Draw("APL");
	fa2->Draw("same");
	c1->Print("plot/fit_down.pdf");

}

