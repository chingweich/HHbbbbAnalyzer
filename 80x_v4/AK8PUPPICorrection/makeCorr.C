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

void makeCorr(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	
	TFile *f;
	f=TFile::Open("corr/corr.root");
		
	double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	double ptBinsCenter[14]={250,350,450,550,650,750,850,950,1125,1375,1625,1875,2250,2750};
	double ptBinsError[14]={0};
	
	double mean[4][15];
	double sigma[4][15];
	
	for(int i=0;i<14;i++){
		TH1D* th1=(TH1D*)f->Get(Form("genBarelMass%.0f",ptBins[i]));
		TF1 *tf1[4];
	tf1[0]=new TF1("fa1","gaus(25000)",50,150);
		th1->Fit(tf1[0],"","",50,150);
		mean[0][i]=125/tf1[0]->GetParameter(1);
		sigma[0][i]=tf1[0]->GetParameter(2)/tf1[0]->GetParameter(1);
		th1->Draw();
		th1->SetTitle(Form("%.0f",ptBins[i]));
		tf1[0]->Draw("same");
		if(i==0)c1->Print("plots/genBarel.pdf(");
		else if(i==13)c1->Print("plots/genBarel.pdf)");
		else  c1->Print("plots/genBarel.pdf");
	}
	
	for(int i=0;i<14;i++){
		TH1D* th1=(TH1D*)f->Get(Form("genEndcapMass%.0f",ptBins[i]));
		TF1 *tf1[4];
	tf1[0]=new TF1("fa1","gaus(25000)",50,150);
		th1->Fit(tf1[0],"","",50,150);
		mean[1][i]=125/tf1[0]->GetParameter(1);
		sigma[1][i]=tf1[0]->GetParameter(2)/tf1[0]->GetParameter(1);
			th1->Draw();
			th1->SetTitle(Form("%.0f",ptBins[i]));
		tf1[0]->Draw("same");
		if(i==0)c1->Print("plots/genEndcap.pdf(");
		else if(i==13)c1->Print("plots/genEndcap.pdf)");
		else  c1->Print("plots/genEndcap.pdf");
	}
	
	for(int i=0;i<14;i++){
		TH1D* th1=(TH1D*)f->Get(Form("recoBarelMass%.0f",ptBins[i]));
		TF1 *tf1[4];
		tf1[0]=new TF1("fa1","gaus(25000)",50,150);
		th1->Fit(tf1[0],"","",50,150);
		mean[2][i]=125/(tf1[0]->GetParameter(1)*mean[0][i]);
		sigma[2][i]=sqrt(pow((tf1[0]->GetParameter(2)/tf1[0]->GetParameter(1)),2)+pow(sigma[0][i],2));
			th1->Draw();
			th1->SetTitle(Form("%.0f",ptBins[i]));
		tf1[0]->Draw("same");
		if(i==0)c1->Print("plots/recoBarel.pdf(");
		else if(i==13)c1->Print("plots/recoBarel.pdf)");
		else  c1->Print("plots/recoBarel.pdf");
	}
	
	for(int i=0;i<14;i++){
		TH1D* th1=(TH1D*)f->Get(Form("recoEndcapMass%.0f",ptBins[i]));
		TF1 *tf1[4];
		tf1[0]=new TF1("fa1","gaus(25000)",50,150);
		th1->Fit(tf1[0],"","",50,150);
		mean[3][i]=125/(tf1[0]->GetParameter(1)*mean[1][i]);
		sigma[3][i]=sqrt(pow((tf1[0]->GetParameter(2)/tf1[0]->GetParameter(1)),2)+pow(sigma[1][i],2));
			th1->Draw();
			th1->SetTitle(Form("%.0f",ptBins[i]));
		tf1[0]->Draw("same");
		if(i==0)c1->Print("plots/recoEndcap.pdf(");
		else if(i==13)c1->Print("plots/recoEndcap.pdf)");
		else  c1->Print("plots/recoEndcap.pdf");
	}
	
	
	TGraphErrors* tg1[4];
	tg1[0]=new TGraphErrors(14,ptBinsCenter,mean[0],ptBinsError,sigma[0]);
	tg1[1]=new TGraphErrors(11,ptBinsCenter,mean[1],ptBinsError,sigma[1]);
	tg1[2]=new TGraphErrors(14,ptBinsCenter,mean[2],ptBinsError,sigma[2]);
	tg1[3]=new TGraphErrors(11,ptBinsCenter,mean[3],ptBinsError,sigma[3]);
	
	tg1[0]->GetXaxis()->SetTitle("jet Pt");
	tg1[0]->GetYaxis()->SetTitle("M_{PDG}/M_{Gen}");
	tg1[0]->SetTitle("Gen Correction");
	tg1[0]->Draw("APL");
	tg1[0]->SetFillColor(0);
	tg1[1]->SetFillColor(0);
	tg1[1]->SetLineColor(2);
	tg1[1]->SetMarkerColor(2);
	tg1[1]->Draw("PLsame");
	TF1* genBarel = new TF1("genBarel","[0]+[1]*pow(x*[2],-[3])");
	  genBarel->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
	TF1* genEndcap = new TF1("genEndcap","[0]+[1]*pow(x*[2],-[3])");
	  genEndcap->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
	genBarel->SetLineColor(1);
	genEndcap->SetLineColor(2);
	tg1[0]->Fit(genBarel);
	tg1[1]->Fit(genEndcap);
	
	TLegend *leg = new TLegend(0.68, 0.65, 0.94, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  
  leg->AddEntry(tg1[0],"gen barel");
  leg->AddEntry(tg1[1],"gen endcap");

	leg->Draw("same");
	genBarel->Draw("same");
	genEndcap->Draw("same");
	
	c1->Print("plots/gen.pdf");
	
	tg1[2]->GetXaxis()->SetTitle("jet Pt");
	tg1[2]->GetYaxis()->SetTitle("M_{Gen}/M_{Reco}");
	tg1[2]->SetTitle("Reco Correction");
	tg1[2]->Draw("APL");
	tg1[2]->SetFillColor(0);
	tg1[3]->SetFillColor(0);
	tg1[3]->SetLineColor(2);
	tg1[3]->SetMarkerColor(2);
	tg1[3]->Draw("PLsame");
	
	  TF1* recoBarel = new TF1("recoBarel","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  recoBarel->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18
   				      );
		  TF1* recoEndcap= new TF1("recoEndcap","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  recoEndcap->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18
   				      );

	recoBarel->SetLineColor(1);
	recoEndcap->SetLineColor(2);
	tg1[2]->Fit(recoBarel);
	tg1[3]->Fit(recoEndcap);
	
	leg->Clear();
	 leg->AddEntry(tg1[0],"reco barel");
  leg->AddEntry(tg1[1],"reco endcap");

	leg->Draw("same");
	recoBarel->Draw("same");
	recoEndcap->Draw("same");
	c1->Print("plots/reco.pdf");
	
}
