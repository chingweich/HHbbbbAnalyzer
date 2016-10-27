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

void makeCorr3(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	
	TFile *f;
	f=TFile::Open("corr2/corr.root");
	
	TFile *f2;
	f2=TFile::Open("pdgToReco2.root");
	
	TFile* tf1[9];
	int masspoint[9]={1200,1400,1600,1800,2000,2500,3000,4000,4500};
	for(int i=0;i<9;i++){
		tf1[i]=TFile::Open(Form("corr2/B%d.root",masspoint[i]));
		//f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	double xsec[9]={1.90,0.763,0.33,0.155,7.65e-2,1.58e-2,3.73e-3,2.08e-4,4.91e-5};
		
	double ptBins[14]={300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	double ptBinsCenter[14]={350,450,550,650,750,850,950,1125,1375,1625,1875,2250,2750};
	double ptBinsCenterE[14]={350,450,550,650,750,850,950,1125,1375,1625,1875,2250,2750};
	double ptBinsError[14]={0};
	double ptBinsErrorE[14]={0};
	
	double mean[6][15];
	double sigma[6][15];
	
	
	TLegend *leg = new TLegend(0.68, 0.45, 0.94, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	
	TGraphErrors* tg1[9];
	
	for(int k=0;k<9;k++){
		for(int i=0;i<13;i++){
		TH1D* th1=(TH1D	*)tf1[k]->Get(Form("ptBarel%.0f",ptBins[i]));
		ptBinsCenter[i]=th1->GetMean();
		ptBinsError[i]=th1->GetRMS();
		cout<<i<<","<<ptBinsCenter[i]<<","<<ptBinsError[i]<<endl;
	}
	
	
		
		for(int i=0;i<13;i++){
		TH1D* th1=(TH1D*)tf1[k]->Get(Form("recoBarelMass%.0f",ptBins[i]));
		TF1 *tf1[4];
	if(i<3){
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
		}
		else {
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-20,th1->GetMaximumBin()+20);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-20,th1->GetMaximumBin()+20);
		}
		mean[4][i]=125/tf1[0]->GetParameter(1);
		sigma[4][i]=tf1[0]->GetParError(1)*(125/tf1[0]->GetParameter(1))/tf1[0]->GetParameter(1);
		double temp=tf1[0]->GetParError(1)/tf1[0]->GetParameter(1);
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		
		//cout<<i<<"="<<mean[4][i]<<endl;
		//sigma[4][i]=0;
		th1->SetTitle(Form("%.0f",ptBins[i]));
		if(i==0)c1->Print(Form("plots/recoBarel%d.pdf(",masspoint[k]));
		else if(i==12)c1->Print(Form("plots/recoBarel%d.pdf)",masspoint[k]));
		else  c1->Print(Form("plots/recoBarel%d.pdf",masspoint[k]));
		if(th1->GetEntries()<50|| mean[4][i]<0 || mean[4][i]>2.5||temp>0.3){
			mean[4][i]=0;
			sigma[4][i]=0;
		}
	}
	
	
	
	tg1[k]=new TGraphErrors(12,ptBinsCenter,mean[4],ptBinsError,sigma[4]);
	
	
	
	tg1[k]->GetXaxis()->SetTitle("jet Pt");
	tg1[k]->GetXaxis()->SetRangeUser(0,1800);
	tg1[k]->GetYaxis()->SetTitle("M_{PDG}/M_{Reco}");
	tg1[k]->SetTitle("Gen Correction");
	tg1[k]->SetMinimum(1);
	tg1[k]->SetMaximum(2.6);
	tg1[k]->SetFillColor(0);
	tg1[k]->SetLineColor(97-5*k);
	tg1[k]->SetMarkerColor(97-5*k);
	tg1[k]->SetMarkerStyle(20+k);
	tg1[k]->SetMarkerSize(2);
	tg1[k]->SetLineWidth(3);
		
		
		leg->AddEntry(tg1[k],Form("M=%d",masspoint[k]));
	
	}
	for(int k=0;k<9;k++){
		tg1[k]->GetXaxis()->SetRangeUser(200,2500);
		if(k==0)tg1[k]->Draw("APL");
	else tg1[k]->Draw("samePL");
	
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		cout<<"k="<<k<<",j="<<j<<","<<x<<","<<y<<endl;
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
		tg1[k]->GetPoint(j,x,y);
		cout<<"k="<<k<<",j="<<j<<","<<x<<","<<y<<endl;
	}
	
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
	}
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
	}
	}
	
	
	TGraphErrors* te1=(TGraphErrors*)f2->Get("barel");
	leg->AddEntry(te1,"avg");
	te1->Draw("samePL");
	te1->SetLineColor(1);
	te1->SetLineWidth(3);
	te1->SetMarkerColor(1);
	te1->SetMarkerSize(1.5);
	te1->SetMarkerStyle(20);
	leg->Draw("same");
	c1->Print("plots/recoOneBarel.pdf");
	
	
	leg->Clear();
	
		for(int k=0;k<9;k++){

	
	for(int i=0;i<11;i++){
		TH1D* th1=(TH1D*)tf1[k]->Get(Form("ptEndcap%.0f",ptBins[i]));
		ptBinsCenterE[i]=th1->GetMean();
		ptBinsErrorE[i]=th1->GetRMS();
		cout<<i<<","<<ptBinsCenterE[i]<<","<<ptBinsErrorE[i]<<endl;
	}
		
		
		for(int i=0;i<13;i++){
		TH1D* th1=(TH1D*)tf1[k]->Get(Form("recoEndcapMass%.0f",ptBins[i]));
		TF1 *tf1[4];
	if(i<3){
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-40,th1->GetMaximumBin()+40);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-40,th1->GetMaximumBin()+40);
		}
		else {
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
			tf1[0]->SetParameter(1,100);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
		
		}
		mean[4][i]=125/tf1[0]->GetParameter(1);
		sigma[4][i]=tf1[0]->GetParError(1)*(125/tf1[0]->GetParameter(1))/tf1[0]->GetParameter(1);
		double temp=tf1[0]->GetParError(1)/tf1[0]->GetParameter(1);
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		
		//cout<<i<<"="<<mean[4][i]<<endl;
		//sigma[4][i]=0;
		th1->SetTitle(Form("%.0f",ptBins[i]));
		if(i==0)c1->Print(Form("plots/recoEndcap%d.pdf(",masspoint[k]));
		else if(i==12)c1->Print(Form("plots/recoEndcap%d.pdf)",masspoint[k]));
		else  c1->Print(Form("plots/recoEndcap%d.pdf",masspoint[k]));
		if(th1->GetEntries()<50|| mean[4][i]<1 || mean[4][i]>2.5||temp>0.3){
			mean[4][i]=0;
			sigma[4][i]=0;
		}
	}
	
	
	
	tg1[k]=new TGraphErrors(9,ptBinsCenter,mean[4],ptBinsError,sigma[4]);
	
	
	
	tg1[k]->GetXaxis()->SetTitle("jet Pt");
	tg1[k]->GetXaxis()->SetRangeUser(0,1800);
	tg1[k]->GetYaxis()->SetTitle("M_{PDG}/M_{Reco}");
	tg1[k]->SetTitle("Gen Correction");
	tg1[k]->SetMinimum(1);
	tg1[k]->SetMaximum(1.8);
	tg1[k]->SetFillColor(0);
	tg1[k]->SetLineColor(97-5*k);
	tg1[k]->SetMarkerColor(97-5*k);
	tg1[k]->SetMarkerStyle(20+k);
	tg1[k]->SetMarkerSize(2);
	tg1[k]->SetLineWidth(3);
		
		
		leg->AddEntry(tg1[k],Form("M=%d",masspoint[k]));
	
	}
	for(int k=0;k<9;k++){
		tg1[k]->GetXaxis()->SetRangeUser(200,2500);
		if(k==0)tg1[k]->Draw("APL");
	else tg1[k]->Draw("samePL");
	
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		cout<<"k="<<k<<",j="<<j<<","<<x<<","<<y<<endl;
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
		tg1[k]->GetPoint(j,x,y);
		cout<<"k="<<k<<",j="<<j<<","<<x<<","<<y<<endl;
	}
	
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
	}
	for(int j=0;j<12;j++){
		double x,y;
		tg1[k]->GetPoint(j,x,y);
		
		if(y<0.1){
			tg1[k]->RemovePoint (j);
			//j--;
		}
	}
	}
	 te1=(TGraphErrors*)f2->Get("endcap");
	leg->AddEntry(te1,"avg");
	te1->Draw("samePL");
	te1->SetLineColor(1);
	te1->SetLineWidth(3);
	te1->SetMarkerColor(1);
	te1->SetMarkerSize(1.5);
	te1->SetMarkerStyle(20);
	leg->Draw("same");
	c1->Print("plots/recoOneEndcap.pdf");
	
	/*
	
	
	for(int i=0;i<13;i++){
		TH1D* th1=(TH1D*)f->Get(Form("recoEndcapMass%.0f",ptBins[i]));
		TF1 *tf1[4];
		if(i<2){
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-30,th1->GetMaximumBin()+30);
		}
		else {
			tf1[0]=new TF1("fa1","gaus(25000)",th1->GetMaximumBin()-20,th1->GetMaximumBin()+20);
		th1->Fit(tf1[0],"","",th1->GetMaximumBin()-20,th1->GetMaximumBin()+20);
		}
	
		mean[5][i]=125/tf1[0]->GetParameter(1);
		sigma[5][i]=tf1[0]->GetParError(1)*(125/tf1[0]->GetParameter(1)	)/tf1[0]->GetParameter(1);
		//mean[5][i]=125/th1->GetMean();
		//sigma[5][i]=th1->GetMeanError()/th1->GetMean();
			th1->Draw();
			th1->SetTitle(Form("%.0f",ptBins[i]));
		tf1[0]->Draw("same");
	//cout<<i<<"="<<mean[5][i]<<endl;
	
  
	}
	
	
	
	
	TGraphErrors* tg1[6];
	
	tg1[4]=new TGraphErrors(12,ptBinsCenter,mean[4],ptBinsError,sigma[4]);
	tg1[5]=new TGraphErrors(9,ptBinsCenter,mean[5],ptBinsError,sigma[5]);
	
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[4][i]<<endl;
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[5][i]<<endl;
	

	
	
  

	tg1[4]->Draw("APL");
	tg1[4]->SetFillColor(0);
	tg1[5]->SetFillColor(0);
	tg1[5]->SetLineColor(2);
	tg1[5]->SetMarkerColor(2);
	tg1[5]->Draw("PLsame");
	

leg->Clear();


  
  leg->AddEntry(tg1[5],"reco endcap");
  leg->AddEntry(puppisd_corrRECO_cen,"Thea barel");
  leg->AddEntry(puppisd_corrRECO_for,"Thea endcap");

	leg->Draw("same");
	genBarel->Draw("same");
	genEndcap->Draw("same");
	puppisd_corrRECO_cen->SetLineColor(3);
	puppisd_corrRECO_for->SetLineColor(4);
	//tg1[2]->Draw("APL");
	puppisd_corrRECO_cen->Draw("same");
	puppisd_corrRECO_for->Draw("same");
	
	c1->Print("plots/recoOne.pdf");
	*/
}
